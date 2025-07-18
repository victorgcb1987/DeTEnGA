#!/usr/bin/env python

import argparse
import os
import sys

from pathlib import Path


from src.parsers import (parse_fof, get_pfams_from_db, get_pfams_from_interpro_query, 
                         parse_TEsort_output, classify_pfams, create_summary, write_summary,
                         get_stats)
from src.run import run_gffread, run_TEsorter, remove_stop_codons, run_interpro, run_agat

REXDB_PFAMS = {"rexdb-plant": Path(os.path.dirname(os.path.realpath(__file__))) / "data" / "Viridiplantae_2.0_pfams.txt",
               "rexdb-metazoa": Path(os.path.dirname(os.path.realpath(__file__))) / "data" / "Metazoa_3.1_pfams.txt",
               "rexdb": Path(os.path.dirname(os.path.realpath(__file__))) / "data" / "Combined_pfams.txt"}


CATEGORIES = {"No_TE(PcpM0)": "PcpM0", "Protein_TE_only(PteM0)": "PteM0",
              "Chimeric_Protein_Only(PchM0)": "PchM0", "mRNA_TE_Only(PcpMte)": "PcpMte",
              "Protein_and_mRNA_TE(PteMte)": "PteMte", "Chimeric_Protein_and_mRNA_TE(PchMte)": "PchMte",
              "No_Protein_Domains_mRNA_TE(P0Mte)": "P0Mte"}

#Generating program options
def parse_arguments():
    desc = "Pipeline to identify Transposable Elments (TE) in annotated genes"
    parser = argparse.ArgumentParser(description=desc)
    
    
    help_input = '''(Required) File of Files with the following format:
                    "NAME   FASTA   GFF'''
    parser.add_argument("--input", "-i", type=str,
                        help=help_input,
                        required=True)
    
    help_output_dir = '''(Required) Output dir'''
    parser.add_argument("--output", "-out", type=str,
                        help=help_output_dir,
                        required=True)
    
    help_threads = "(Optional) number of threads. 1 by default"
    parser.add_argument("--threads", "-t", type=int,
                        help=help_threads, default=1)
    
    help_database = "(Optional) database for TEsorter. rexdb-plant by default"
    parser.add_argument("--tesorter_database", "-d", type=str,
                        help=help_database, default="rexdb-plant")
    
    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()

def get_arguments():
    parser = parse_arguments()
    output = Path(parser.output)
    if not output.exists():
        output.mkdir(parents=True)
    return {"input": parser.input,
            "out": output,
            "threads": parser.threads,
            "tesorter_database": parser.tesorter_database}


def create_header():
    header = ["Run", "Genome", "Annotation", "Annotated_transcripts"]
    for key in CATEGORIES:
        header.append(f"{key}_N")
    for key in CATEGORIES:
        header.append(f"{key}_%")
    header += ["Summary_N", "Summary_%"]
    return "\t".join(header)+"\n"


def get_row(label, genome, annotation, stats):
    inverse_categories = {value: key for key, value in CATEGORIES.items()}
    categories = ["T"] + [key for key in inverse_categories]
    values = [str(stats["num_transcripts"])] + [str(stats[key]) for key in inverse_categories]
    per_values = [str(stats["num_transcripts"])] + [str(round(float(stats[key]/stats["num_transcripts"])*100, 2)) for key in inverse_categories]
    summary = "{0}: {1};{2}: {3};{4}: {5};{6}: {7};{8}: {9};{10}: {11};{12}: {13}"
    row = [label, genome, annotation]
    row += values
    row += per_values[1:]
    row += [summary.format(*[item for pair in zip(categories, values) for item in pair])]
    row += [summary.format(*[item for pair in zip(categories, per_values) for item in pair])]    
    return "\t".join(row)+"\n"


def main():
    args = get_arguments()
    files = parse_fof(args["input"])
    out_dir = args["out"]
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)
    #Create log file
    log = out_dir / "log.txt"
    log_fhand = open(log, "a")
    log_fhand.write("#Command used: {}\n".format(" ".join(sys.argv)))
    

    #Retrieve sequences
    msg = "##STEP 1: Retrive sequences with gffread\n"
    print(msg)
    log_fhand.write(msg)
    sequences = run_gffread(files, args["out"])
    failed_runs = []
    for label, values in sequences.items():
        for kind in ["mrna", "protein"]:
            log_fhand.write("{} | {}\n".format(values["command"][kind], values["msg"][kind]))
            log_fhand.flush()
        if values["returncode"]["mrna"] == 1 or values["returncode"]["protein"] == 1:
            failed_runs.append(label)
            log_fhand.write("Removed {} from pipeline, please check the error message\n\n".format(label))
            log_fhand.flush()
    for key in failed_runs:
        sequences.pop(key)

    #Create TEsorter input
    msg = "##STEP 2: Analyze mRNA transposable elements with TEsorter\n"
    print(msg)
    log_fhand.write(msg)
    log_fhand.flush()
    TEsorter_results = run_TEsorter(sequences, args["tesorter_database"], args["threads"])
    failed_runs = []
    for label, values in TEsorter_results.items():
        log_fhand.write("{} | {}\n".format(values["command"], values["msg"]))
        log_fhand.flush()
        if values["returncode"] == 1:
            failed_runs.append(label)
            log_fhand.write("Removed {} from pipeline, please check the error message\n\n".format(label))
            log_fhand.flush()

    for key in failed_runs:
        sequences.pop(key)
        TEsorter_results.pop(key)
    #Trim Sequences with internal stop codons
    msg = "##STEP 3: Remove internal stop codons from proteins\n"
    print(msg)
    log_fhand.write(msg)
    log_fhand.flush()
    no_stop_codons_sequences = {}
    for key, values in sequences.items():
        protein_seqs = values["out_fpath"]["protein"]
        results = remove_stop_codons(protein_seqs)
        log_fhand.write("{} | {}\n".format(results["command"], results["msg"]))
        no_stop_codons_sequences[key] = results
    

    #Run interproscan
    msg = "##STEP 4: Analyze protein transposable elements with interproscan\n"
    print(msg)
    log_fhand.write(msg)
    log_fhand.flush()
    interpro_results = run_interpro(no_stop_codons_sequences, args["threads"])
    failed_runs = []
    for label, values in interpro_results.items():
         if values["returncode"] == 1:
            failed_runs.append(label)
            log_fhand.write("Removed {} from pipeline, please check the error message\n\n".format(label))
         else:
            log_fhand.write("{} | {}\n".format(values["command"], values["msg"]))
    for label in failed_runs:
        interpro_results.pop(label)


    msg = "##STEP 5: merging evidences from interpro and TEsorter\n"
    print(msg)
    log_fhand.write(msg)
    log_fhand.flush()
    database = REXDB_PFAMS[args["tesorter_database"]]
    TE_pfams = get_pfams_from_db(database)
    summaries = {}
    for label in sequences:
        if label in interpro_results and label in TEsorter_results:
            with open(TEsorter_results[label]["out_fpath"]) as TEsorter_fhand:
                te_sorter_output = parse_TEsort_output(TEsorter_fhand)
        
            with open(interpro_results[label]["out_fpath"]) as interpro_fhand:
                interpro = get_pfams_from_interpro_query(interpro_fhand)
                classified_pfams = classify_pfams(interpro, TE_pfams)
    
            te_summary = create_summary(classified_pfams, te_sorter_output)
    
            out_fpath = Path(out_dir / label / "{}_TE_summary.csv".format(label))
            with open(out_fpath, "w") as out_fhand:
                write_summary(te_summary, out_fhand)
                summaries[label] = out_fpath
                msg = "TE Summary for {} written in {}\n".format(label, out_fpath)
                log_fhand.write(msg)
                log_fhand.flush()
    
    msg = "##STEP 6: Running stats on annotation files\n"
    print(msg)
    log_fhand.write(msg)
    log_fhand.flush()
    agat_results = run_agat(summaries, files)
    with open(args["out"]/ "combined_summaries.tsv", "w") as combined_summaries_fhand:
        header = create_header()
        combined_summaries_fhand.write(header)
        for label, results in agat_results.items():
            stats = get_stats(results["out_fpath"], summaries[label])
            genome = files[label]["assembly"].stem
            annotation = files[label]["annotation"].stem
            row = get_row(label, genome, annotation, stats)
            combined_summaries_fhand.write(row)


if __name__ == "__main__":
    main()