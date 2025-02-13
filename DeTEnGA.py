#!/usr/bin/env python

import argparse
import os
import sys

from pathlib import Path


from src.parsers import (parse_fof, get_pfams_from_db, get_pfams_from_interpro_query, 
                         parse_TEsort_output, classify_pfams, create_summary, write_summary,
                         get_stats)
from src.run import run_gffread, run_TEsorter, remove_stop_codons, run_interpro, run_agat

REXDB_PFAMS = Path(os.path.dirname(os.path.realpath(__file__))) / "data" / "rexdb_Viridiplantae_4.0_pfams.txt"


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
    parser.add_argument("--tesorter_database", "-s", type=str,
                        help=help_database, default="rexdb-plant")
    
    help_combine = "(Optional) combine all summaries. False by default"
    parser.add_argument("--combine", "-c", action="store_true",
                        help=help_combine)
    
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
            "tesorter_database": parser.tesorter_database,
            "combine": parser.combine}


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
    print(failed_runs)
    for label in failed_runs:
        interpro_results.pop(label)


    msg = "##STEP 5: merging evidences from interpro and TEsorter\n"
    print(msg)
    log_fhand.write(msg)
    log_fhand.flush()
    TE_pfams = get_pfams_from_db(REXDB_PFAMS)
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
    
    if args["combine"]:
        msg = "##STEP 6: Running stats on annotation files\n"
        print(msg)
        log_fhand.write(msg)
        log_fhand.flush()
        agat_results = run_agat(summaries, files)
        with open(args["out"]/ "combined_summaries.tsv", "w") as combined_summaries_fhand:
            header = "Run\tAnnotated_transcripts(N)\tCoding_proteins(N)"
            header += "\tTE_proteins(N)\tMixed_Proteins(N)"
            header += "\tTE_mRNA(N)\tNonTE_mRNA(N)"
            header += "\tTE_detected_in_both(N)"
            header += "\tCoding_proteins(%)"
            header += "\tTE_proteins(%)\tMixed_Proteins(%)"
            header += "\tTE_mRNA(%)\tNonTE_mRNA(%)"
            header += "\tTE_detected_in_both(%)\n"
            combined_summaries_fhand.write(header)
            for label, results in agat_results.items():
                stats = get_stats(results["out_fpath"], summaries[label])
                line = f'{label}\t{stats["num_transcripts"]}'
                line += f'\t{stats["coding_protein"]}\t{stats["te_protein"]}'
                line += f'\t{stats["mixed_protein"]}\t{stats["te_mrna"]}'
                line += f'\t{stats["nonte_mrna"]}\t{stats["both"]}'
                line += f'\t{float(stats["coding_protein"]/stats["num_transcripts"])}'
                line += f'\t{float(stats["te_protein"]/stats["num_transcripts"])}'
                line += f'\t{float(stats["mixed_protein"]/stats["num_transcripts"])}'
                line += f'\t{float(stats["te_mrna"]/stats["num_transcripts"])}'
                line += f'\t{float(stats["nonte_mrna"]/stats["num_transcripts"])}'
                line += f'\t{float(stats["both"]/stats["num_transcripts"])}'
                line += "\n"
                combined_summaries_fhand.write(line)

    print("Program finished. Exiting")

if __name__ == "__main__":
    main()