import os
import re

from collections import defaultdict
from csv import DictReader
from pathlib import Path



def parse_fof(input):
    fof = {}
    with open(input) as fhand:
        for line in fhand:
            label, fasta, gff = line.rstrip().split()
            fof[label] = {"assembly": Path(fasta), 
                          "annotation": Path(gff)}
    return fof

def get_pfams_from_db(fpath):
    pfams = {}
    with open(fpath) as fhand:
        for line in fhand:
            if line.startswith("#"):
                continue
            line = line.split()
            pfams[line[0]] = " ".join(line[1:])
    return pfams

def get_pfams_from_interpro_query(fhand):
    pattern = r"^(\S+)\t(\S+)\t(\S+)\tPfam\t([^\t]+)\t([^\t]+)\t(\d+)\t(\d+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)"

    text = fhand.read()
    genes = defaultdict(list)
    matches = re.findall(pattern, text, re.MULTILINE)
    for match in matches:
        gen, code, description, start, end = match[0], match[3], match[4], match[5], match[6]
        genes[gen].append([code, description, start, end])
    sorted_genes = {
                    key: sorted(value, key=lambda x: int(x[2]))  # Ordena por el tercer valor (convertido a entero)
                    for key, value in genes.items()}
    return sorted_genes

def classify_pfams(interpro, te_pfams):
    for gene, pfams in interpro.items():
        for pfam in pfams:
            if pfam[0] in te_pfams:
                pfam.append("TE")
            else:
                pfam.append("NT")
    return interpro

def parse_TEsort_output(fhand):
    output = defaultdict(list)
    for line in DictReader(fhand,delimiter="\t"):
        if line["Domains"] != "none":
            output[line["#TE"]] = {"domains": line["Domains"], 
                                   "complete": line["Complete"],
                                   "classification": "{}|{}|{}".format(line["Order"],
                                                                       line["Superfamily"],
                                                                       line["Clade"]),
                                    "strand": line["Strand"]}
    return output


def create_summary(interpro_classified, tesort_output):
    summary = []
    for transcript, values in interpro_classified.items():
        row = {"transcript": transcript}
        transposable = False
        no_transposable = False
        pfams_ids = []
        pfams_descriptions = []
        for value in values:
            pfams_ids.append(value[0])
            pfams_descriptions.append(value[1])
            if "TE" in value:
                transposable = True
            if "NT" in value:
                no_transposable = True
        if transposable and not no_transposable:
            status = "transposable_element"
        if not transposable and no_transposable:
            status = "coding_sequence"
        if transposable and no_transposable:
            status = "mixed"
        row["interpro_status"] = status
        row["pfams_ids"] = "|".join(pfams_ids)
        row["pfams_descriptions"] = "|".join(pfams_descriptions)
        transcript_tesort = tesort_output.get(transcript, None)
        if transcript_tesort is not None:
            row["tesort_domains"] = transcript_tesort["domains"]
            row["tesort_complete"] = transcript_tesort["complete"]
            row["tesort_class"] = transcript_tesort["classification"]
            row["tesort_strand"] = transcript_tesort["strand"]
        else:
            row["tesort_domains"] = "NA"
            row["tesort_complete"] = "NA"
            row["tesort_class"] = "NA"
            row["tesort_strand"] = "NA"
        summary.append(row)
    for transcript, values in tesort_output.items():
        if transcript not in interpro_classified:
            row = {"transcript": transcript,
                   "interpro_status": "NA",
                   "pfams_ids": "NA",
                   "pfams_descriptions": "NA",
                   "tesort_domains": values["domains"],
                   "tesort_complete": values["complete"],
                   "tesort_class": values["classification"],
                   "tesort_strand": values["strand"]} 
            summary.append(row)
    return summary
        

def write_summary(summary, out_fhand):
    out_fhand.write("Transcript_ID;Interpro_status;TEsort_class;PFAM_domains;")
    out_fhand.write("PFAM_descriptions;TEsort_domains;TEsort_completness;")
    out_fhand.write("TEsort_strand\n")
    out_fhand.flush()
    for row in summary:
        line_total = "" 
        line_total += "{};{};{};{};".format(row["transcript"], row["interpro_status"],
                                            row["tesort_class"], row["pfams_ids"])
        line_total += "{};{};{};{}\n".format(row["pfams_descriptions"], row["tesort_domains"],
                                             row["tesort_complete"], row["tesort_strand"])
        out_fhand.write(line_total)
        out_fhand.flush()
