from pathlib import Path

def parse_fof(input):
    fof = {}
    with open(input) as fhand:
        for line in fhand:
            label, fasta, gff = line.rstrip().split()
            fof[label] = {"assembly": Path(fasta), 
                          "annotation": Path(gff)}
    return fof