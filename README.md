# DeTEnGA

Detection of Transposable Elements (TEs) on Gene Annotations (**DeTEnGA**) is a a tool to detect TEs in sequences annotated as gene ORFs. It uses the outputs of TEsorter for detect TEs in mRNA sequences and interproscan to detect domains associated with TEs in proteins and merge results in a summary file.

## Requirements

**TEsorter**: we recommend installing TEsorter using conda as follows: `conda create -n TEsorter -c bioconda tesorter`. Then, update python from this conda installation (python v3.6 uses some python deprecated functions): `conda install python=3.12`.

**Interproscan**: we recommend using the github version instead of any conda installation (https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html)
Then, add interproscan.sh to your PATH variable:

`export PATH=$PATH:/path/to/interproscan.sh`

**DeTEnGA**: just clone this respository `git clone https://github.com/victorgcb1987/DeTEnGA.git`


## How to use
fof.txt

label path_to_fasta  path_to_annotation
