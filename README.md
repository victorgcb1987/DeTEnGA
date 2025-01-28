# DeTEnGA

Detection of Transposable Elements (TEs) on Gene Annotations (DeTEnGA) is a a tool to detect TEs on sequences annotated as gene ORFs. It uses the outputs of TEsorter for detect TEs in mRNA sequences and interproscan to detect domains associated with TEs in proteins and merge results in a summary file.

## Requirements

**TEsorter** : we recommend installing TEsorter using conda as follows:
`conda create -n TEsorter -c bioconda tesorter`

Then, update python from this conda installation (python v3.6 uses some python deprecated functions):
`conda install python=3.12`

**Interproscan**: we recommend using the github version instead of any conda installation (https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html)
Then, add interproscan.sh to your PATH variable:

`export PATH=$PATH:/path/to/interproscan.sh`



**smudgeplot v0.2.5** : We recommend to use conda (`conda install smudgeplot==0.2.5`) to install this module.  

kmc https://github.com/refresh-bio/KMC. You will need to add kmc binaries to PATH variable for example using `export PATH=$PATH:kmc/install/bin/`

clone KATULU respository: `git clone https://github.com/victorgcb1987/KATULU.git`

conda create -n TEsorter -c bioconda tesorter
conda activate TEsorter
conda install python=3.12

export PATH=$PATH:~/soft/interproscan-5.72-103.0/

fof.txt

label path_to_fasta  path_to_annotation
