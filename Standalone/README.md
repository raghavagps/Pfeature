# Standalone Package of Pfeature
## Introduction
Pfeature is developed for computing wide range of protein and peptides features from their amino acid sequences, and structures. More information on Pfeature is available from its web server https://webs.iiitd.edu.in/raghava/pfeature. This page provide information about standalone version of Pfeature. This standalone contains three scripts, their description is as follows:
  - pfeature_comp.py
  - pfeature_bin.py
  - pfeature_pssm.py
  
  - Minimum USAGE
```sh
Minimum ussage is "pfeature_comp.py -i protein.fa" where protein.fa is a input fasta file. This will calculate the amino acid composition of the seqeunces provided in the fasta file. It will use other parameters by default. It will save output in "pfeature_result.csv" in CSV (comma seperated variables).
 ```

## Pfeature Package Files
It contantain following files, brief description of these files are given below:

* `LICENSE`                  : License information

* `README.txt`               : This file provide information about this package

* `pfeature_comp.py`         : Python program to calculate composition based features

* `pfeature_bin.py`          : Python program to calculate binary profile based features

* `pfeature_pssm.py`         : Python program to calculate pssm profile based features

* `protein.seq`              : Example file contain protein sequences in simple format

* `protein.fa`               : Example file contain protein sequences in FASTA format

* `Data `                    : This folder contains the files required to calcuate the composition and binary profile based features.

* `envfile`                  : This file contains the path information required to run the pfeature_pssm.py script.

* `Pfeature_Descriptors.pdf` : This file comprises of description of the header of output files, which is generated using the aforementioned scripts.

* `Requirement.txt`          : This file consists of commands and pre-requisite to run the aforementioned scripts.
