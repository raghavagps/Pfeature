# Pfeature
# Introduction
Pfeature is developed for computing wide range of protein and peptides features from their amino acid sequences. More information on Pfeature is abvailble from its web server https://webs.iiitd.edu.in/raghava/pfeature. This page provide information about standalone version of Pfeature. This standalone contains three scripts, their description is as follows:

#############################################################################################################################################################################################################################

1: Standalone for calculating composition based features:

**Important: To run this script 'Data' folder should be in the same directory.**

Minimum USAGE: Minimum ussage is "pfeature_comp.py -i protein.fa" where protein.fa is a input fasta file. This will calculate the amino acid composition of the seqeunces provided in the fasta file. It will use other parameters by default. It will save output in "pfeature_result.csv" in CSV (comma seperated variables).

#Full Usage: Following is complete list of all options, you may get these options by "pfeature_comp.py -h"
usage: pfeature_comp.py [-h] -i INPUT [-o OUTPUT]
                        [-j {AAC,DPC,TPC,ATC,BTC,PCP,AAI,RRI,PRI,DDR,SEP,SER,SPC,ACR,CTC,CeTD,PAAC,APAAC,QSO,SOC,ALLCOMP}]
                        [-n N_TERMINAL] [-c C_TERMINAL] [-nct NC_TERMINAL]
                        [-rn REST_N] [-rc REST_C] [-s SPLIT] [-d LAG]
                        [-w WEIGHT] [-t PWEIGHT]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: protein or protein.sequence in FASTA format or single sequence per line in single letter code
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default pfeature_result.csv
  -j {AAC,DPC,TPC,ATC,BTC,PCP,AAI,RRI,PRI,DDR,SEP,SER,SPC,ACR,CTC,CeTD,PAAC,APAAC,QSO,SOC,ALLCOMP}, --job {AAC,DPC,TPC,ATC,BTC,PCP,AAI,RRI,PRI,DDR,SEP,SER,SPC,ACR,CTC,CeTD,PAAC,APAAC,QSO,SOC,ALLCOMP}
                        Job Type:
                        AAC: Amino acid composition
                        DPC: Dipeptide composition
                        TPC: Tripeptide composition
                        ATC: Atomic composition
                        BTC: Bond composition
                        PCP: Physico-chemical properties composition
                        AAI: Amino-acid indices composition
                        RRI: Residue repeat information
                        PRI: Physico-chemical properties repeat information
                        DDR: Distance distribution of residues
                        SEP: Shannon entropy of protein
                        SER: Shannon entropy of residues
                        SPC: Shannon entropy of physico-chemical properties
                        ACR: Autocorrelation descriptors
                        CTC: Conjoint triad descriptors
                        CeTD: Composition enhanced transition distribution
                        PAAC: Pseudo amino acid composition
                        APAAC: Amphiphilic pseudo amino acid composition
                        QSO: Quasi sequence order
                        SOC: Sequence order coupling number
                        ALLCOMP:All composition features together except ACR and AAI
                        by default AAC
  -n N_TERMINAL, --n_terminal N_TERMINAL
                        Window Length from N-terminal: by default 0
  -c C_TERMINAL, --c_terminal C_TERMINAL
                        Window Length from C-terminal: by default 0
  -nct NC_TERMINAL, --nc_terminal NC_TERMINAL
                        Residues from N- and C-terminal: by default 0
  -rn REST_N, --rest_n REST_N
                        Number of residues removed from N-terminal, by default 0
  -rc REST_C, --rest_c REST_C
                        Number of residues removed from C-terminal, by default 0
  -s SPLIT, --split SPLIT
                        Number of splits a sequence divided into, by default 0
  -d LAG, --lag LAG     This represents the order of gap, lag or dipeptide, by default 1
  -w WEIGHT, --weight WEIGHT
                        Weighting Factor for QSO: Value between 0 to 1, by default 0.1
  -t PWEIGHT, --pweight PWEIGHT
                        Weighting factor for pseudo and amphiphlic pseudo amino acid composition: Value between 0 to 1, by default 0.05

#Parameters Description:

Input File: It allow users to provide input in two format; 
		i) FASTA format (standard) (e.g. protein.fa)
		ii) Simple Format, in this case, file should have sequences in a single line in single letter code (eg. protein.seq).

Output File: Program will save result in CSV format, in the provided filename.
		In case user do not provide output file name, it will be stored in pfeature_results.csv.
		In case user want to calculate all the features except AAI and ACR, the job name will be 'ALLCOMP'. Reason to leave AAI and ACR is, the feature calculation takes long time for longer sequences.

Job name: It allows users to choose the type of composition, the user want to calculate, such as AAC which stands for Amino Acid composition.
		In case user do not provide any job name, it will choose AAC by default.

N-terminal: It allows user to cut the specific number of residues from the N-terminal of the sequences.
		

C-terminal: It allows user to cut the specific number of residues from the C-terminal of the sequences.

NCT-terminal: It allows user to cut the specific number of residues from the N- and C-terminal of the sequences, and join them.

Rest_N : It allow users to drop the specific number of residues from N-terminal, and perform operations on the rest.

Rest_C : It allow users to drop the specific number of residues from C-terminal, and perform operations on the rest.

Split: It allow users to divided the sequence into number of sequences.

Lag : It defines the value for order of lag, lambda, gap or dipeptide, to calculate certain features.

Weight: It defines the weight factor to calculate the quasi-sequence order, by default it is set at 0.1.

Pweight: It defines the weight factor to calculate the pseudo and amphiphlic pseudo amino acid composition, by default it is set at 0.05.

#############################################################################################################################################################################################################################

2: Standalone for calculating binary profiles based features:

**Important: To run this script 'Data' folder should be in the same directory.**

Minimum USAGE: Minimum ussage is "pfeature_bin.py -i protein.fa" where protein.fa is a input fasta file. This will calculate the amino acid binary profile of the seqeunces provided in the fasta file. It will use other parameters by default. It will save output in "pfeature_result.csv" in CSV (comma seperated variables).

usage: pfeature_bin.py [-h] -i INPUT [-o OUTPUT]
                       [-j {AAB,DPB,ATB,BTB,PCB,AIB,ALLBIN}] [-n N_TERMINAL]
                       [-c C_TERMINAL] [-nct NC_TERMINAL] [-rn REST_N]
                       [-rc REST_C] [-s SPLIT] [-d LAG]

Please provide following arguments

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: protein or protein.sequence in FASTA format or single sequence per line in single letter code
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default pfeature_result.csv
  -j {AAB,DPB,ATB,BTB,PCB,AIB,ALLBIN}, --job {AAB,DPB,ATB,BTB,PCB,AIB,ALLBIN}
                        Job Type:
                        AAB: Amino acid based binary profile
                        DPB: Dipeptide based binary profile
                        ATB: Atom based binary profile
                        BTB: Bond based binary profile
                        PCB: Physico-chemical properties based binary profile
                        AIB: Amino-acid indices based binary profile
                        ALLBIN:All binary profiles together except ATB and BTB
                        by default AAB
  -n N_TERMINAL, --n_terminal N_TERMINAL
                        Window Length from N-terminal: by default 0
  -c C_TERMINAL, --c_terminal C_TERMINAL
                        Window Length from C-terminal: by default 0
  -nct NC_TERMINAL, --nc_terminal NC_TERMINAL
                        Residues from N- and C-terminal: by default 0
  -rn REST_N, --rest_n REST_N
                        Number of residues removed from N-terminal, by default 0
  -rc REST_C, --rest_c REST_C
                        Number of residues removed from C-terminal, by default 0
  -s SPLIT, --split SPLIT
                        Number of splits a sequence divided into, by default 0
  -d LAG, --lag LAG     This represents the order of gap, lag or dipeptide, by default 1

#Parameters Description:

Input File: It allow users to provide input in two format;
                i) FASTA format (standard) (e.g. protein.fa)
                ii) Simple Format, in this case, file should have sequences in a single line in single letter code (eg. protein.seq).

Output File: Program will save result in CSV format, in the provided filename.
                In case user do not provide output file name, it will be stored in pfeature_results.csv.
                In case user want to calculate all the features except ATB and BTB, the job name will be 'ALLBIN'. Reason to leave ATB and BTB is, the number of atoms and bonds are not equal in all amino acid residues.

Job name: It allows users to choose the type of composition, the user want to calculate, such as AAB which stands for Amino Acid based binary profile.
                In case user do not provide any job name, it will choose AAB by default.

N-terminal: It allows user to cut the specific number of residues from the N-terminal of the sequences.


C-terminal: It allows user to cut the specific number of residues from the C-terminal of the sequences.

NCT-terminal: It allows user to cut the specific number of residues from the N- and C-terminal of the sequences, and join them.

Rest_N : It allow users to drop the specific number of residues from N-terminal, and perform operations on the rest.

Rest_C : It allow users to drop the specific number of residues from C-terminal, and perform operations on the rest.

Split: It allow users to divided the sequence into number of sequences.

Lag : It defines the value for order of dipeptide, to calculate the dipeptide based binary profiles.

#############################################################################################################################################################################################################################

3: Standalone for calculating PSSM profile

**Important: To run this script a file with the file name 'envfile' is required.**

*This envfile contains paths for the following scripts/data:*
    i)   Path for blastpgp
    ii)  Path for blast database
    iii) Path for makemat


Minimum USAGE: Minimum ussage is "pfeature_pssm.py -i protein.fa" where protein.fa is a input fasta file. This will calculate the PSSM profile of the seqeunces provided in the fasta file. It will use other parameters by default. It will save output in "pssm_profile.csv" in CSV (comma seperated variables).

  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: protein or peptide sequence in FASTA format or single sequence per line in single letter code
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default pssm_profile.csv
  -n {N0,N1,N2,N3,N4}, --normalization_method {N0,N1,N2,N3,N4}
                        Normalization Method:
                        N0: It provides pssm profile without any normalization
                        N1: It normalizes pssm profile based on 1/(1+e^-x) formula
                        N2: It normalizes pssm profile based on (x-min)/(max-min) formula
                        N3: It normalizes pssm profile based on ((x-min)/(max-min))*100 formula
                        N4: It normalizes pssm profile based on 1/(1+e^-(x/100) formula
                        By default it is N0

#Parameters Description:

Input File: It allow users to provide input in two format;
                i) FASTA format (standard) (e.g. protein.fa)
                ii) Simple Format, in this case, file should have sequences in a single line in single letter code (eg. protein.seq).

Output File: Program will save result in CSV format, in the provided filename.
                In case user do not provide output file name, it will be stored in pssm_profile.csv

Normalization methods: It allows user to normalize the PSSM profiles using four different formula. The description is as follows:
                        N0: It provides pssm profile without any normalization
                        N1: It normalizes pssm profile based on 1/(1+e^-x) formula
                        N2: It normalizes pssm profile based on (x-min)/(max-min) formula
                        N3: It normalizes pssm profile based on ((x-min)/(max-min))*100 formula
                        N4: It normalizes pssm profile based on 1/(1+e^-(x/100) formula
                        By default it is N0 

Pfeature Packakage Files
=======================
It contantain following files, brief description of these files are given below:

LICENSE                  : License information

README.md                : This file provide information about this package

pfeature_comp.py         : Python program to calculate composition based features

pfeature_bin.py          : Python program to calculate binary profile based features

pfeature_pssm.py         : Python program to calculate pssm profile based features

protein.seq              : Example file contain protein sequences in simple format

protein.fa               : Example file contain protein sequences in FASTA format

Data                     : This folder contains the files required to calcuate the composition and binary profile based features.

envfile                  : This file contains the path information required to run the pfeature_pssm.py script.

Pfeature_Descriptors.pdf : This file comprises of description of the header of output files, which is generated using the aforementioned scripts.

Requirement.txt          : This file consists of commands and pre-requisite to run the aforementioned scripts.
