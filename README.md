# Pfeature: Computation of features of peptides and proteins

## Introduction
Pfeature is a comprehensive software developed for computing wide range of protein/peptide features that have been discovered over the past decades. It has the following five major modules for computing protein features based on; i) Composition, ii) Binary profiles, iii) Evolutionary information iv) Structure and v) Pattern. The composition based module allows user to compute; i) Simple compositions like amino acid, dipeptide, tripeptide; ii) Physicochemical properties based compositions; iii) Repeats and distribution of amino acids; iv) Shannon entropy to measure the low complexity regions; iv) Miscellaneous compositions like pseudo amino acid, autocorrelation, conjoint triad, quasi-sequence order. Binary profile of amino acid sequence provides complete information including order of residues or type of residues, which is not possible with composition based features. Thus, binary profile can be used to annotate protein at residue level. It is well established in literature that sequence profile based on evolutionary information provides more information then sequence itself.

We have developed number of isoforms of Pfeature that include: i) A web server that uses Pfeature functions via web interface from https://webs.iiitd.edu.in/raghava/pfeature/ ; ii) Standalone version of Pfeature; iii) Library of python for Pfeature and iv) Python scripts for computing features.

## Documentation
One can read more about subroutines developed under Pfeature to compute wide range of proteins and peptide features from https://github.com/raghavagps/Pfeature/blob/master/Pfeature_Man.pdf . Further information is available from help page of web site https://webs.iiitd.edu.in/raghava/pfeature/help.php

## Web Service for Pfeature
A web server for computing wide range of protein and peptides features from their amino acid sequences. Following are main menus for computing features; i) Composition-based features, ii) Binary profile of sequences, iii) evolutionary information based features, iv) structural descriptors,and v) pattern based descriptors, for a group of protein/peptide sequences. Additionally, users will also be able to generate these features for sub-parts of protein/peptide sequences. Pfeature will be helpful to annotate structure, function and therapeutic properties of proteins/peptides.

**Available from URL: https://webs.iiitd.edu.in/raghava/pfeature/**

### Installation of Pfeature Library

### Prerequisite
The prerequisite to run the python library is pandas, numpy and python version above 3.6
pandas can be installed using following command: pip3 install pandas
numpy can be installed using following command: pip3 install numpy

### Steps for setting library
It has been tested on wide range of platforms that include Apple MAC, Windows and Linux (Ubuntu,Fedora). After installing pandas and numbpy user can install using following commands<br>

1.Download Pfeature from https://github.com/raghavagps/Pfeature/blob/master/PyLib/Pfeature.zip <br>
2. Extract or uncompress Pfeature.zip <br>
3. cd Pfeature <br>
4. python setup.py install <br>
<br>
## Installation of Pfeature Executables
In order to facilitate users, we created a single program of Pfeature which computes all possible descriptors for a protein/peptide sequence. We converted Python program in "C++" using Cython then "C++" codes are compiled to create executables. Pfeature executables have been created for number of platforms/operating systems like Apple MAC, Windows, Linux.
