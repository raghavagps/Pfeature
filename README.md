## Pfeature:Computation of features of peptides and protein
## Introduction
Pfeature is a comprehensive software developed for computing wide range of protein/discovered features that have been discovered over the past decades. It has following four major modules for computing protein features based on; i) Composition, ii) Binary profiles, iii) Evolutionary information and iv) Structure.  The composition based module allows user to compute; i) Simple compositions like amino acid, dipeptide, tripeptide; ii) Physicochemical properties based compositions; iii) Repeats and distribution of amino acids; iv) Shannon entropy to measure the low complexity regions; iv) Miscellaneous compositions like pseudo amino acid, autocorrelation, conjoint triad, quasi-sequence order. Binary profile of amino acid sequence provides complete information including order of residues or type of residues, which is not possible with composition based features. Thus, binary profile can be used to annotate protein at residue level. It is well established in literature that sequence profile based on evolutionary information provides more information then sequence itself.
We have developed number of isofoms of Pfeature that includes: i) A web server has been developed that uses Pfeature functions via web interface from https://webs.iiitd.edu.in/raghava/pfeature/ ; ii) Standalone version of Pfeature; iii) Library of python for Pfeature and iv) Python scripts for computing features. 
## Installation of Pfeature Library
### Prerequest
The prerequisite to run the python library is pandas, numpy and python version above 3.6<br>
pandas can be installed using following command: pip3 install pandas<br>
numpy can be installed using following command: pip3 install numpy<br>
### Steps for setting library
It has been tested on wide range of platforms that include Apple MAC, Windows and Linux (Ubuntu). After installing pandas and numbpy user can install using following commands<br>
1) Download Pfeature from https://github.com/raghavagps/Pfeature/blob/master/PyLib/Pfeature.zip 
2) Extract or uncompress Pfeature.zip
3) cd Prefature
4) python setup.py install

## Installation Executables of Pfeature
In order to facilitate users, we create single program of Pfeature which compute all possible descriptors for a protein/peptide sequence. We converted Python program in "C++" in using Cython then "C++" codes are compiled to create executables. Pfeature executables have been created for number of platforms/operating systems like Apple MAC, Windows, Linux. 
 
