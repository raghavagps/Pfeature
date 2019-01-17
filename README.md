## Pfeature:Computation of features of peptides and protein
## Introduction
Pfeature is a comprehensive software developed for computing wide range of protein/discovered features that have been discovered over the past decades. It has following four major modules for computing protein features based on; i) Composition, ii) Binary profiles, iii) Evolutionary information and iv) Structure.  The composition based module allows user to compute; i) Simple compositions like amino acid, dipeptide, tripeptide; ii) Physicochemical properties based compositions; iii) Repeats and distribution of amino acids; iv) Shannon entropy to measure the low complexity regions; iv) Miscellaneous compositions like pseudo amino acid, autocorrelation, conjoint triad, quasi-sequence order. Binary profile of amino acid sequence provides complete information including order of residues or type of residues, which is not possible with composition based features. Thus, binary profile can be used to annotate protein at residue level. It is well established in literature that sequence profile based on evolutionary information provides more information then sequence itself.
We have developed number of isofoms of Pfeature that includes: i) A web server has been developed that uses Pfeature functions via web interface from https://webs.iiitd.edu.in/raghava/pfeature/ ; ii) Standalone version of Pfeature; iii) Library of python for Pfeature and iv) Python scripts for computing features. 
## Installation of Python Library
## Installation of Standalone (executables)
The prerequisite to run the standalone is pandas, numpy and python version above 3.6.<br/>

If pandas is not already installed, install it with command : pip3 install pandas<br/>

If numpy is not already installed, install it with command : pip3 install numpy<br/>

Install Pfeature<br/>

on Windows<br/>
&nbsp;&nbsp;(1) Download pfeature_win.zip<br/>
&nbsp;&nbsp;(2) Uncompress pfeature_win.zip file<br/>
&nbsp;&nbsp;(3) cd pfeature_win<br/>
&nbsp;&nbsp;(4) run pfeature.exe -i inputfile -o outputfile -m method<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;where the methods could be as follows:<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1 for Composition<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2 for Binary profiles<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3 for All<br/>
          
 On Mac<br/>
   (1) Download pfeature_mac.zip<br/>
   (2) Uncompress pfeature_win.zip file<br/>
   (3) cd pfeature_win<br/>
   (4) run ./pfeature -i inputfile -o outputfile -m method<br/>
       where the methods could be as follows:<br/>
          1 for Composition<br/>
          2 for Binary profiles<br/>
          3 for All<br/>

