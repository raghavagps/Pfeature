## Standalone Package of Pfeature
### Introduction:
Pfeature is a standalone software package for computing wide range of protein and peptides features from their amino acid
sequence. It has the following five major modules for computing protein features based on <br> i) Composition, ii) Binary profiles,
iii) Evolutionary information iv) Structure and v) Pattern.  We have developed number of forms of Pfeature that include: i) A web server
that uses Pfeature functions via web interface from https://webs.iiitd.edu.in/raghava/pfeature/ ; ii) Standalone version of Pfeature;
iii) Library of python for Pfeature and iv) Python scripts for computing features.

This is a standalone version of Pfeature (executables), we have generate executables of Pfeature that can be executed on different platform.

### Installation
Installation of Pfeature is simple as executables are available for  different operating systems (Windows, Ubuntu, Mac, Fedora, Centos). Following are main steps to install Pfeatures on different operating systems.
<br>
#### On Microsoft Windows:
1.      Download Pfeature_win.zip  from https://github.com/raghavagps/Pfeature/tree/master/exec/PFeature_win.zip
2.      unzip Pfeature_win.zip
3.      change directory to Pfeature_win
4.      Run the command: pfeature_win.exe -i <input_file> -o <output_file> -m <options>

#### On MacOs:<br>
1.      Download Pfeature_mac.zip  from https://github.com/raghavagps/Pfeature/tree/master/exec/PFeature_mac.zip <br>
2.      unzip Pfeature_mac.zip <br>
3.      change directory to Pfeature_mac <br>
4.      Run the command: ./pfeature_mac -i <input_file> -o <output_file> -m <options> <br>

#### On Ubuntu:<br>
1.      Download Pfeature_ubuntu.zip  from https://github.com/raghavagps/Pfeature/tree/master/exec/PFeature_ubuntu.zip <br>
2.      unzip Pfeature_ubuntu.zip <br>
3.      change directory to Pfeature_ubuntu <br>
4.      Run the command: ./pfeature_ubuntu -i <input_file> -o <output_file> -m <options> <br>

#### On Fedora:<br>
1.      Download Pfeature_fedora.zip  from https://github.com/raghavagps/Pfeature/tree/master/exec/PFeature_fedora.zip <br>
2.      unzip Pfeature_fedora.zip <br>
3.      change directory to Pfeature_fedora <br>
4.      Run the command: ./pfeature_fedora -i <input_file> -o <output_file> -m <options> <br>

#### On Centos: <br>
1.      Download Pfeature_cantos.zip  from https://github.com/raghavagps/Pfeature/tree/master/exec/PFeature_centos.zip <br>
2.      unzip Pfeature_centos.zip <br>
3.      change directory to Pfeature_centos <br>
4.      Run the command: ./pfeature_centos -i <input_file> -o <output_file> -m <options> <br>

### Folders & Files
Following is brief description of folders/files in folder Pfeature_OS.
Data: This folder contain csv files for different parameters required to run Pfeature <br>
example.seq : It is an example input file contain sequence of peptides (one sequence per line) <br>
example.out: An example output file corresponding to example.seq <br>
Features_Table.pdf : File provides description of features in different columns in output file <br>
lib: contain libraries required for package <br>
pfeature_os : Executable file, where os is operating system like win, mac, ubuntu. <br>
README : This file <br>
