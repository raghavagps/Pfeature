import io
import os
import shlex
import subprocess
import glob
import pandas as pd
import numpy as np
import sys
import uuid
import re
import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(description='Please provide following arguments',formatter_class=RawTextHelpFormatter)
parser.add_argument("-i","-I","--input", type=str, required=True, help="Input: protein or peptide sequence in FASTA format or single sequence per line in single letter code")
parser.add_argument("-o","-O","--output",type=str, help="Output: File for saving results by default pssm_profile.csv")
parser.add_argument("-n","-N","--normalization_method",type=str.upper, choices = ['N0','N1','N2','N3','N4'], help="Normalization Method:\nN0: It provides pssm profile without any normalization\nN1: It normalizes pssm profile based on 1/(1+e^-x) formula\nN2: It normalizes pssm profile based on (x-min)/(max-min) formula\nN3: It normalizes pssm profile based on ((x-min)/(max-min))*100 formula\nN4: It normalizes pssm profile based on 1/(1+e^-(x/100) formula\nBy default it is N0")
args = parser.parse_args()

# Input variable
sequencefile= args.input
# Output file
if args.output == None:
    result_filename= "pssm_profile.csv"
else:
    result_filename = args.output
if args.normalization_method == None:
    nm= "N0"
else:
    nm = args.normalization_method

def pssm_2(file):
    pd.options.display.float_format = '{:.2e}'.format
    df=file
    df1 = df
    df2 = df1.applymap(pssm_1)
    return df2
def pssm_1(x):
    if nm == 'N0':
        if type(x) is str:
            return x
        elif x:
            return x
        else:
            return
    if nm == 'N1':
        if type(x) is str:
            return x
        elif x < -700:
            return 0
        elif x:
            return (1/(1+(2.7182)**(-x)))
        else:
            return x
    if nm == 'N2':
        a = 1000
        b = -1000
        if type(x) is str:
            return x
        elif x:
            return (x-a)/(b - a)
        else:
            return
    if nm == 'N3':
        a = 1000
        b = -1000
        if type(x) is str:
            return x
        elif x:
            return ((x-a)*100)/(b - a)
        else:
            return
    if nm == 'N4':
        if type(x) is str:
            return x
        elif x:
            return (1/(1+(2.7182)**(-x/100)))
        else:
            return 
def readseq(file,out):
    with open(file) as f:
        records = f.read()
    records = records.split('>')[1:]
    seqid = []
    seq = []
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '', ''.join(array[1:]).upper())
        seqid.append(name)
        seq.append(sequence)
    if len(seqid) == 0:
        f=open(file,"r")
        data1 = f.readlines()
        for each in data1:
            seq.append(each.replace('\n',''))
        for i in range (1,len(seq)+1):
            seqid.append("Seq_"+str(i))
    for i in seq:
        if 'B' in i:
            print('\nError: The input sequences contain non-natural amino acids. Kindly check the sequence.\n')
            sys.exit()
        if 'J' in i:
            print('\nError: The input sequences contain non-natural amino acids. Kindly check the sequence.\n')
            sys.exit()
        if 'O' in i:
            print('\nError: The input sequences contain non-natural amino acids. Kindly check the sequence.\n')
            sys.exit()
        if 'U' in i:
            print('\nError: The input sequences contain non-natural amino acids. Kindly check the sequence.\n')
            sys.exit()
        if 'Z' in i:
            print('\nError: The input sequences contain non-natural amino acids. Kindly check the sequence.\n')
            sys.exit()
        if 'X' in i:
            print('\nError: The input sequences contain non-natural amino acids. Kindly check the sequence.\n')
            sys.exit()
    df4 = pd.DataFrame(seq)
    df4.to_csv(out,index=None,header=False)
def pssm(inputfile,outputfile):
    if os.path.exists('envfile'):
        with open('envfile', 'r') as file:
            data = file.readlines()
        output = []
        for line in data:
            if not "#" in line:
                output.append(line)
        if len(output)==3:
            paths = []
            for i in range (0,len(output)):
                paths.append(output[i].split(':')[1].replace('\n',''))
            blastpgp = paths[0]
            blastdb = paths[1]
            makemat = paths[2]
        if os.path.isfile(blastpgp) and os.access(blastpgp, os.R_OK):
            print('The provided directory for blastpgp is correct and readable.')
        else:
            print("########################################################################################################################")
            print("Error: Either 'blastbgp' file is missing from the provided directory in the 'envfile', or not readable. Kindly check.", file=sys.stderr)
            print("########################################################################################################################")
            sys.exit()
        if os.path.isfile(makemat) and os.access(makemat, os.R_OK):
            print('The provided directory for makemat is correct and readable.')
        else:
            print("########################################################################################################################")
            print("Error: Either 'makemat' file is missing from the provided directory in the 'envfile', or not readable. Kindly check.", file=sys.stderr)
            print("########################################################################################################################")
            sys.exit()
        if (glob.glob(blastdb+".pin")) and (glob.glob(blastdb+".psq")) and (glob.glob(blastdb+".phr")):
            print('The provided directory for blast database is correct and readable.')
        else:
            dbfiles = blastdb.split('/')[-1]
            print("##############################################################################################################################################################################################")
            print("Error: Either the files for BLAST database are missing from the provided directory in the 'envfile', or not readable. Please provide the files with extension of", dbfiles+".pin,", dbfiles+".phr,", dbfiles+".psq in the provided directory.", "Kindly check.", file=sys.stderr)
            print("##############################################################################################################################################################################################")
            sys.exit()

    else:
        print("####################################################################################")
        print("Error: Please provide the '{}', which comprises paths for BLASTPGP and MAKEMAT".format('envfile'), file=sys.stderr)
        print("####################################################################################")
        sys.exit()
    ss = []
    file_list = []
    readseq(inputfile,'temp.readseq')
    df1 = pd.read_csv('temp.readseq', header=None)
    aa = []
    for i in df1[0]:
        aa.append(len(i))
    for i in range(0,len(df1)):
        ss.append('>seq_'+str(i+1))
    df1['seq'] = ss
    df1 = df1[['seq',0]]
    for ii in range(0,len(df1)):
        name_file = df1['seq'][ii].replace('>','')+'.fasta'
        file_out = df1['seq'][ii].replace('>','')+'.pssmout'
        file_list.append(df1['seq'][ii].replace('>','')+'.mtx')
        df1.iloc[ii].to_csv(name_file, index=None, header=False, sep='\n')
        filename, file_extension = os.path.splitext(name_file)
        filename_o, file_extension_o = os.path.splitext(file_out)
        S =  blastpgp + ' -d ' +  blastdb + ' -i ' + name_file + ' -j 3 -C ' + file_out
        os.system(S)
        os.rename(file_out, filename_o + '.chk')
        outputfile1 = filename_o + ".chk"
        temp1 =filename_o + ".sn"
        C2 = 'echo {} > {}'.format(name_file,temp1)
        os.system(C2)
        temp2 = filename_o + ".pn"
        C1 = 'echo {} > {}'.format(outputfile1,temp2)
        os.system(C1)
        P = makemat + ' -P ' + filename_o
        os.system(P)
    dir = '.'
    tt = []
    for i in file_list:
        ss = []
        fp = open(i)
        all_line = fp.readlines()
        ss.append(all_line[14:])
        uu = []
        for j in all_line[14:]:
            uu.append(j.replace('\n','').replace('  ',','))
        np.savetxt(i+'_temp',uu, fmt="%s")
        df11 = pd.read_csv(i+'_temp', header=None, sep=",")
        col = [1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22]
        df13 = df11.iloc[:,col].reset_index(drop=True)
        df12 = pssm_2(df13)
        for i in range(0,len(df12)):
            tt.extend(df12.loc[i])
    bb = []
    cc = 0
    for i in range(0,len(aa)):
        bb.append(tt[cc:cc+20*(aa[i])])
        cc = cc+20*(aa[i])
    df123 = pd.DataFrame(bb)
    df456 = df123.fillna('NA')
    df456.to_csv(outputfile,index=None,header=False,float_format='%.2e')
    allfiles = os.listdir(dir)
    for item in allfiles :
        if item.endswith(".aux") or item.endswith(".sn") or item.endswith(".pn") or item.endswith(".mn") or item.endswith(".mtx_temp") or item.endswith(".mtx") or item.endswith(".chk") or item.endswith(".fasta") or item.endswith(".readseq"):
            os.remove(os.path.join(dir, item))

pssm(sequencefile,result_filename)
