#!/usr/bin/python
from __future__ import print_function
import getopt
import sys
import os
import numpy as np
import pandas as pd
import math
import itertools
from itertools import repeat
import argparse
import csv
from collections import Counter
import re
import glob
import time
from time import sleep
from tqdm import tqdm
from argparse import RawTextHelpFormatter
import uuid
import warnings
warnings.filterwarnings("ignore") 
std = list("ACDEFGHIKLMNPQRSTVWY")
PCP= pd.read_csv('Data/PhysicoChemical.csv', header=None)
AAindices = 'Data/aaind.txt'
AAIndex = pd.read_csv('Data/aaindex.csv',index_col='INDEX');
AAIndexNames = pd.read_csv('Data/AAIndexNames.csv',header=None);

parser = argparse.ArgumentParser(description='Please provide following arguments',formatter_class=RawTextHelpFormatter)
## Read Arguments from command
parser.add_argument("-i","-I","--input", type=str, required=True, help="Input: protein or peptide sequence in FASTA format or single sequence per line in single letter code")
parser.add_argument("-o","-O", "--output",type=str, help="Output: File for saving results by default pfeature_result.csv")
parser.add_argument("-j","-J", "--job",type=str.upper, choices = ['AAC','DPC','TPC','ATC','BTC','PCP','AAI','RRI','PRI','DDR','SEP','SER','SPC','ACR','CTC','CETD','PAAC','APAAC','QSO','SOC','ALLCOMP'], help="Job Type:\nAAC: Amino acid composition\nDPC: Dipeptide composition\nTPC: Tripeptide composition\nATC: Atomic composition\nBTC: Bond composition\nPCP: Physico-chemical properties composition\nAAI: Amino-acid indices composition\nRRI: Residue repeat information\nPRI: Physico-chemical properties repeat information\nDDR: Distance distribution of residues\nSEP: Shannon entropy of protein\nSER: Shannon entropy of residues\nSPC: Shannon entropy of physico-chemical properties\nACR: Autocorrelation descriptors\nCTC: Conjoint triad descriptors\nCeTD: Composition enhanced transition distribution\nPAAC: Pseudo amino acid composition\nAPAAC: Amphiphilic pseudo amino acid composition\nQSO: Quasi sequence order\nSOC: Sequence order coupling number\nALLCOMP:All composition features together except ACR and AAI\nby default AAC")
parser.add_argument("-n","-N","--n_terminal", type=int, help="Window Length from N-terminal: by default 0")
parser.add_argument("-c","-C","--c_terminal", type=int, help="Window Length from C-terminal: by default 0")
parser.add_argument("-nct","-Nct","-NCt","-NCT","-nCt","-ncT","-nCT","-NcT","--nc_terminal", type=int, help="Residues from N- and C-terminal: by default 0")
parser.add_argument("-rn","-Rn","-RN","-rN","--rest_n", type=int, help="Number of residues removed from N-terminal, by default 0")
parser.add_argument("-rc","-Rc","-RC","-rC","--rest_c", type=int, help="Number of residues removed from C-terminal, by default 0")
parser.add_argument("-rnc","-RNC","-Rnc","-rNc","-rnC","-RNc","-rNC","-RnC","--rest_nc", type=int, help="Number of residues removed from N- and C-terminal, by default 0")
parser.add_argument("-s","-S","--split", type=int, help="Number of splits a sequence divided into, by default 0")
parser.add_argument("-d","-D","--lag", type=int, help="This represents the order of gap, lag or dipeptide, by default 1")
parser.add_argument("-w","-W","--weight", type=float, help="Weighting Factor for QSO: Value between 0 to 1, by default 0.1")
parser.add_argument("-t","-T","--pweight", type=float, help="Weighting factor for pseudo and amphiphlic pseudo amino acid composition: Value between 0 to 1, by default 0.05")
args = parser.parse_args()

def aac_comp(file,out):
    filename, file_extension = os.path.splitext(file)
    f = open(out, 'w')
    sys.stdout = f
    df = pd.read_csv(file, header = None)
    zz = df.iloc[:,0]
    print("AAC_A,AAC_C,AAC_D,AAC_E,AAC_F,AAC_G,AAC_H,AAC_I,AAC_K,AAC_L,AAC_M,AAC_N,AAC_P,AAC_Q,AAC_R,AAC_S,AAC_T,AAC_V,AAC_W,AAC_Y,")
    for j in zz:
        for i in std:
            count = 0
            for k in j:
                temp1 = k
                if temp1 == i:
                    count += 1
                composition = (count/len(j))*100
            print("%.2f"%composition, end = ",")
        print("")
    f.truncate()
	
def aac_split(file,v,out):
    file1 = split(file,v)
    file2 = file1.iloc[:,0]
    file2.to_csv('sam_input.csv', index=None, header=None)
    aac_comp('sam_input.csv','tempfile_out')
    df4_1 = pd.read_csv("tempfile_out")
    df4 = df4_1.iloc[:,:-1]
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),(v*20)):
        pp.append(ss[i:i+(v*20)])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')
    os.remove('sam_input.csv')

	
def dpc_comp(file,q,out):
    filename, file_extension = os.path.splitext(file)
    f = open(out, 'w')
    sys.stdout = f
    df = pd.read_csv(file, header = None)
    df1 = pd.DataFrame(df[0].str.upper())
    zz = df1.iloc[:,0]
    for s in std:
        for u in std:
            print("DPC"+str(q)+"_"+s+u, end=',')
    print("")
    for i in range(0,len(zz)):
        for j in std:
            for k in std:
                count = 0
                temp = j+k
                for m3 in range(0,len(zz[i])-q):
                    b = zz[i][m3:m3+q+1:q]
                   # b.upper()

                    if b == temp:
                        count += 1
                    composition = (count/(len(zz[i])-(q)))*100
                print("%.2f" %composition, end = ',')
        print("")
    f.truncate()
def dpc_split(file,q,n,out):
    filename,file_ext = os.path.splitext(file)
    df1 = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df1[0].str.upper())
    k1 = []
    for w in range(0,len(df2)):
        s = 0
        k2 = []
        r = 0
        if len(df2[0][w])%n == 0:
            k2.extend(repeat(int(len(df2[0][w])/n),n))
        else:
            r = int(len(df2[0][w])%n)
            k2.extend(repeat(int(len(df2[0][w])/n),n-1))
            k2.append((int(len(df2[0][w])/n)+r))
        for j in k2:
            df3 = df2[0][w][s:j+s]
            k1.append(df3)
            s = j+s
    f = open(out, 'w')
    sys.stdout = f
    for h in range(1,n+1):
        for e in std:
            for r in std:
                print('DPC'+str(q)+'_'+e+r+'_s'+str(h), end=",")
    print("")
    for i in range(0,len(k1),n):
        k4 = k1[i:i+n]
        for j in k4:
            for k in std:
                for l in std:
                    count = 0
                    temp = k+l
                    for m3 in range(0,(len(j)-q)):
                        b = j[m3:m3+q+1:q]
                        if b == temp:
                            count += 1
                        composition = (count/(len(j)-q))*100
                    print("%.2f" %composition, end = ',')
        print("")
    f.truncate()

#############################TPC_COMP##############################
def tpc_comp(file,out):
    filename, file_extension = os.path.splitext(file)
    f = open(out, 'w')
    sys.stdout = f
    df = pd.read_csv(file, header = None)
    zz = df.iloc[:,0]
    for s in std:
       for u in std:
           for m in std:
               print('TPC_'+s+u+m, end=",")
    print("")
    for i in range(0,len(zz)):
        for j in std:
            for k in std:
                for m1 in std:
                    count = 0
                    temp = j+k+m1
                    for m3 in range(0,len(zz[i])):
                        b = zz[i][m3:m3+3]
                        if b == temp:
                            count += 1
                        composition = (count/(len(zz[i])-2))*100
                    print("%.2f" %composition, end = ',')
        print("")
    f.truncate()
def tpc_split(file,v,out):
    file1 = split(file,v)
    file2 = file1.iloc[:,0]
    file2.to_csv('sam_input.csv', index=None, header=None)
    tpc_comp('sam_input.csv','tempfile_out')
    df4_1 = pd.read_csv("tempfile_out")
    df4 = df4_1.iloc[:,:-1]
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),(v*8000)):
        pp.append(ss[i:i+(v*8000)])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')
    os.remove('sam_input.csv')
###########################atom###############
def atc(file,out):
    filename,file_ext = os.path.splitext(file)
    atom=pd.read_csv("Data/atom.csv",header=None)
    at=pd.DataFrame()
    i = 0
    C_atom = []
    H_atom = []
    N_atom = []
    O_atom = []
    S_atom = []

    while i < len(atom):
        C_atom.append(atom.iloc[i,1].count("C"))
        H_atom.append(atom.iloc[i,1].count("H"))
        N_atom.append(atom.iloc[i,1].count("N"))
        O_atom.append(atom.iloc[i,1].count("O"))
        S_atom.append(atom.iloc[i,1].count("S"))
        i += 1
    atom["C_atom"]=C_atom
    atom["O_atom"]=O_atom
    atom["H_atom"]=H_atom
    atom["N_atom"]=N_atom
    atom["S_atom"]=S_atom
##############read file ##########
    test1 = pd.read_csv(file,header=None)
    dd = []
    for i in range(0, len(test1)):
        dd.append(test1[0][i].upper())
    test = pd.DataFrame(dd)
    count_C = 0
    count_H = 0
    count_N = 0
    count_O = 0
    count_S = 0
    count = 0
    i1 = 0
    j = 0
    k = 0
    C_ct = []
    H_ct = []
    N_ct = []
    O_ct = []
    S_ct = []
    while i1 < len(test) :
        while j < len(test[0][i1]) :
            while k < len(atom) :
                if test.iloc[i1,0][j]==atom.iloc[k,0].replace(" ","") :
                    count_C = count_C + atom.iloc[k,2]
                    count_H = count_H + atom.iloc[k,3]
                    count_N = count_N + atom.iloc[k,4]
                    count_O = count_O + atom.iloc[k,5]
                    count_S = count_S + atom.iloc[k,6]
                #count = count_C + count_H + count_S + count_N + count_O
                k += 1
            k = 0
            j += 1
        C_ct.append(count_C)
        H_ct.append(count_H)
        N_ct.append(count_N)
        O_ct.append(count_O)
        S_ct.append(count_S)
        count_C = 0
        count_H = 0
        count_N = 0
        count_O = 0
        count_S = 0
        j = 0
        i1 += 1
    test["C_count"]=C_ct
    test["H_count"]=H_ct
    test["N_count"]=N_ct
    test["O_count"]=O_ct
    test["S_count"]=S_ct

    ct_total = []
    m = 0
    while m < len(test) :
        ct_total.append(test.iloc[m,1] + test.iloc[m,2] + test.iloc[m,3] + test.iloc[m,4] + test.iloc[m,5])
        m += 1
    test["count"]=ct_total
##########final output#####
    final = pd.DataFrame()
    n = 0
    p = 0
    C_p = []
    H_p = []
    N_p = []
    O_p = []
    S_p = []
    while n < len(test):
        C_p.append((test.iloc[n,1]/test.iloc[n,6])*100)
        H_p.append((test.iloc[n,2]/test.iloc[n,6])*100)
        N_p.append((test.iloc[n,3]/test.iloc[n,6])*100)
        O_p.append((test.iloc[n,4]/test.iloc[n,6])*100)
        S_p.append((test.iloc[n,5]/test.iloc[n,6])*100)
        n += 1
    final["ATC_C"] = C_p
    final["ATC_H"] = H_p
    final["ATC_N"] = N_p
    final["ATC_O"] = O_p
    final["ATC_S"] = S_p

    (final.round(2)).to_csv(out, index = None, encoding = 'utf-8')
########################################atom#################
def atc_split(file,N,out):
    import pandas as pd
    filename,file_ext = os.path.splitext(file)
    atom=pd.read_csv("Data/atom.csv",header=None)
    #atom=pd.read_csv("atom.csv",header=None)
    at=pd.DataFrame()
    i = 0
    C_atom = []
    H_atom = []
    N_atom = []
    O_atom = []
    S_atom = []

    while i < len(atom):
        C_atom.append(atom.iloc[i,1].count("C"))
        H_atom.append(atom.iloc[i,1].count("H"))
        N_atom.append(atom.iloc[i,1].count("N"))
        O_atom.append(atom.iloc[i,1].count("O"))
        S_atom.append(atom.iloc[i,1].count("S"))
        i += 1
    atom["C_atom"]=C_atom
    atom["H_atom"]=H_atom
    atom["N_atom"]=N_atom
    atom["O_atom"]=O_atom
    atom["S_atom"]=S_atom

##############read file ##########
    data1 = pd.read_csv(file,header=None)
    data = pd.DataFrame(data1[0].str.upper())
    k1 = []
    for e in range(0,len(data)):
        s = 0
        k2 = []
        r = 0
        if len(data[0][e])%N == 0:
            k2.extend(repeat(int(len(data[0][e])/N),N))
        else:
            r = int(len(data[0][e])%N)
            k2.extend(repeat(int(len(data[0][e])/N),N-1))
            k2.append((int(len(data[0][e])/N))+r)
        for j in k2:
            df3 = data[0][e][s:j+s]
            k1.append(df3)
            s = j+s
    test = pd.DataFrame(k1)
    count_C = 0
    count_H = 0
    count_N = 0
    count_O = 0
    count_S = 0
    count = 0
    i1 = 0
    j = 0
    k = 0
    C_ct = []
    H_ct = []
    N_ct = []
    O_ct = []
    S_ct = []
    while i1 < len(test) :
        while j < len(test[0][i1]) :
            while k < len(atom) :
                if test.iloc[i1,0][j]==atom.iloc[k,0].replace(" ","") :
                    count_C = count_C + atom.iloc[k,2]
                    count_H = count_H + atom.iloc[k,3]
                    count_N = count_N + atom.iloc[k,4]
                    count_O = count_O + atom.iloc[k,5]
                    count_S = count_S + atom.iloc[k,6]
                k += 1
            k = 0
            j += 1
        C_ct.append(count_C)
        H_ct.append(count_H)
        N_ct.append(count_N)
        O_ct.append(count_O)
        S_ct.append(count_S)
        count_C = 0
        count_H = 0
        count_N = 0
        count_O = 0
        count_S = 0
        j = 0
        i1 += 1
    test["C_count"]=C_ct
    test["H_count"]=H_ct
    test["N_count"]=N_ct
    test["O_count"]=O_ct
    test["S_count"]=S_ct

    ct_total = []
    m = 0
    while m < len(test) :
        ct_total.append(test.iloc[m,1] + test.iloc[m,2] + test.iloc[m,3] + test.iloc[m,4] + test.iloc[m,5])
        m += 1
    test["count"]=ct_total
##########final output#####
    final = pd.DataFrame()
    n = 0
    p = 0
    C_p = []
    H_p = []
    N_p = []
    O_p = []
    S_p = []
    while n < len(test):
        C_p.append((test.iloc[n,1]/test.iloc[n,6])*100)
        H_p.append((test.iloc[n,2]/test.iloc[n,6])*100)
        N_p.append((test.iloc[n,3]/test.iloc[n,6])*100)
        O_p.append((test.iloc[n,4]/test.iloc[n,6])*100)
        S_p.append((test.iloc[n,5]/test.iloc[n,6])*100)
        n += 1
    final["C"] = C_p
    final["H"] = H_p
    final["N"] = N_p
    final["O"] = O_p
    final["S"] = S_p

    df3 = final
    bb = []
    for i in range(0,len(df3),N):
        aa = []
        for j in range(N):
            aa.extend(df3.loc[i+j])
        bb.append(aa)
    zz = pd.DataFrame(bb)
    header = []
    head = ['ATC_C_s','ATC_H_s','ATC_N_s','ATC_O_s','ATC_S_s']
    for e in range(1,N+1):
        for t in head:
            header.append(t+str(e))
    zz.columns=header
    (zz.round(2)).to_csv(out, index=None)

#########################################bond#######################################################
def bond(file,out) :
    tota = []
    hy = []
    Si = []
    Du = []
    b1 = []
    b2 = []
    b3 = []
    b4 = []
    bb = pd.DataFrame()
    filename, file_extension = os.path.splitext(file)
    df = pd.read_csv(file, header = None)
    bonds=pd.read_csv("Data/bonds.csv", sep = ",")
    for i in range(0,len(df)) :
        tot = 0
        h = 0
        S = 0
        D = 0
        tota.append([i])
        hy.append([i])
        Si.append([i])
        Du.append([i])
        for j in range(0,len(df[0][i])) :
            temp = df[0][i][j]
            for k in range(0,len(bonds)) :
                if bonds.iloc[:,0][k] == temp :
                    tot = tot + bonds.iloc[:,1][k]
                    h = h + bonds.iloc[:,2][k]
                    S = S + bonds.iloc[:,3][k]
                    D = D + bonds.iloc[:,4][k]
        tota[i].append(tot)
        hy[i].append(h)
        Si[i].append(S)
        Du[i].append(D)
    for m in range(0,len(df)) :
        b1.append(tota[m][1])
        b2.append(hy[m][1])
        b3.append(Si[m][1])
        b4.append(Du[m][1])

    bb["BTC_T"] = b1
    bb["BTC_H"] = b2
    bb["BTC_S"] = b3
    bb["BTC_D"] = b4

    bb.to_csv(out, index=None, encoding="utf-8")
#########################################bond_split#########################
def bond_split(file,n,out) :
    tota = []
    hy = []
    Si = []
    Du = []
    b1 = []
    b2 = []
    b3 = []
    b4 = []
    bb = pd.DataFrame()
    filename, file_extension = os.path.splitext(file)
    #df = pd.read_csv(file, header = None)
    df = split(file,n)
    bonds=pd.read_csv("Data/bonds.csv")
    for i in range(0,len(df)) :
        tot = 0
        h = 0
        S = 0
        D = 0
        tota.append([i])
        hy.append([i])
        Si.append([i])
        Du.append([i])
        for j in range(0,len(df[0][i])) :
            temp = df[0][i][j]
            for k in range(0,len(bonds)) :
                if bonds.iloc[:,0][k] == temp :
                    tot = tot + bonds.iloc[:,1][k]
                    h = h + bonds.iloc[:,2][k]
                    S = S + bonds.iloc[:,3][k]
                    D = D + bonds.iloc[:,4][k]
        tota[i].append(tot)
        hy[i].append(h)
        Si[i].append(S)
        Du[i].append(D)
    for m in range(0,len(df)) :
        b1.append(tota[m][1])
        b2.append(hy[m][1])
        b3.append(Si[m][1])
        b4.append(Du[m][1])

    bb["total number of bonds"] = b1
    bb["hydrogen bonds"] = b2
    bb["single bonds"] = b3
    bb["double bonds"] = b4
    header = []
    header1 = ('BTC_T','BTC_H','BTC_S','BTC_D')
    for i in range(1,n+1):
        for j in header1:
            header.append(j+"_s"+str(i))
    qq = []
    for i in range(0,len(bb),n):
        aa = []
        for j in range(n):
            aa.extend(bb.loc[i+j])
        qq.append(aa)
    zz = pd.DataFrame(qq)
    zz.columns = header
    zz.to_csv(out, index=None)

############################PhysicoChemical Properties###################################
PCP= pd.read_csv('Data/PhysicoChemical.csv', header=None)

headers = ['PCP_PC','PCP_NC','PCP_NE','PCP_PO','PCP_NP','PCP_AL','PCP_CY','PCP_AR','PCP_AC','PCP_BS','PCP_NE_pH','PCP_HB','PCP_HL','PCP_NT','PCP_HX','PCP_SC','PCP_SS_HE','PCP_SS_ST','PCP_SS_CO','PCP_SA_BU','PCP_SA_EX','PCP_SA_IN','PCP_TN','PCP_SM','PCP_LR','PCP_Z1','PCP_Z2','PCP_Z3','PCP_Z4','PCP_Z5'];

def encode(peptide):
    l=len(peptide);
    encoded=np.zeros(l);
    for i in range(l):
        if(peptide[i]=='A'):
            encoded[i] = 0;
        elif(peptide[i]=='C'):
            encoded[i] = 1;
        elif(peptide[i]=='D'):
            encoded[i] = 2;
        elif(peptide[i]=='E'):
            encoded[i] = 3;
        elif(peptide[i]=='F'):
            encoded[i] = 4;
        elif(peptide[i]=='G'):
            encoded[i] = 5;
        elif(peptide[i]=='H'):
            encoded[i] = 6;
        elif(peptide[i]=='I'):
            encoded[i] = 7;
        elif(peptide[i]=='K'):
            encoded[i] = 8;
        elif(peptide[i]=='L'):
            encoded[i] = 9;
        elif(peptide[i]=='M'):
            encoded[i] = 10;
        elif(peptide[i]=='N'):
            encoded[i] = 11;
        elif(peptide[i]=='P'):
            encoded[i] = 12;
        elif(peptide[i]=='Q'):
            encoded[i] = 13;
        elif(peptide[i]=='R'):
            encoded[i] = 14;
        elif(peptide[i]=='S'):
            encoded[i] = 15;
        elif(peptide[i]=='T'):
            encoded[i] = 16;
        elif(peptide[i]=='V'):
            encoded[i] = 17;
        elif(peptide[i]=='W'):
            encoded[i] = 18;
        elif(peptide[i]=='Y'):
            encoded[i] = 19;
        else:
            print('Wrong residue!');
    return encoded;
def lookup(peptide,featureNum):
    l=len(peptide);
    peptide = list(peptide);
    out=np.zeros(l);
    peptide_num = encode(peptide);

    for i in range(l):
        out[i] = PCP[peptide_num[i]][featureNum];
    return sum(out);
def pcp_1(file,out123):

    if(type(file) == str):
        seq = pd.read_csv(file,header=None);
        #seq=seq.T
        seq[0].values.tolist()
        seq=seq[0];
    else:
        seq  = file;

    l = len(seq);

    rows = PCP.shape[0]; # Number of features in our reference table
    col = 20 ; # Denotes the 20 amino acids

    seq=[seq[i].upper() for i in range(l)]
    sequenceFeature = [];
    sequenceFeature.append(headers); #To put property name in output csv

    for i in range(l): # Loop to iterate over each sequence
        nfeatures = rows;
        sequenceFeatureTemp = [];
        for j in range(nfeatures): #Loop to iterate over each feature
            featureVal = lookup(seq[i],j)
            if(len(seq[i])!=0):
                sequenceFeatureTemp.append(round(featureVal/len(seq[i]),3));
            else:
                sequenceFeatureTemp.append('NaN')

        sequenceFeature.append(sequenceFeatureTemp);

    out = pd.DataFrame(sequenceFeature);
    file = open(out123,'w')
    with file:
        writer = csv.writer(file);
        writer.writerows(sequenceFeature);
    return sequenceFeature;

###############################################################################################################################################
def pcp_split(file,v,out):
    file1 = split(file,v)
    file2 = file1.iloc[:,0]
    pcp_1(file2,'tempfile_out')
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),(v*30)):
        pp.append(ss[i:i+(v*30)])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')
###############################RRI#################################
def RAAC(file,out):
    filename, file_extension = os.path.splitext(file)
    df = pd.read_csv(file, header = None)
    df1 = pd.DataFrame(df[0].str.upper())
    count = 0
    cc = []
    i = 0
    x = 0
    temp = pd.DataFrame()
    f = open(out,'w')
    sys.stdout = f
    print("RRI_A,RRI_C,RRI_D,RRI_E,RRI_F,RRI_G,RRI_H,RRI_I,RRI_K,RRI_L,RRI_M,RRI_N,RRI_P,RRI_Q,RRI_R,RRI_S,RRI_T,RRI_V,RRI_W,RRI_Y,")
    for q in range(0,len(df1)):
        while i < len(std):
            cc = []
            for j in df1[0][q]:
                if j == std[i]:
                    count += 1
                    cc.append(count)
                else:
                    count = 0
            while x < len(cc) :
                if x+1 < len(cc) :
                    if cc[x]!=cc[x+1] :
                        if cc[x] < cc[x+1] :
                            cc[x]=0
                x += 1
            cc1 = [e for e in cc if e!= 0]
            cc = [e*e for e in cc if e != 0]
            zz= sum(cc)
            zz1 = sum(cc1)
            if zz1 != 0:
                zz2 = zz/zz1
            else:
                zz2 = 0
            print("%.2f"%zz2,end=',')
            i += 1
        i = 0
        print(" ")
    f.truncate()		

def RAAC_split(file,n,out):
    filename, file_extension = os.path.splitext(file)
    file1 = split(file,n)
    data = file1
    count = 0
    cc = []
    i = 0
    x = 0
    temp = pd.DataFrame()
    f = open("sam.raac_split",'w')
    sys.stdout = f
    print("RRI_A,RRI_C,RRI_D,RRI_E,RRI_F,RRI_G,RRI_H,RRI_I,RRI_K,RRI_L,RRI_M,RRI_N,RRI_P,RRI_Q,RRI_R,RRI_S,RRI_T,RRI_V,RRI_W,RRI_Y,")
    for q in range(0,len(data)):
        mm = data[0][q]
        while i < len(std):
            cc = []
            for j in mm:
                if j == std[i]:
                    count += 1
                    cc.append(count)
                else:
                    count = 0
            while x < len(cc) :
                if x+1 < len(cc) :
                    if cc[x]!=cc[x+1] :
                        if cc[x] < cc[x+1] :
                            cc[x]=0
                x += 1
            cc1 = [e for e in cc if e!= 0]
            cc = [e*e for e in cc if e != 0]
            zz= sum(cc)
            zz1 = sum(cc1)
            if zz1 != 0:
                zz2 = zz/zz1
            else:
                zz2 = 0
            print("%.2f"%zz2,end=',')
            i += 1
        i = 0
        print(" ")
    f.truncate()
    df = pd.read_csv("sam.raac_split")
    df1 = df.iloc[:,:-1]
    header= []
    for i in range(1,n+1):
        for j in std:
            header.append('RRI_'+j+'_s'+str(i))
    ss = []
    for i in range(0,len(df1)):
        ss.extend(df1.loc[i])
    pp = []
    for i in range(0,len(ss),(n*20)):
        pp.append(ss[i:i+(n*20)])
    df5= pd.DataFrame(pp)
    df5.columns = header
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('sam.raac_split')
##########################################PRI###########################
PCP= pd.read_csv('Data/PhysicoChemical.csv', header=None) #Our reference table for properties
headers_1 = ['PRI_PC','PRI_NC','PRI_NE','PRI_PO','PRI_NP','PRI_AL','PRI_CY','PRI_AR','PRI_AC','PRI_BS','PRI_NE_pH','PRI_HB','PRI_HL','PRI_NT','PRI_HX','PRI_SC','PRI_SS_HE','PRI_SS_ST','PRI_SS_CO','PRI_SA_BU','PRI_SA_EX','PRI_SA_IN','PRI_TN','PRI_SM','PRI_LR'];
def encode(peptide):
    l=len(peptide);
    encoded=np.zeros(l);
    for i in range(l):
        if(peptide[i]=='A'):
            encoded[i] = 0;
        elif(peptide[i]=='C'):
            encoded[i] = 1;
        elif(peptide[i]=='D'):
            encoded[i] = 2;
        elif(peptide[i]=='E'):
            encoded[i] = 3;
        elif(peptide[i]=='F'):
            encoded[i] = 4;
        elif(peptide[i]=='G'):
            encoded[i] = 5;
        elif(peptide[i]=='H'):
            encoded[i] = 6;
        elif(peptide[i]=='I'):
            encoded[i] = 7;
        elif(peptide[i]=='K'):
            encoded[i] = 8;
        elif(peptide[i]=='L'):
            encoded[i] = 9;
        elif(peptide[i]=='M'):
            encoded[i] = 10;
        elif(peptide[i]=='N'):
            encoded[i] = 11;
        elif(peptide[i]=='P'):
            encoded[i] = 12;
        elif(peptide[i]=='Q'):
            encoded[i] = 13;
        elif(peptide[i]=='R'):
            encoded[i] = 14;
        elif(peptide[i]=='S'):
            encoded[i] = 15;
        elif(peptide[i]=='T'):
            encoded[i] = 16;
        elif(peptide[i]=='V'):
            encoded[i] = 17;
        elif(peptide[i]=='W'):
            encoded[i] = 18;
        elif(peptide[i]=='Y'):
            encoded[i] = 19;
        else:
            print(peptide[i], ' is a wrong residue!');
    return encoded;
def lookup_1(peptide,featureNum):
    l=len(peptide);
    peptide = list(peptide);
    out=[];
    peptide_num = encode(peptide);

    for i in range(l):
        out.append(PCP[peptide_num[i]][featureNum]);
    return out;
def binary_profile_1(file,featureNumb):
    if(type(file) == str):
        seq = pd.read_csv(file,header=None, sep=',');
        seq=seq.T
        seq[0].values.tolist()
        seq=seq[0];
    else:
        seq  = file;
    l = len(seq);
    bin_prof = [];
    for i in range(0,l):
        temp = lookup_1(seq[i],featureNumb);
        bin_prof.append(temp);
    return bin_prof;
def repeats(file,out123):
    if(type(file) == str):
        seq = pd.read_csv(file,header=None, sep=',');
        #seq=seq.T
        seq[0].values.tolist()
        seq=seq[0];
    else:
        seq  = file;
    seq=[seq[i].upper() for i in range(len(seq))]
    dist =[];
    dist.append(headers_1);
    l = len(seq);
    for i in range(l):
        temp=[];
        for j in range(25):
            bin_prof = binary_profile_1(seq, j);
            if(j>=25):
                print('Error! Feature Number must be between 0-24');
                break;
            k=0;
            num=0;
            denom=0;
            ones=0;
            zeros=0;
            for j in range(len(bin_prof[i])):
                if(bin_prof[i][j]==0):
                    num+=k*k;
                    denom+=k;
                    k=0;
                    zeros+=1;
                elif(j==len(bin_prof[i])-1):
                    k+=1;
                    num+=k*k;
                    denom+=k;
                else:
                    k+=1;
                    ones+=1;
            if(ones!=0):
                answer = num/(ones*ones)
                temp.append(round(num/(ones*ones),2));
            elif(ones==0):
                temp.append(0);
        dist.append(temp)
    out = pd.DataFrame(dist)
    file1 = open(out123,'w')
    with file1:
        writer = csv.writer(file1);
        writer.writerows(dist);
    return out

def repeats_split(file,v,out):
    file1 = split(file,v)
    file2 = file1.iloc[:,0]
    repeats(file2,'tempfile_out')
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),(v*25)):
        pp.append(ss[i:i+(v*25)])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')
	
#################################DDOR####################################
def DDOR(file,out) :
    df = pd.read_csv(file, header = None)
    df1 = pd.DataFrame(df[0].str.upper())
    f = open(out,'w')
    sys.stdout = f
    for i in std:
        print('DDR_'+i, end=",")
    print("")
    for i in range(0,len(df1)):
        s = df1[0][i]
        p = s[::-1]
        for j in std:
            zz = ([pos for pos, char in enumerate(s) if char == j])
            pp = ([pos for pos, char in enumerate(p) if char == j])
            ss = []
            for i in range(0,(len(zz)-1)):
                ss.append(zz[i+1] - zz[i]-1)
            if zz == []:
                ss = []
            else:
                ss.insert(0,zz[0])
                ss.insert(len(ss),pp[0])
            cc1=  (sum([e for e in ss])+1)
            cc = sum([e*e for e in ss])
            zz2 = cc/cc1
            print("%.2f"%zz2,end=",")
        print("")
    f.truncate()

def DDOR_split(file,v,out):
    file1 = split(file,v)
    df = file1
    df1 = pd.DataFrame(df[0].str.upper())
    f = open('tempfile_out','w')
    sys.stdout = f
    for i in std:
        print('DDR_'+i, end=",")
    print("")
    for i in range(0,len(df1)):
        s = df1[0][i]
        p = s[::-1]
        for j in std:
            zz = ([pos for pos, char in enumerate(s) if char == j])
            pp = ([pos for pos, char in enumerate(p) if char == j])
            ss = []
            for i in range(0,(len(zz)-1)):
                ss.append(zz[i+1] - zz[i]-1)
            if zz == []:
                ss = []
            else:
                ss.insert(0,zz[0])
                ss.insert(len(ss),pp[0])
            cc1=  (sum([e for e in ss])+1)
            cc = sum([e*e for e in ss])
            zz2 = cc/cc1
            print("%.2f"%zz2,end=",")
        print("")
    f.truncate()
    df3 = pd.read_csv("tempfile_out")
    df4 = df3.iloc[:,:-1]
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),(v*20)):
        pp.append(ss[i:i+(v*20)])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')

########################Shannon_Entropy Whole protein######################################
def entropy_single(seq):
    seq=seq.upper()
    num, length = Counter(seq), len(seq)
    return -sum( freq/length * math.log(freq/length, 2) for freq in num.values())

def SE(filename,out):
    data=list((pd.read_csv(filename,sep=',',header=None)).iloc[:,0])
#     print(data)
    Val=[]
    header=["SEP"]
    for i in range(len(data)):
        data1=''
        data1=str(data[i])
        data1=data1.upper()
        allowed = set(('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'))
        is_data_invalid = set(data1).issubset(allowed)
        if is_data_invalid==False:
            print("Error: Please check for invalid inputs in the sequence.","\nError in: ","Sequence number=",i+1,",","Sequence = ",data[i],",","\nNOTE: Spaces, Special characters('[@_!#$%^&*()<>?/\|}{~:]') and Extra characters(BJOUXZ) should not be there.")
            return
        Val.append(round((entropy_single(str(data[i]))),3))
        #print(Val[i])
        file= open(out,'w', newline='\n')#output file
        with file:
            writer=csv.writer(file,delimiter='\n');
            writer.writerow(header)
            writer.writerow(Val);
    return Val

def SE_split(file,v,out):
    SE(file,'tempfile_out')
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),v):
        pp.append(ss[i:i+v])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')

################################Shannon_Entropy residue#############################################
def SE_residue_level(filename,out):
    data=list((pd.read_csv(filename,sep=',',header=None)).iloc[:,0])
    data2=list((pd.read_csv(filename,sep=',',header=None)).iloc[:,0])
    Val=np.zeros(len(data))
    GH=[]
    for i in range(len(data)):
        my_list={'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'W':0,'Y':0}
        data1=''
        data1=str(data[i])
        data1=data1.upper()
        allowed = set(('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'))
        is_data_invalid = set(data1).issubset(allowed)
        if is_data_invalid==False:
            print("Error: Please check for invalid inputs in the sequence.","\nError in: ","Sequence number=",i+1,",","Sequence = ",data[i],",","\nNOTE: Spaces, Special characters('[@_!#$%^&*()<>?/\|}{~:]') and Extra characters(BJOUXZ) should not be there.")
            return
        seq=data[i]
        seq=seq.upper()
        num, length = Counter(seq), len(seq)
        num=dict(sorted(num.items()))
        C=list(num.keys())
        F=list(num.values())
        for key, value in my_list.items():
             for j in range(len(C)):
                if key == C[j]:
                    my_list[key] = round(((F[j]/length)* math.log(F[j]/length, 2)),3)
        GH.append(list(my_list.values()))
    file= open(out,'w', newline='')#output file
    with file:
        writer=csv.writer(file);
        writer.writerow(('SER_A','SER_C','SER_D','SER_E','SER_F','SER_G','SER_H','SER_I','SER_K','SER_L','SER_M','SER_N','SER_P','SER_Q','SER_R','SER_S','SER_T','SER_V','SER_W','SER_Y'));
        writer.writerows(GH);
    return(GH)

def SE_residue_level_split(file,v,out):
    SE_residue_level(file,'tempfile_out')
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),(v*20)):
        pp.append(ss[i:i+(v*20)])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')

#########################Shanon entropy for PCP#################################
def lookup(peptide,featureNum):
    l=len(peptide);
    peptide = list(peptide);
    out=np.zeros(l);
    peptide_num = encode(peptide);
    for i in range(l):
        out[i] = PCP[peptide_num[i]][featureNum];
    return sum(out);
def pcp(file):
    SEP_headers = ['SEP_PC','SEP_NC','SEP_NE','SEP_PO','SEP_NP','SEP_AL','SEP_CY','SEP_AR','SEP_AC','SEP_BS','SEP_NE_pH','SEP_HB','SEP_HL','SEP_NT','SEP_HX','SEP_SC','SEP_SS_HE','SEP_SS_ST','SEP_SS_CO','SEP_SA_BU','SEP_SA_EX','SEP_SA_IN','SEP_TN','SEP_SM','SEP_LR']
    if(type(file) == str):
        seq = pd.read_csv(file,header=None, sep=',');
        seq=seq.T
        seq[0].values.tolist()
        seq=seq[0];
    else:
        seq  = file;
    l = len(seq);
    rows = PCP.shape[0]; # Number of features in our reference table
    col = 20 ; # Denotes the 20 amino acids
    seq=[seq[i].upper() for i in range(l)]
    sequenceFeature = [];
    sequenceFeature.append(SEP_headers); #To put property name in output csv
    
    for i in range(l): # Loop to iterate over each sequence
        nfeatures = rows;
        sequenceFeatureTemp = [];
        for j in range(nfeatures): #Loop to iterate over each feature
            featureVal = lookup(seq[i],j)   
            if(len(seq[i])!=0):
                sequenceFeatureTemp.append(featureVal/len(seq[i]))
            else:
                sequenceFeatureTemp.append('NaN')
        sequenceFeature.append(sequenceFeatureTemp);
    out = pd.DataFrame(sequenceFeature);
    return sequenceFeature;
def phyChem(file,mode='all',m=0,n=0):
    if(type(file) == str):
        seq1 = pd.read_csv(file,header=None, sep=',');
        seq1 = pd.DataFrame(seq1[0].str.upper())
        seq=[]
        [seq.append(seq1.iloc[i][0]) for i in range(len(seq1))]
    else:
        seq  = file;
    l = len(seq);
    newseq = [""]*l; # To store the n-terminal sequence
    for i in range(0,l):
        l = len(seq[i]);
        if(mode=='NT'):
            n=m;
            if(n!=0):
                newseq[i] = seq[i][0:n];
            elif(n>l):
                print('Warning! Sequence',i,"'s size is less than n. The output table would have NaN for this sequence");
            else:
                print('Value of n is mandatory, it cannot be 0')
                break;
        elif(mode=='CT'):
            n=m;
            if(n!=0):
                newseq[i] = seq[i][(len(seq[i])-n):]
            elif(n>l):
                print('WARNING: Sequence',i+1,"'s size is less than the value of n given. The output table would have NaN for this sequence");
            else:
                print('Value of n is mandatory, it cannot be 0')
                break;
        elif(mode=='all'):
            newseq = seq;
        elif(mode=='rest'):
            if(m==0):
                print('Kindly provide start index for rest, it cannot be 0');
                break;
            else:
                if(n<=len(seq[i])):
                    newseq[i] = seq[i][m-1:n+1]
                elif(n>len(seq[i])):
                    newseq[i] = seq[i][m-1:len(seq[i])]
                    print('WARNING: Since input value of n for sequence',i+1,'is greater than length of the protein, entire sequence starting from m has been considered')
        else:
            print("Wrong Mode. Enter 'NT', 'CT','all' or 'rest'");        
    output = pcp(newseq);
    return output
def shannons(filename,out123):
    SEP_headers = ['SEP_PC','SEP_NC','SEP_NE','SEP_PO','SEP_NP','SEP_AL','SEP_CY','SEP_AR','SEP_AC','SEP_BS','SEP_NE_pH','SEP_HB','SEP_HL','SEP_NT','SEP_HX','SEP_SC','SEP_SS_HE','SEP_SS_ST','SEP_SS_CO','SEP_SA_BU','SEP_SA_EX','SEP_SA_IN','SEP_TN','SEP_SM','SEP_LR']
    if(type(filename) == str):
        seq1 = pd.read_csv(filename,header=None, sep=',');
        seq1 = pd.DataFrame(seq1[0].str.upper())
    else:
        seq1  = filename;
    seq=[]
    [seq.append(seq1.iloc[i][0]) for i in range(len(seq1))]
    comp = phyChem(seq);
    new = [comp[i][0:25] for i in range(len(comp))]
    entropy  = [];
    entropy.append(SEP_headers[0:25])
    for i in range(1,len(new)):
        seqEntropy = [];
        for j in range(len(new[i])):
            p = new[i][j]; 
            if((1-p) == 0. or p==0.):
                temp = 0;#to store entropy of each sequence
            else:
                temp = -(p*math.log2(p)+(1-p)*math.log2(1-p));
            seqEntropy.append(round(temp,3));
        entropy.append(seqEntropy);
    out = pd.DataFrame(entropy);
    file = open(out123,'w')
    with file:
        writer = csv.writer(file);
        writer.writerows(entropy);
    return entropy;
def shannons_split(file,v,out):
    shannons(file,'tempfile_out')
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),(v*25)):
        pp.append(ss[i:i+(v*25)])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')	
	
####################autocorr####################
def p_aa(prop,a):
    if ((a.upper()=='A') or (a.upper()=='C') or (a.upper()=='D') or (a.upper()=='E') or (a.upper()=='F') or (a.upper()=='G') or (a.upper()=='H') or (a.upper()=='I') or (a.upper()=='K') or (a.upper()=='L') or (a.upper()=='M') or (a.upper()=='N') or (a.upper()=='P') or (a.upper()=='Q') or (a.upper()=='R') or (a.upper()=='S') or (a.upper()=='T') or (a.upper()=='V') or (a.upper()=='W') or (a.upper()=='Y')):
        data=pd.read_table('Data/z_aaindex.csv',sep=',',index_col='INDEX' )
        p=data.loc[prop][a.upper()]
        return p
    else:
        print("Error!: Invalid sequence. Special character(s)/invalid amino acid letter(s) present in input.")
        return
def NMB(prop,seq,d):
    if (d<=30):
        sum=0
        for i in range(len(seq)-d):
            sum=sum+p_aa(prop,seq[i])*p_aa(prop,seq[i+d])
        ats=sum/(len(seq)-d)
        return ats
    else:
        print("Error!: d should be less than equal to 30")
        return
def pav(prop,seq):
    av=0
    for i in range(len(seq)):
        av=av+p_aa(prop,seq[i])
    av=av/len(seq)
    return av
def moran(prop,seq,d):
    if (d<=30):
        s1=0
        s2=0
        p_bar=pav(prop,seq)
        for i in range(len(seq)-d):
            s1=s1+(p_aa(prop,seq[i])-p_bar)*(p_aa(prop,seq[i+d])-p_bar)
        for i in range(len(seq)):
            s2=s2+(p_aa(prop,seq[i])-p_bar)*(p_aa(prop,seq[i])-p_bar)
        return (s1/(len(seq)-d))/(s2/len(seq))
    else:
        print("Error!: d should be less than equal to 30")
        return
def geary(prop,seq,d):
    if (d<=30):
        s1=0
        s2=0
        p_bar=pav(prop,seq)
        for i in range(len(seq)-d):
            s1=s1+(p_aa(prop,seq[i])-p_aa(prop,seq[i+d]))*(p_aa(prop,seq[i])-p_aa(prop,seq[i+d]))
        for i in range(len(seq)):
            s2=s2+(p_aa(prop,seq[i])-p_bar)*(p_aa(prop,seq[i])-p_bar)
        return (s1/(2*(len(seq)-d)))/(s2/(len(seq)-1))
    else:
        print("Error!: d should be less than equal to 30")
        return

def autocorr_full_aa(filename,d,out):
    if (d<=30):
        seq_data=list((pd.read_csv(filename,sep=',',header=None)).iloc[:,0])
        prop=list((pd.read_csv('Data/aaindices.csv',sep=',',header=None)).iloc[0,:])
        output=[[]]
        for k in range(len(prop)):
            output[0]=output[0]+['ACR'+str(d)+'_'+'MB','ACR'+str(d)+'_'+'MO','ACR'+str(d)+'_'+'GE']
        temp=[]
        for i in range(len(seq_data)):
            for j in range(len(prop)):
                temp=temp+[round(NMB(prop[j],seq_data[i],d),3),round(moran(prop[j],seq_data[i],d),3),round(geary(prop[j],seq_data[i],d),3)]
            output.append(temp)
            temp=[]
        file = open(out,'w')
        with file:
            writer = csv.writer(file);
            writer.writerows(output);
        return output
    else:
        print("Error!: d should be less than equal to 30")
        return

def autocorr_split(file,v,lg,out):
    autocorr_full_aa(file,lg,'tempfile_out')
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),(v*3)):
        pp.append(ss[i:i+(v*3)])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')

##################paac####################
def val(AA_1, AA_2, aa, mat):
    return sum([(mat[i][aa[AA_1]] - mat[i][aa[AA_2]]) ** 2 for i in range(len(mat))]) / len(mat)
def paac_1(file,lambdaval,w=0.05):
    data1 = pd.read_csv("Data/data", sep = "\t")
    filename, file_extension = os.path.splitext(file)
    df = pd.read_csv(file, header = None)
    df1 = pd.DataFrame(df[0].str.upper())
    dd = []
    cc = []
    pseudo = []
    aa = {}
    for i in range(len(std)):
        aa[std[i]] = i
    for i in range(0,3):
        mean = sum(data1.iloc[i][1:])/20
        rr = math.sqrt(sum([(p-mean)**2 for p in data1.iloc[i][1:]])/20)
        dd.append([(p-mean)/rr for p in data1.iloc[i][1:]])
        zz = pd.DataFrame(dd)
    head = []
    for n in range(1, lambdaval + 1):
        head.append('_lam' + str(n))
    head = ['PAAC'+str(lambdaval)+sam for sam in head]
    pp = pd.DataFrame()
    ee = []
    for k in range(0,len(df1)):
        cc = []
        pseudo1 = [] 
        for n in range(1,lambdaval+1):
            cc.append(sum([val(df1[0][k][p], df1[0][k][p + n], aa, dd) for p in range(len(df1[0][k]) - n)]) / (len(df1[0][k]) - n))
            qq = pd.DataFrame(cc)
        pseudo = pseudo1 + [(w * p) / (1 + w * sum(cc)) for p in cc]
        ee.append(pseudo)
        ii = round(pd.DataFrame(ee, columns = head),4)
    ii.to_csv(filename+".lam",index = None)
		
def paac(file,lambdaval,out,w=0.05):
    filename, file_extension = os.path.splitext(file)
    paac_1(file,lambdaval,w=0.05)
    aac_comp(file,filename+".aac")
    data1 = pd.read_csv(filename+".aac")
    header = ['PAAC'+str(lambdaval)+'_A','PAAC'+str(lambdaval)+'_C','PAAC'+str(lambdaval)+'_D','PAAC'+str(lambdaval)+'_E','PAAC'+str(lambdaval)+'_F','PAAC'+str(lambdaval)+'_G','PAAC'+str(lambdaval)+'_H','PAAC'+str(lambdaval)+'_I','PAAC'+str(lambdaval)+'_K','PAAC'+str(lambdaval)+'_L','PAAC'+str(lambdaval)+'_M','PAAC'+str(lambdaval)+'_N','PAAC'+str(lambdaval)+'_P','PAAC'+str(lambdaval)+'_Q','PAAC'+str(lambdaval)+'_R','PAAC'+str(lambdaval)+'_S','PAAC'+str(lambdaval)+'_T','PAAC'+str(lambdaval)+'_V','PAAC'+str(lambdaval)+'_W','PAAC'+str(lambdaval)+'_Y','Un']	
    data1.columns = header    
    data2 = pd.read_csv(filename+".lam")
    data3 = pd.concat([data1.iloc[:,:-1],data2], axis = 1).reset_index(drop=True)
    data3.to_csv(out,index=None)
    os.remove(filename+".lam")
    os.remove(filename+".aac")

def paac_split(file,v,lg,out,w):
    paac(file,lg,'tempfile_out',w)
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),int(v*(len(ss)/len(df4)))):
        pp.append(ss[i:i+int((v*(len(ss)/len(df4))))])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')
	
######################apaac############################	
def apaac_1(file,lambdaval,w=0.05):
    data1 = pd.read_csv("Data/data", sep = "\t")
    filename, file_extension = os.path.splitext(file)
    df = pd.read_csv(file, header = None)
    df1 = pd.DataFrame(df[0].str.upper())
    dd = []
    cc = []
    pseudo = []
    aa = {}
    for i in range(len(std)):
        aa[std[i]] = i
    for i in range(0,3):
        mean = sum(data1.iloc[i][1:])/20
        rr = math.sqrt(sum([(p-mean)**2 for p in data1.iloc[i][1:]])/20)
        dd.append([(p-mean)/rr for p in data1.iloc[i][1:]])
        zz = pd.DataFrame(dd)
    head = []
    for n in range(1, lambdaval + 1):
        for e in ('HB','HL','SC'):
            head.append(e+'_lam' + str(n))
    head = ['APAAC'+str(lambdaval)+'_'+sam for sam in head]
    pp = pd.DataFrame()
    ee = []
    for k in range(0,len(df1)):
        cc = [] 
        for n in range(1,lambdaval+1):
            for b in range(0,len(zz)):
                cc.append(sum([zz.loc[b][aa[df1[0][k][p]]] * zz.loc[b][aa[df1[0][k][p + n]]] for p in range(len(df1[0][k]) - n)]) / (len(df1[0][k]) - n))
                qq = pd.DataFrame(cc)
        pseudo = [(w * p) / (1 + w * sum(cc)) for p in cc]
        ee.append(pseudo)
        ii = round(pd.DataFrame(ee, columns = head),4)
    ii.to_csv(filename+".plam",index = None)
		
def apaac(file,lambdaval,out,w=0.05):
    filename, file_extension = os.path.splitext(file)
    apaac_1(file,lambdaval,w=0.05)
    aac_comp(file,filename+".aac")
    data1 = pd.read_csv(filename+".aac")
    headaac = []
    for i in std:
        headaac.append('APAAC'+str(lambdaval)+'_'+i)
    headaac.insert(len(headaac),0)
    data1.columns = headaac
    data2 = pd.read_csv(filename+".plam")
    data3 = pd.concat([data1.iloc[:,:-1],data2], axis = 1).reset_index(drop=True)
    data3.to_csv(out, index = None)
    os.remove(filename+".plam")
    os.remove(filename+".aac")
def apaac_split(file,v,lg,out,w):
    apaac(file,lg,'tempfile_out',w)
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),int(v*(len(ss)/len(df4)))):
        pp.append(ss[i:i+int((v*(len(ss)/len(df4))))])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')
###################################qos#######################################
def qos(file,gap,out,w=0.1):
    ff = []
    filename, file_extension = os.path.splitext(file)
    df = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df[0].str.upper())
    for i in range(0,len(df2)):
        ff.append(len(df2[0][i]))
    if min(ff) < gap:
        print("Error: All sequences' length should be higher than :", gap)
    else:
        mat1 = pd.read_csv("Data/Schneider-Wrede.csv", index_col = 'Name')
        mat2 = pd.read_csv("Data/Grantham.csv", index_col = 'Name')
        s1 = []
        s2 = []
        for i in range(0,len(df2)):
            for n in range(1, gap+1):
                sum1 = 0
                sum2 = 0
                for j in range(0,(len(df2[0][i])-n)):
                    sum1 = sum1 + (mat1[df2[0][i][j]][df2[0][i][j+n]])**2
                    sum2 = sum2 + (mat2[df2[0][i][j]][df2[0][i][j+n]])**2
                s1.append(sum1)
                s2.append(sum2)
        zz = pd.DataFrame(np.array(s1).reshape(len(df2),gap))
        zz["sum"] = zz.sum(axis=1)
        zz2 = pd.DataFrame(np.array(s2).reshape(len(df2),gap))
        zz2["sum"] = zz2.sum(axis=1)
        c1 = []
        c2 = []
        c3 = []
        c4 = []
        h1 = []
        h2 = []
        h3 = []
        h4 = []
        for aa in std:
            h1.append('QSO'+str(gap)+'_SC_' + aa)
        for aa in std:
            h2.append('QSO'+str(gap)+'_G_' + aa)
        for n in range(1, gap+1):
            h3.append('SC' + str(n))
        h3 = ['QSO'+str(gap)+'_'+sam for sam in h3]
        for n in range(1, gap+1):
            h4.append('G' + str(n))
        h4 = ['QSO'+str(gap)+'_'+sam for sam in h4]
        for i in range(0,len(df2)):
            AA = {}
            for j in std:
                AA[j] = df2[0][i].count(j)
                c1.append(AA[j] / (1 + w * zz['sum'][i]))
                c2.append(AA[j] / (1 + w * zz2['sum'][i]))
            for k in range(0,gap):
                c3.append((w * zz[k][i]) / (1 + w * zz['sum'][i]))
                c4.append((w * zz[k][i]) / (1 + w * zz['sum'][i]))
        pp1 = np.array(c1).reshape(len(df2),len(std))
        pp2 = np.array(c2).reshape(len(df2),len(std))
        pp3 = np.array(c3).reshape(len(df2),gap)
        pp4 = np.array(c4).reshape(len(df2),gap)
        zz5 = round(pd.concat([pd.DataFrame(pp1, columns = h1),pd.DataFrame(pp2,columns = h2),pd.DataFrame(pp3, columns = h3),pd.DataFrame(pp4, columns = h4)], axis = 1),4)
        zz5.to_csv(out, index = None, encoding = 'utf-8')	
def qos_split(file,v,lg,out,w):
    qos(file,lg,'tempfile_out',w)
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),int(v*(len(ss)/len(df4)))):
        pp.append(ss[i:i+int((v*(len(ss)/len(df4))))])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')
##########################soc################
def soc(file,gap,out):
    ff = []
    filename, file_extension = os.path.splitext(file)
    df = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df[0].str.upper())
    for i in range(0,len(df2)):
        ff.append(len(df2[0][i]))
    if min(ff) < gap:
        print("Error: All sequences' length should be higher than :", gap)
        return 0
    mat1 = pd.read_csv("Data/Schneider-Wrede.csv", index_col = 'Name')
    mat2 = pd.read_csv("Data/Grantham.csv", index_col = 'Name')
    h1 = []
    h2 = []
    for n in range(1, gap+1):
        h1.append('SC' + str(n))
    for n in range(1, gap + 1):
        h2.append('G' + str(n))
    h1 = ['SOC'+str(gap)+'_'+sam for sam in h1]
    h2 = ['SOC'+str(gap)+'_'+sam for sam in h2]
    s1 = []
    s2 = []
    for i in range(0,len(df2)):
        for n in range(1, gap+1):
            sum = 0
            sum1 =0
            sum2 =0
            sum3 =0
            for j in range(0,(len(df2[0][i])-n)):
                sum = sum + (mat1[df2[0][i][j]][df2[0][i][j+n]])**2
                sum1 = sum/(len(df2[0][i])-n)
                sum2 = sum2 + (mat2[df2[0][i][j]][df2[0][i][j+n]])**2
                sum3 = sum2/(len(df2[0][i])-n)
            s1.append(sum1)
            s2.append(sum3)
    zz = np.array(s1).reshape(len(df2),gap)
    zz2 = np.array(s2).reshape(len(df2),gap)
    zz3 = round(pd.concat([pd.DataFrame(zz, columns = h1),pd.DataFrame(zz2,columns = h2)], axis = 1),4)
    zz3.to_csv(out, index = None, encoding = 'utf-8') 

def soc_split(file,v,lg,out):
    soc(file,lg,'tempfile_out')
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),int(v*(len(ss)/len(df4)))):
        pp.append(ss[i:i+int((v*(len(ss)/len(df4))))])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')
##########################################CTC###################################
x = [1, 2, 3, 4, 5, 6,7]
p=[]
Y=[]
LS=[]


for i in range(len(x)):
    p=itertools.product(x,repeat=3)
    p=list(p)

def concatenate_list_data(list):
    result= ''
    for element in list:
        result += str(element)
    return result

for i in range(len(p)):
    LS.append(concatenate_list_data(p[i]))

def repstring(string):
    string=string.upper()
    char={"A":"1","G":"1","V":"1","I":"2","L":"2","F":"2","P":"2","Y":"3","M":"3","T":"3","S":"3","H":"4","N":"4","Q":"4","W":"4","R":"5","K":"5","D":"6","E":"6","C":"7"}
    string=list(string)
    for index,item in enumerate(string):
        for key,value in char.items():
            if item==key:
                string[index]=value
    return("".join(string))

def occurrences(string, sub_string):
    count=0
    beg=0
    while(string.find(sub_string,beg)!=-1) :
        count=count+1
        beg=string.find(sub_string,beg)
        beg=beg+1
    return count


def CTC(filename,out):
    df = pd.DataFrame(columns=['Sequence','Triad:Frequency'])
    data=list((pd.read_csv(filename,sep=',',header=None)).iloc[:,0])
    for i in range(len(data)):
        data1=''
        data1=str(data[i])
        data1=data1.upper()
        allowed = set(('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'))
        is_data_invalid = set(data1).issubset(allowed)
        if is_data_invalid==False:
            print("Errror: Please check for invalid inputs in the sequence.","\nError in: ","Sequence number=",i+1,",","Sequence = ",data[i],",","\nNOTE: Spaces, Special characters('[@_!#$%^&*()<>?/\|}{~:]') and Extra characters(BJOUXZ) should not be there.")
            return
        df.at[i,'Sequence'] = data[i]
        Y.append("".join(repstring(str(data[i]))))
    val2=[[]]
    for f in range(len(LS)):
        val2[0]=val2[0]+["CTC_"+str(LS[f])]
    for j in range(len(data)):
        MM=[]
        for m in range(len(LS)):
            MM=MM+[occurrences(Y[j],LS[m])]
        Min_MM=min(MM)
        Max_MM=max(MM)
        if (Max_MM==0):
            print("Errror: Splits/ Sequence length should be greater than equal to 3")
            return
        val=[]
#         val.append(data[j])
        for k in range(len(LS)):
            val=val+[round(((occurrences(Y[j],LS[k])-Min_MM)/Max_MM),3)]
        val2.append(val)
#     print(val2)
    #file= open(sys.argv[2],'w', newline='')#output file
    file= open(out,'w', newline='')
    with file:
        writer=csv.writer(file);
        writer.writerows(val2);
    return val2

def ctc_split(file,v,out):
    CTC(file,'tempfile_out')
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),int(v*(len(ss)/len(df4)))):
        pp.append(ss[i:i+int((v*(len(ss)/len(df4))))])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')

######################################AAC###################################
def ctd(file,out):
    attr=pd.read_csv("Data/aa_attr_group.csv", sep="\t")
    filename, file_extension = os.path.splitext(file)
    df1 = pd.read_csv(file, header = None)
    df = pd.DataFrame(df1[0].str.upper())
    n = 0
    stt1 = []
    m = 1
    for i in range(0,len(attr)) :
        st =[]
        stt1.append([])
        for j in range(0,len(df)) :
            stt1[i].append([])
            for k in range(0,len(df[0][j])) :
                while m < 4 :
                    while n < len(attr.iloc[i,m]) :
                        if df[0][j][k] == attr.iloc[i,m][n] :
                            st.append(m)
                            stt1[i][j].append(m)
                        n += 2
                    n = 0
                    m += 1
                m = 1
#####################Composition######################
    f = open("compout_1", 'w')
    sys.stdout = f
    std = [1,2,3]
    print("1,2,3,")
    for p in range (0,len(df)) :
        for ii in range(0,len(stt1)) :
            #for jj in stt1[ii][p]:
            for pp in std :
                count = 0
                for kk in stt1[ii][p] :
                    temp1 = kk
                    if temp1 == pp :
                        count += 1
                    composition = (count/len(stt1[ii][p]))*100
                print("%.2f"%composition, end = ",")
            print("")
    f.truncate()

#################################Transition#############
    tt = []
    tr=[]
    kk =0
    for ii in range(0,len(stt1)) :
        tt = []
        tr.append([])
        for p in range (0,len(df)) :
            tr[ii].append([])
            while kk < len(stt1[ii][p]) :
                if kk+1 <len(stt1[ii][p]):
                #if  stt1[ii][p][kk] < stt1[ii][p][kk+1] or stt1[ii][p][kk] > stt1[ii][p][kk+1]: # condition for adjacent values
                    tt.append(stt1[ii][p][kk])
                    tt.append(stt1[ii][p][kk+1])
                    tr[ii][p].append(stt1[ii][p][kk])
                    tr[ii][p].append(stt1[ii][p][kk+1])

                kk += 1
            kk = 0

    pp = 0
    xx = []
    xxx = []
    for mm in range(0,len(tr)) :
        xx = []
        xxx.append([])
        for nn in range(0,len(tr[mm])):
            xxx[mm].append([])
            while pp < len(tr[mm][nn]) :
                xx .append(tr[mm][nn][pp:pp+2])
                xxx[mm][nn].append(tr[mm][nn][pp:pp+2])
                pp+=2
            pp = 0

    f1 = open("compout_2", 'w')
    sys.stdout = f1
    std1 = [[1,1],[1,2],[1,3],[2,1],[2,2],[2,3],[3,1],[3,2],[3,3]]
    print("1->1,1->2,1->3,2->1,2->2,2->3,3->1,3->2,3->3,")
    for rr in range(0,len(df)) :
        for qq in range(0,len(xxx)):
            for tt in std1 :
                count = 0
                for ss in xxx[qq][rr] :
                    temp2 = ss
                    if temp2 == tt :
                        count += 1
                print(count, end = ",")
            print("")
    f1.truncate()

    #################################Distribution#############
    c_11 = []
    c_22 = []
    c_33 = []
    zz = []
    #print("0% 25% 50% 75% 100%")
    for x in range(0,len(stt1)) :
        #c_11.append([])
        c_22.append([])
        #c_33.append([])
        yy_c_1 = []
        yy_c_2 = []
        yy_c_3 = []
        ccc = []

        k = 0
        j = 0
        for y in range(0,len(stt1[x])):
            #c_11[x].append([])
            c_22[x].append([])
            for i in range(1,4) :
                cc = []
                c1 = [index for index,value in enumerate(stt1[x][y]) if value == i]
                c_22[x][y].append(c1)
    cc = []
    for ss in range(0,len(df)):
        for uu in range(0,len(c_22)):
            for mm in range(0,3):
                for ee in range(0,101,25):
                    k = (ee*(len(c_22[uu][ss][mm])))/100
                    cc.append(math.floor(k))
    f2 = open('compout_3', 'w')
    sys.stdout = f2
    print("0% 25% 50% 75% 100%")
    for i in range (0,len(cc),5):
        print(*cc[i:i+5])
    f2.truncate()
    head = []
    header1 = ['CeTD_HB','CeTD_VW','CeTD_PO','CeTD_PZ','CeTD_CH','CeTD_SS','CeTD_SA']
    for i in header1:
        for j in range(1,4):
            head.append(i+str(j))
    df11 = pd.read_csv("compout_1")
    df_1 = df11.iloc[:,:-1]
    zz = pd.DataFrame()
    for i in range(0,len(df_1),7):
        zz = zz.append(pd.DataFrame(pd.concat([df_1.loc[i],df_1.loc[i+1],df_1.loc[i+2],df_1.loc[i+3],df_1.loc[i+4],df_1.loc[i+5],df_1.loc[i+6]],axis=0)).transpose()).reset_index(drop=True)
    zz.columns = head
    #zz.to_csv(filename+".ctd_comp", index=None, encoding='utf-8')
    head2 = []
    header2 = ['CeTD_11','CeTD_12','CeTD_13','CeTD_21','CeTD_22','CeTD_23','CeTD_31','CeTD_32','CeTD_33']
    for i in header2:
        for j in ('HB','VW','PO','PZ','CH','SS','SA'):
            head2.append(i+'_'+str(j))
    df12 = pd.read_csv("compout_2")
    df_2 = df12.iloc[:,:-1]
    ss = pd.DataFrame()
    for i in range(0,len(df_2),7):
        ss = ss.append(pd.DataFrame(pd.concat([df_2.loc[i],df_2.loc[i+1],df_2.loc[i+2],df_2.loc[i+3],df_2.loc[i+4],df_2.loc[i+5],df_2.loc[i+6]],axis=0)).transpose()).reset_index(drop=True)
    ss.columns = head2
    #ss.to_csv(filename+".ctd_trans", index=None, encoding='utf-8')
    head3 = []
    header3 = ['CeTD_0_p','CeTD_25_p','CeTD_50_p','CeTD_75_p','CeTD_100_p']
    header4 = ['HB','VW','PO','PZ','CH','SS','SA']
    for j in range(1,4):
        for k in header4:
            for i in header3:
                head3.append(i+'_'+k+str(j))
    df_3 = pd.read_csv("compout_3", sep=" ")
    rr = pd.DataFrame()
    for i in range(0,len(df_3),21):
        rr = rr.append(pd.DataFrame(pd.concat([df_3.loc[i],df_3.loc[i+1],df_3.loc[i+2],df_3.loc[i+3],df_3.loc[i+4],df_3.loc[i+5],df_3.loc[i+6],df_3.loc[i+7],df_3.loc[i+8],df_3.loc[i+9],df_3.loc[i+10],df_3.loc[i+11],df_3.loc[i+12],df_3.loc[i+13],df_3.loc[i+14],df_3.loc[i+15],df_3.loc[i+16],df_3.loc[i+17],df_3.loc[i+18],df_3.loc[i+19],df_3.loc[i+20]],axis=0)).transpose()).reset_index(drop=True)
    rr.columns = head3
    cotrdi= pd.concat([zz,ss,rr],axis=1)
    cotrdi.to_csv(out, index=None, encoding='utf-8')
    os.remove('compout_1')
    os.remove('compout_2')
    os.remove('compout_3')
    #rr.to_csv(filename+".ctd_dist", index=None, encoding='utf-8')

def ctd_split(file,v,out):
    ctd(file,'tempfile_out')
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),int(v*(len(ss)/len(df4)))):
        pp.append(ss[i:i+int((v*(len(ss)/len(df4))))])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')
################################################################################
class ProgressBar(object):
    DEFAULT = 'Progress: %(bar)s %(percent)3d%%'
    FULL = "Status: %(bar)s %(percent)3d%% Completed...."

    def __init__(self, total, width=40, fmt=DEFAULT, symbol='=',
                 output=sys.stderr):
        assert len(symbol) == 1

        self.total = total
        self.width = width
        self.symbol = symbol
        self.output = output
        self.fmt = re.sub(r'(?P<name>%\(.+?\))d',
            r'\g<name>%dd' % len(str(total)), fmt)

        self.current = 0

    def __call__(self):
        percent = self.current / float(self.total)
        size = int(self.width * percent)
        remaining = self.total - self.current
        bar = '[' + self.symbol * size + ' ' * (self.width - size) + ']'

        args = {
            'total': self.total,
            'bar': bar,
            'current': self.current,
            'percent': percent * 100,
            'remaining': remaining
        }
        print('\r' + self.fmt % args, file=self.output, end='')

    def done(self):
        self.current = self.total
        self()
        print('', file=self.output)
############################################################################	
def searchAAIndex(AAIndex):
    found = -1;
    for i in range(len(AAIndexNames)):
        if(str(AAIndex) == AAIndexNames.iloc[i][0]):
            found = i;
    return found;
	
def phychem_AAI(file,AAIn,mode): 
    if(type(file) == str):
        seq = pd.read_csv(file,header=None, sep=',');
        seq=seq.T
        seq[0].values.tolist()
        seq=seq[0];
    else:
        seq = file;
    if(type(AAIn) == str):
        AAI = pd.read_csv(AAIn,header=None, sep=',');
        print (AAI.shape)
        AAI = AAI.values.tolist()
        print(AAI)
    else:
        AAI = AAIn;
    l2 = len(AAI)
    header  = AAI[0:l2];
    final=[];
    final.append(AAI);
    l1 = len(seq);
    seq=[seq[i].upper() for i in range(l1)];
    for i in range(l1):
        coded = encode(seq[i]);
        temp=[];
        for j in range(l2):
            pos = searchAAIndex(AAI[j]);
            sum=0;
            for k in range(len(coded)):
                val = AAIndex.iloc[pos,int(coded[k])]
                sum=sum+val;
            avg = round(sum/len(seq[i]),3);
            temp.append(avg);
        final.append(temp);
    if mode == 'all' :        
        file = open('AAIndex_all','w')        
    if mode == 'NT' :        
        file = open('AAIndex_NT','w')
    if mode == 'CT' :        
        file = open('AAIndex_CT','w')	
    if mode == 'rest' :        
        file = open('AAIndex_rest','w')		
    with file:
        writer = csv.writer(file);
        writer.writerows(final);
    return final;
def AAIndex_Phychem(filename,mode='all',m=0,n=0):
    time.sleep(1.5)
    seq = [];    
    if(type(filename) == str):
        seq1 = pd.read_csv(filename,header=None, sep=',');
        seq1 = pd.DataFrame(seq1[0].str.upper());
        [seq.append(seq1.iloc[i][0]) for i in range(len(seq1))];
    else:
        seq  = filename;
    AAIn = AAindices;
    if(type(AAIn) == str):
        AAI = pd.read_csv(AAIn,header=None, sep=',');
        l2 = AAI.shape[1]
        AAI = AAI.values.tolist()
        AAI  = list(itertools.chain.from_iterable(AAI))
    else:
        AAI = AAIn;
    l2 = len(AAI)
    l=len(seq)
    newseq=[];
    for i in range(0,l):
            l = len(seq[i]);
            if(mode=='NT'):
                n=m;
                print('Inside NT, m=',m,'n=',n)
                if(n!=0):
                    new = seq[i][0:n];
                    newseq.append(new);
                elif(n>l):
                    print('Warning! Sequence',i,"'s size is less than n. The output table would have NaN for this sequence");
                else:
                    print('Value of n is mandatory, it cannot be 0')
                    break;
            elif(mode=='CT'):
                n=m;
                if(n!=0):
                    new = seq[i][(len(seq[i])-n):]
                    newseq.append(new);
                elif(n>l):
                    print('WARNING: Sequence',i+1,"'s size is less than the value of n given. The output table would have NaN for this sequence");
                else:
                    print('Value of n is mandatory, it cannot be 0')
                    break;
            elif(mode=='all'):
                newseq = seq;
            elif(mode=='rest'):
                if(m==0):
                    print('Kindly provide start index for rest, it cannot be 0');
                    break;
                else:
                    #if(n<=len(seq[i])):
                    new = seq[i][m:len(seq[i])-n]
                    newseq.append(new)
            else:
                print("Wrong Mode. Enter 'NT', 'CT','all' or 'rest'");
    if(mode=='split'):
        newseq  = list(itertools.chain.from_iterable(newseq));
    phychem_AAI(newseq,AAI,mode);

def AAI_split(file,v,out):
    file1 = split(file,v)
    file2 = file1.iloc[:,0]
    AAIndex_Phychem(file2,'all')
    df4 = pd.read_csv('AAIndex_all')
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),int(v*(len(ss)/len(df4)))):
        pp.append(ss[i:i+int((v*(len(ss)/len(df4))))])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5.columns = 'AAI_'+df5.columns
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('AAIndex_all')

##################Parts of Seqeunces#########################
def nt(file,n):
    filename, file_extension = os.path.splitext(file)
    df1 = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df1[0].str.upper())
    df3 = []
    for i in range(0,len(df2)):
        df3.append(df2[0][i][0:n])
        df4 = pd.DataFrame(df3)
        #df4.to_csv(filename+".nt", index = None, header = False, encoding = 'utf-8')
    for i in range(0,len(df4)):
        ss = len(df4[0][i])
        if ss < n:
            print('\nSequence number',i+1,'has length of',ss,'which is less than the provided value of N-terminal, that is',n,'. Kindly check the sequences.')
            os.remove(file_output)
            sys.exit()
    return df4
def ct(file,n):
    filename, file_extension = os.path.splitext(file)
    df1 = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df1[0].str.upper())
    df3 = []
    for i in range(0,len(df2)):
        df3.append(df2[0][i][-n:])
        df4 = pd.DataFrame(df3)
        #df4.to_csv(filename+".ct", index = None, header = False, encoding = 'utf-8')
    for i in range(0,len(df4)):
        ss = len(df4[0][i])
        if ss < n:
            print('\nSequence number',i+1,'has length of',ss,'which is less than the provided value of C-terminal, that is',n,'. Kindly check the sequences.')
            os.remove(file_output)
            sys.exit()
    return df4
def rest(file,n,c):
    filename, file_extension = os.path.splitext(file)
    df1 = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df1[0].str.upper())
    df3 = []
    for i in range(0,len(df2)):
        if c == 0:
            df3.append(df2[0][i][n:])
        else:
            df3.append(df2[0][i][n:-c])
    df4 = pd.DataFrame(df3)
    for i in range(0,len(df4)):
        ss = len(df4[0][i])
        if ss == 0:
            print('\nSequence number',i+1,'has length of',ss,'after removing provided N- and C-terminal residues, that is',str(n)+','+str(c),'. Kindly check the sequences.')
            os.remove(file_output)
            sys.exit()
        #df4.to_csv(filename+".rest", index = None, header = False, encoding = 'utf-8')
    return df4
def restnc(file,n):
    df1 = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df1[0].str.upper())
    df3 = []
    for i in range(0,len(df2)):
        if n == 0:
            df3.append(df2[0][i][n:])
        else:
            df3.append(df2[0][i][n:-n])
    df4 = pd.DataFrame(df3)
    df4 = pd.DataFrame(df3)
    for i in range(0,len(df4)):
        ss = len(df4[0][i])
        if ss == 0:
            print('\nSequence number',i+1,'has length of',ss,'after removing provided',str(n),' residues from N- and C-terminal. Kindly check the sequences.')
            os.remove(file_output)
            sys.exit()
        #df4.to_csv(filename+".rest", index = None, header = False, encoding = 'utf-8')
    return df4
def nct(file,n):
    filename, file_extension = os.path.splitext(file)
    df1 = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df1[0].str.upper())
    df3 = []
    for i in range(0,len(df2)):
        df3.append(df2[0][i][:n]+df2[0][i][::-1][:n])
        df4 = pd.DataFrame(df3)
        #df4.to_csv(filename+".ct", index = None, header = False, encoding = 'utf-8')
    for i in range(0,len(df4)):
        ss = len(df4[0][i])
        if ss/2 < n:
            print('\nSequence number',i+1,'has length of',int(ss/2),'which is less than the provided value of N- and C-terminal, that is',n,'. Kindly check the sequences.')
            os.remove(file_output)
            sys.exit()
    return df4
def split(file,v):
    filename, file_extension = os.path.splitext(file)
    df1 = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df1[0].str.upper())
    for i in range(0,len(df2)):
        ss = len(df2[0][i])
        if ss == 0:
            print('\nSequence number',i+1,'has length of',ss,'. Hence, the number of splits should be between 2 to',ss,'. Kindly provide number of splits in the suggested range.')
            os.remove(file_output)
            sys.exit()
    k1 = []
    for e in range(0,len(df2)):
        s = 0
        k2 = []
        r = 0
        if len(df2[0][e])%v == 0:
            k2.extend(repeat(int(len(df2[0][e])/v),v))
        else:
            r = int(len(df2[0][e])%v)
            k2.extend(repeat(int(len(df2[0][e])/v),v-1))
            k2.append((int(len(df2[0][e])/v))+r)
        for j in k2:
            df3 = df2[0][e][s:j+s]
            k1.append(df3)
            s = j+s
    df4 = pd.DataFrame(k1)
    #df4.to_csv(filename+".split", index = None, header = False, encoding = 'utf-8')
    return df4
####################################File reading###################################
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
##############################################################################
#Seqeunce file
Sequence_1 = args.input
file_output = str(uuid.uuid4())
readseq(Sequence_1,file_output) 
Sequence = file_output
seq = file_output
#Function to be called
if args.job == None:
    Job = 'AAC'
else:
    Job = args.job
#Number of residue from N-terminal
if args.n_terminal == None:
    nter = int(0)
else:
    nter = int(args.n_terminal)
#Number of residues from C-terminal
if args.c_terminal == None:
    cter = int(0)
else:
    cter = int(args.c_terminal)
#Number of residues from C-terminal
if args.nc_terminal == None:
    ncter = int(0)
else:
    ncter = int(args.nc_terminal)
#Number of residues removed from N-terminal
if args.rest_n == None:
    nrest = int(0)
else:
    nrest=int(args.rest_n)
#Number of residues removed from C-terminal
if args.rest_c == None:
    crest = int(0)
else:
    crest=int(args.rest_c)
#Number of residues after removing residues from N- and  C-terminal
if args.rest_nc == None:
    ncrest = int(0)
else:
    ncrest=int(args.rest_nc)
#Number of splits a sequence is divided into
if args.split == None:
    sp = int(0)
else:
    sp = int(args.split)
#Lambda/gap/lag
if args.lag == None:
    lg = int(1)
else:
    lg = int(args.lag)
#weighting factor for QSO
if args.weight == None:
    wq = float(0.1)
else:
    wq = float(args.weight)
#weighting factor for PAAC
if args.pweight == None:
    pw = float(0.05)
else:
    pw = float(args.pweight)
#Output_file
if args.output == None:
    result_filename= "pfeature_result.csv"
else:
    result_filename = args.output

if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
    print('Summary of Parameters:''\n')
    print('Input File:',Sequence_1,'; Job:', Job,'; Output File:',result_filename,)
if nter != 0:
    print('Summary of Parameters:''\n')
    print('Input File:',Sequence_1,'; Job:', Job,'; N-terminal:',nter,'; Output File:',result_filename,)
if cter != 0:
    print('Summary of Parameters:''\n')
    print('Input File:',Sequence_1,'; Job:', Job,'; C-terminal:',cter,'; Output File:',result_filename,)
if ncter != 0:
    print('Summary of Parameters:''\n')
    print('Input File:',Sequence_1,'; Job:', Job,'; NC-terminal:',ncter,'; Output File:',result_filename,)
if nrest != 0 or crest != 0:
    print('Summary of Parameters:''\n')
    print('Input File:',Sequence_1,'; Job:', Job,'; Rest:',(nrest,crest),'; Output File:',result_filename,)
if ncrest != 0 :
    print('Summary of Parameters:''\n')
    print('Input File:',Sequence_1,'; Job:', Job,'; NCRest:',ncrest,'; Output File:',result_filename,)
if sp != 0:
    print('Summary of Parameters:''\n')
    print('Input File:',Sequence_1,'; Job:', Job,'; Number of Split:',sp,'; Output File:',result_filename,)
if Job == 'AAC' or Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        aac_comp(seq,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.aac',index=None)
        os.remove('tempfile_out')
    if nter != 0:
        file1 = nt(Sequence,nter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        aac_comp('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.aac_nt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        file1 = ct(Sequence,cter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        aac_comp('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.aac_ct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        file1 = rest(Sequence,nrest,crest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        aac_comp('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.aac_rt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0 :
        file1 = restnc(Sequence,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        aac_comp('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.aac_rnc',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        file1 = nct(Sequence,ncter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        aac_comp('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'NC'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.aac_nct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        aac_split(seq,sp,result_filename)
        aac_split(seq,sp,'sam_allcomp.aac_st')
if Job == 'DPC' or Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        dpc_comp(seq,lg,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        print('Order of Dipeptide:',lg,'\n')
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.dpc',index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.dpc',index=None)
        os.remove('tempfile_out')
    if nter != 0:
        file1 = nt(Sequence,nter)
        print('Order of Dipeptide:',lg,'\n')
        file1.to_csv('sam_input.csv', index=None, header=False)
        dpc_comp('sam_input.csv',lg,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.dpc_nt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        file1 = ct(Sequence,cter)
        print('Order of Dipeptide:',lg,'\n')
        file1.to_csv('sam_input.csv', index=None, header=False)
        dpc_comp('sam_input.csv',lg,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.dpc_ct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        file1 = rest(Sequence,nrest,crest)
        print('Order of Dipeptide:',lg,'\n')
        file1.to_csv('sam_input.csv', index=None, header=False)
        dpc_comp('sam_input.csv',lg,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.dpc_rt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0:
        file1 = restnc(Sequence,ncrest)
        print('Order of Dipeptide:',lg,'\n')
        file1.to_csv('sam_input.csv', index=None, header=False)
        dpc_comp('sam_input.csv',lg,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.dpc_rnc',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        file1 = nct(Sequence,ncter)
        print('Order of Dipeptide:',lg,'\n')
        file1.to_csv('sam_input.csv', index=None, header=False)
        dpc_comp('sam_input.csv',lg,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'NC'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.dpc_nct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        print('Order of Dipeptide:',lg,'\n')
        dpc_split(seq,lg,sp,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.dpc_st',index=None)
        os.remove('tempfile_out')
if Job == 'TPC' or Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        tpc_comp(seq,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.tpc',index=None)
        os.remove('tempfile_out')
    if nter != 0:
        file1 = nt(Sequence,nter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        tpc_comp('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.tpc_nt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        file1 = ct(Sequence,cter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        tpc_comp('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.tpc_ct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        file1 = rest(Sequence,nrest,crest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        tpc_comp('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.tpc_rt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0:
        file1 = restnc(Sequence,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        tpc_comp('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.tpc_rnc',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        file1 = nct(Sequence,ncter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        tpc_comp('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'NC'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.tpc_nct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        tpc_split(seq,sp,result_filename)
        tpc_split(seq,sp,'sam_allcomp.tpc_st')
if Job == 'ATC' or Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        atc(seq,result_filename)
        atc(seq,'sam_allcomp.atc')
    if nter != 0:
        file1 = nt(Sequence,nter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        atc('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.atc_nt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        file1 = ct(Sequence,cter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        atc('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.atc_ct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        file1 = rest(Sequence,nrest,crest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        atc('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.atc_rt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0:
        file1 = restnc(Sequence,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        atc('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.atc_rnc',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        file1 = nct(Sequence,ncter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        atc('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'NC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.atc_nct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        atc_split(seq,sp,result_filename)
        atc_split(seq,sp,'sam_allcomp.atc_st')
if Job == 'BTC' or Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        bond(seq,result_filename)
        bond(seq,'sam_allcomp.btc')
    if nter != 0:
        file1 = nt(Sequence,nter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        bond('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.btc_nt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        file1 = ct(Sequence,cter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        bond('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.btc_ct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        file1 = rest(Sequence,nrest,crest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        bond('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.btc_rt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0:
        file1 = restnc(Sequence,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        bond('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.btc_rnc',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        file1 = nct(Sequence,ncter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        bond('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'NC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.btc_nct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        bond_split(seq,sp,result_filename)
        bond_split(seq,sp,'sam_allcomp.btc_st')
if Job == 'PCP' or Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        pcp_1(seq,result_filename)
        pcp_1(seq,'sam_allcomp.pcp')
    if nter != 0:
        file1 = nt(Sequence,nter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        pcp_1('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.pcp_nt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        file1 = ct(Sequence,cter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        pcp_1('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.pcp_ct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        file1 = rest(Sequence,nrest,crest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        pcp_1('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.pcp_rt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0:
        file1 = restnc(Sequence,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        pcp_1('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.pcp_rnc',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        file1 = nct(Sequence,ncter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        pcp_1('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'NC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.pcp_nct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        pcp_split(seq,sp,result_filename)
        pcp_split(seq,sp,'sam_allcomp.pcp_st')
if Job == 'AAI':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        AAIndex_Phychem(seq,'all')
        df = pd.read_csv('AAIndex_all')
        df.columns = 'AAI_'+df.columns
        df.to_csv(result_filename,index=None)
        os.remove('AAIndex_all')
    if nter != 0:
        AAIndex_Phychem(seq,'NT',0,nter)
        df = pd.read_csv('AAIndex_NT')
        df.columns = 'NAAI_'+df.columns
        df.to_csv(result_filename,index=None)
        os.remove('AAIndex_NT')
    if cter != 0:
        AAIndex_Phychem(seq,'CT',0,cter)
        df = pd.read_csv('AAIndex_CT')
        df.columns = 'CAAI_'+df.columns
        df.to_csv(result_filename,index=None)
        os.remove('AAIndex_CT')
    if nrest != 0 or crest != 0:
        AAIndex_Phychem(seq,'rest',nrest,crest)
        df = pd.read_csv('AAIndex_rest')
        df.columns = 'RAAI_'+df.columns
        df.to_csv(result_filename,index=None)
        os.remove('AAIndex_rest')
    if ncrest != 0:
        file1 = restnc(Sequence,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        AAIndex_Phychem('sam_input.csv','all')
        df = pd.read_csv('AAIndex_all')
        df.columns = 'RNCAAI_'+df.columns
        df.to_csv(result_filename,index=None)
        os.remove('sam_input.csv')
        os.remove('AAIndex_all')
    if ncter != 0:
        file1 = nct(Sequence,ncter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        AAIndex_Phychem('sam_input.csv','all')
        df = pd.read_csv('AAIndex_all')
        df.columns = 'NCAAI_'+df.columns
        df.to_csv(result_filename,index=None)
        os.remove('sam_input.csv')
        os.remove('AAIndex_all')
    if sp != 0:
        AAI_split(seq,sp,result_filename)
if Job=='RRI' or Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        RAAC(seq,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.rri',index=None)
        os.remove('tempfile_out')
    if nter != 0:
        file1 = nt(seq,nter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        RAAC('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.rri_nt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        file1 = ct(seq,cter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        RAAC('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.rri_ct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        file1 = rest(seq,nrest,crest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        RAAC('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.rri_rt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0:
        file1 = restnc(seq,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        RAAC('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.rri_rnc',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        file1 = nct(Sequence,ncter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        RAAC('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'NC'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.rri_nct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        RAAC_split(seq,sp,result_filename)
        RAAC_split(seq,sp,'sam_allcomp.rri_st')
if Job=='PRI' or Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        repeats(seq,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.pri',index=None)
        os.remove('tempfile_out')
    if nter != 0:
        file1 = nt(seq,nter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        repeats('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.pri_nt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        file1 = ct(seq,cter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        repeats('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.pri_ct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        file1 = rest(seq,nrest,crest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        repeats('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.pri_rt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0:
        file1 = restnc(seq,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        repeats('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.pri_rnc',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        file1 = nct(Sequence,ncter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        repeats('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.pri_nct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        repeats_split(seq,sp,result_filename)
        repeats_split(seq,sp,'sam_allcomp.pri_st')
if Job=='DDR' or Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        DDOR(seq,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.ddr',index=None)
        os.remove('tempfile_out')
    if nter != 0:
        file1 = nt(seq,nter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        DDOR('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.ddr_nt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        file1 = ct(seq,cter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        DDOR('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.ddr_ct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        file1 = rest(seq,nrest,crest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        DDOR('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.ddr_rt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0:
        file1 = restnc(seq,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        DDOR('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.ddr_rnc',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        file1 = nct(Sequence,ncter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        DDOR('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'NC'+df.columns
        df.iloc[:,:-1].to_csv(result_filename,index=None)
        df.iloc[:,:-1].to_csv('sam_allcomp.ddr_nct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        DDOR_split(seq,sp,result_filename)
        DDOR_split(seq,sp,'sam_allcomp.ddr_st')
if Job == 'SEP' or Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        SE(seq,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.sep',index=None)
        os.remove('tempfile_out')
    if nter != 0:
        file1 = nt(Sequence,nter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        SE('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.sep_nt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        file1 = ct(Sequence,cter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        SE('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.sep_ct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        file1 = rest(Sequence,nrest,crest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        SE('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.sep_rt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0:
        file1 = restnc(Sequence,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        SE('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.sep_rnc',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        file1 = nct(Sequence,ncter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        SE('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'NC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.sep_nct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        file1 =split(seq,sp)
        file1.to_csv('sam_input.csv', index=None, header=False)
        SE_split('sam_input.csv',sp,result_filename)
        SE_split('sam_input.csv',sp,'sam_allcomp.sep_st')
if Job == 'SER' or Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        SE_residue_level(seq,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.ser',index=None)
        os.remove('tempfile_out')
    if nter != 0:
        file1 = nt(Sequence,nter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        SE_residue_level('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.ser_nt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        file1 = ct(Sequence,cter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        SE_residue_level('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.ser_ct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        file1 = rest(Sequence,nrest,crest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        SE_residue_level('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.ser_rt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0:
        file1 = restnc(Sequence,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        SE_residue_level('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.ser_rnc',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        file1 = nct(Sequence,ncter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        SE_residue_level('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'NC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.ser_nct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        file1 =split(seq,sp)
        file1.to_csv('sam_input.csv', index=None, header=False)
        SE_residue_level_split('sam_input.csv',sp,result_filename)
        SE_residue_level_split('sam_input.csv',sp,'sam_allcomp.ser_st')
if Job == 'SPC' or Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        shannons(seq,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df = df.round(3)
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.spc',index=None)
        os.remove('tempfile_out')
    if nter != 0:
        file1 = nt(Sequence,nter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        shannons('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.spc_nt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        file1 = ct(Sequence,cter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        shannons('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.spc_ct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        file1 = rest(Sequence,nrest,crest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        shannons('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.spc_rt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0:
        file1 = restnc(Sequence,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        shannons('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.spc_rnc',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        file1 = nct(Sequence,ncter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        shannons('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'NC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.spc_nct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        file1 =split(seq,sp)
        file1.to_csv('sam_input.csv', index=None, header=False)
        shannons_split('sam_input.csv',sp,result_filename)
        shannons_split('sam_input.csv',sp,'sam_allcomp.spc_st')
        os.remove('sam_input.csv')
if Job == 'ACR':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        print('Value of lag:',lg,'\n')
        autocorr_full_aa(seq,lg,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.to_csv(result_filename,index=None)
        os.remove('tempfile_out')
    if nter != 0:
        print('Value of lag:',lg,'\n')
        file1 = nt(Sequence,nter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        autocorr_full_aa('sam_input.csv',lg,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.to_csv(result_filename,index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        print('Value of lag:',lg,'\n')
        file1 = ct(Sequence,cter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        autocorr_full_aa('sam_input.csv',lg,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.to_csv(result_filename,index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        print('Value of lag:',lg,'\n')
        file1 = rest(Sequence,nrest,crest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        autocorr_full_aa('sam_input.csv',lg,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.to_csv(result_filename,index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0:
        print('Value of lag:',lg,'\n')
        file1 = restnc(Sequence,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        autocorr_full_aa('sam_input.csv',lg,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.to_csv(result_filename,index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        file1 = nct(Sequence,ncter)
        print('Value of lag:',lg,'\n')
        file1.to_csv('sam_input.csv', index=None, header=False)
        autocorr_full_aa('sam_input.csv',lg,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'NC'+df.columns
        df.to_csv(result_filename,index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        print('Value of lag:',lg,'\n')
        file1 =split(seq,sp)
        file1.to_csv('sam_input.csv', index=None, header=False)
        autocorr_split('sam_input.csv',sp,lg,result_filename)
        os.remove('sam_input.csv')
if Job == 'PAAC' or Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        print('Order of Gap:',lg,' ;Value of weight:',pw,'\n')
        paac(seq,lg,result_filename,pw)
        paac(seq,lg,'sam_allcomp.paac',pw)
    if nter != 0:
        print('Order of Gap:',lg,' ;Value of weight:',pw,'\n')
        file1 = nt(seq,nter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        paac('sam_input.csv',lg,'tempfile_out',pw)
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.paac_nt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        print('Order of Gap:',lg,' ;Value of weight:',pw,'\n')
        file1 = ct(seq,cter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        paac('sam_input.csv',lg,'tempfile_out',pw)
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.paac_ct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        print('Order of Gap:',lg,' ;Value of weight:',pw,'\n')
        file1 = rest(Sequence,nrest,crest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        paac('sam_input.csv',lg,'tempfile_out',pw)
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.paac_rt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0:
        print('Order of Gap:',lg,' ;Value of weight:',pw,'\n')
        file1 = restnc(Sequence,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        paac('sam_input.csv',lg,'tempfile_out',pw)
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.paac_rnc',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        print('Order of Gap:',lg,' ;Value of weight:',pw,'\n')
        file1 = nct(seq,ncter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        paac('sam_input.csv',lg,'tempfile_out',pw)
        df = pd.read_csv('tempfile_out')
        df.columns = 'NC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.paac_nct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        print('Order of Gap:',lg,' ;Value of weight:',pw,'\n')
        file1 =split(seq,sp)
        file1.to_csv('sam_input.csv', index=None, header=False)
        paac_split('sam_input.csv',sp,lg,result_filename,pw)
        paac_split('sam_input.csv',sp,lg,'sam_allcomp.paac_st',pw)
        os.remove('sam_input.csv')
if Job == 'APAAC' or Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        print('Order of Gap:',lg,' ;Value of weight:',pw,'\n')
        apaac(seq,lg,result_filename,pw)
        apaac(seq,lg,'sam_allcomp.apaac',pw)
    if nter != 0:
        print('Order of Gap:',lg,' ;Value of weight:',pw,'\n')
        file1 = nt(seq,nter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        apaac('sam_input.csv',lg,'tempfile_out',pw)
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.apaac_nt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        print('Order of Gap:',lg,' ;Value of weight:',pw,'\n')
        file1 = ct(seq,cter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        apaac('sam_input.csv',lg,'tempfile_out',pw)
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.apaac_ct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        print('Order of Gap:',lg,' ;Value of weight:',pw,'\n')
        file1 = rest(Sequence,nrest,crest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        apaac('sam_input.csv',lg,'tempfile_out',pw)
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.apaac_rt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0:
        print('Order of Gap:',lg,' ;Value of weight:',pw,'\n')
        file1 = restnc(Sequence,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        apaac('sam_input.csv',lg,'tempfile_out',pw)
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.apaac_rnc',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        print('Order of Gap:',lg,' ;Value of weight:',pw,'\n')
        file1 = nct(seq,ncter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        apaac('sam_input.csv',lg,'tempfile_out',pw)
        df = pd.read_csv('tempfile_out')
        df.columns = 'NC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.apaac_nct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        print('Order of Gap:',lg,' ;Value of weight:',pw,'\n')
        file1 =split(seq,sp)
        file1.to_csv('sam_input.csv', index=None, header=False)
        apaac_split('sam_input.csv',sp,lg,result_filename,pw)
        apaac_split('sam_input.csv',sp,lg,'sam_allcomp.apaac_st',pw)
        os.remove('sam_input.csv')
if Job == 'QSO' or Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        print('Order of Gap:',lg,' ;Value of weight:',wq,'\n')
        qos(seq,lg,result_filename,wq)
        qos(seq,lg,'sam_allcomp.qso',wq)
    if nter != 0:
        print('Order of Gap:',lg,' ;Value of weight:',wq,'\n')
        file1 = nt(seq,nter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        qos('sam_input.csv',lg,'tempfile_out',wq)
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.qso_nt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        print('Order of Gap:',lg,' ;Value of weight:',wq,'\n')
        file1 = ct(seq,cter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        qos('sam_input.csv',lg,'tempfile_out',wq)
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.qso_ct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        print('Order of Gap:',lg,' ;Value of weight:',wq,'\n')
        file1 = rest(Sequence,nrest,crest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        qos('sam_input.csv',lg,'tempfile_out',wq)
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.qso_rt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0:
        print('Order of Gap:',lg,' ;Value of weight:',wq,'\n')
        file1 = restnc(Sequence,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        qos('sam_input.csv',lg,'tempfile_out',wq)
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.qso_rnc',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        print('Order of Gap:',lg,' ;Value of weight:',wq,'\n')
        file1 = nct(seq,ncter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        qos('sam_input.csv',lg,'tempfile_out',wq)
        df = pd.read_csv('tempfile_out')
        df.columns = 'NC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.qso_nct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        print('Order of Gap:',lg,' ;Value of weight:',wq,'\n')
        file1 =split(seq,sp)
        file1.to_csv('sam_input.csv', index=None, header=False)
        qos_split('sam_input.csv',sp,lg,result_filename,wq)
        qos_split('sam_input.csv',sp,lg,'sam_allcomp.qso_st',wq)
        os.remove('sam_input.csv')
if Job == 'SOC' or Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        print('Order of Gap:',lg,'\n')
        soc(seq,lg,result_filename)
        soc(seq,lg,'sam_allcomp.soc')
    if nter != 0:
        print('Order of Gap:',lg,'\n')
        file1 = nt(seq,nter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        soc('sam_input.csv',lg,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.soc_nt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        print('Order of Gap:',lg,'\n')
        file1 = ct(seq,cter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        soc('sam_input.csv',lg,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.soc_ct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        print('Order of Gap:',lg,'\n')
        file1 = rest(Sequence,nrest,crest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        soc('sam_input.csv',lg,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.soc_rt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0:
        print('Order of Gap:',lg,'\n')
        file1 = restnc(Sequence,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        soc('sam_input.csv',lg,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.soc_rnc',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        print('Order of Gap:',lg,'\n')
        file1 = nct(seq,ncter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        soc('sam_input.csv',lg,'tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'NC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.soc_nct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        print('Order of Gap:',lg,'\n')
        file1 =split(seq,sp)
        file1.to_csv('sam_input.csv', index=None, header=False)
        soc_split('sam_input.csv',sp,lg,result_filename)
        soc_split('sam_input.csv',sp,lg,'sam_allcomp.soc_st')
        os.remove('sam_input.csv')
if Job == 'CTC' or Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        CTC(seq,result_filename)
        CTC(seq,'sam_allcomp.ctc')
    if nter != 0:
        print('Order of Gap:',lg,'\n')
        file1 = nt(seq,nter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        CTC('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.ctc_nt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        print('Order of Gap:',lg,'\n')
        file1 = ct(seq,cter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        CTC('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.ctc_ct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        print('Order of Gap:',lg,'\n')
        file1 = rest(Sequence,nrest,crest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        CTC('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.ctc_rt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0:
        print('Order of Gap:',lg,'\n')
        file1 = restnc(Sequence,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        CTC('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.ctc_rnc',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        print('Order of Gap:',lg,'\n')
        file1 = nct(seq,ncter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        CTC('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'NC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.ctc_nct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        print('Order of Gap:',lg,'\n')
        file1 =split(seq,sp)
        file1.to_csv('sam_input.csv', index=None, header=False)
        ctc_split('sam_input.csv',sp,result_filename)
        ctc_split('sam_input.csv',sp,'sam_allcomp.ctc_st')
        os.remove('sam_input.csv')
if Job == 'CETD' or Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        ctd(seq,result_filename)
        ctd(seq,'sam_allcomp.ctd')
    if nter != 0:
        print('Order of Gap:',lg,'\n')
        file1 = nt(seq,nter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        ctd('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'N'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.ctd_nt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if cter != 0:
        print('Order of Gap:',lg,'\n')
        file1 = ct(seq,cter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        ctd('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'C'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.ctd_ct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if nrest != 0 or crest != 0:
        print('Order of Gap:',lg,'\n')
        file1 = rest(Sequence,nrest,crest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        ctd('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'R'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.ctd_rt',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncrest != 0:
        print('Order of Gap:',lg,'\n')
        file1 = restnc(Sequence,ncrest)
        file1.to_csv('sam_input.csv', index=None, header=False)
        ctd('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'RNC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.ctd_rnc',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if ncter != 0:
        print('Order of Gap:',lg,'\n')
        file1 = nct(seq,ncter)
        file1.to_csv('sam_input.csv', index=None, header=False)
        ctd('sam_input.csv','tempfile_out')
        df = pd.read_csv('tempfile_out')
        df.columns = 'NC'+df.columns
        df.to_csv(result_filename,index=None)
        df.to_csv('sam_allcomp.ctd_nct',index=None)
        os.remove('sam_input.csv')
        os.remove('tempfile_out')
    if sp != 0:
        print('Order of Gap:',lg,'\n')
        file1 =split(seq,sp)
        file1.to_csv('sam_input.csv', index=None, header=False)
        ctd_split('sam_input.csv',sp,result_filename)
        ctd_split('sam_input.csv',sp,'sam_allcomp.ctd_st')
        os.remove('sam_input.csv')
if Job == 'ALLCOMP':
    if (nter == 0) and (cter == 0) and (nrest == 0) and (crest == 0) and (sp == 0) and (ncter == 0) and (ncrest == 0):
        df1 = pd.read_csv("sam_allcomp.aac")
        df2 = pd.read_csv("sam_allcomp.dpc")
        df3 = pd.read_csv("sam_allcomp.tpc")
        df4 = pd.read_csv("sam_allcomp.atc")
        df5 = pd.read_csv("sam_allcomp.btc")
        df6 = pd.read_csv("sam_allcomp.pcp")
        df7 = pd.read_csv("sam_allcomp.rri")
        df8 = pd.read_csv("sam_allcomp.pri")
        df9 = pd.read_csv("sam_allcomp.ddr")
        df10 = pd.read_csv("sam_allcomp.sep")
        df11 = pd.read_csv("sam_allcomp.ser")
        df12 = pd.read_csv("sam_allcomp.spc")
        df13 = pd.read_csv("sam_allcomp.ctc")
        df14 = pd.read_csv("sam_allcomp.ctd")
        df15 = pd.read_csv("sam_allcomp.paac")
        df16 = pd.read_csv("sam_allcomp.apaac")
        df17 = pd.read_csv("sam_allcomp.qso")
        df18 = pd.read_csv("sam_allcomp.soc")
        df19 = pd.concat([df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17,df18],axis=1)
        df19.to_csv(result_filename, index=None)
    if nter != 0:
        df1 = pd.read_csv("sam_allcomp.aac_nt")
        df2 = pd.read_csv("sam_allcomp.dpc_nt")
        df3 = pd.read_csv("sam_allcomp.tpc_nt")
        df4 = pd.read_csv("sam_allcomp.atc_nt")
        df5 = pd.read_csv("sam_allcomp.btc_nt")
        df6 = pd.read_csv("sam_allcomp.pcp_nt")
        df7 = pd.read_csv("sam_allcomp.rri_nt")
        df8 = pd.read_csv("sam_allcomp.pri_nt")
        df9 = pd.read_csv("sam_allcomp.ddr_nt")
        df10 = pd.read_csv("sam_allcomp.sep_nt")
        df11 = pd.read_csv("sam_allcomp.ser_nt")
        df12 = pd.read_csv("sam_allcomp.spc_nt")
        df13 = pd.read_csv("sam_allcomp.ctc_nt")
        df14 = pd.read_csv("sam_allcomp.ctd_nt")
        df15 = pd.read_csv("sam_allcomp.paac_nt")
        df16 = pd.read_csv("sam_allcomp.apaac_nt")
        df17 = pd.read_csv("sam_allcomp.qso_nt")
        df18 = pd.read_csv("sam_allcomp.soc_nt")
        df19 = pd.concat([df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17,df18],axis=1)
        df19.to_csv(result_filename, index=None)
    if cter != 0:
        df1 = pd.read_csv("sam_allcomp.aac_ct")
        df2 = pd.read_csv("sam_allcomp.dpc_ct")
        df3 = pd.read_csv("sam_allcomp.tpc_ct")
        df4 = pd.read_csv("sam_allcomp.atc_ct")
        df5 = pd.read_csv("sam_allcomp.btc_ct")
        df6 = pd.read_csv("sam_allcomp.pcp_ct")
        df7 = pd.read_csv("sam_allcomp.rri_ct")
        df8 = pd.read_csv("sam_allcomp.pri_ct")
        df9 = pd.read_csv("sam_allcomp.ddr_ct")
        df10 = pd.read_csv("sam_allcomp.sep_ct")
        df11 = pd.read_csv("sam_allcomp.ser_ct")
        df12 = pd.read_csv("sam_allcomp.spc_ct")
        df13 = pd.read_csv("sam_allcomp.ctc_ct")
        df14 = pd.read_csv("sam_allcomp.ctd_ct")
        df15 = pd.read_csv("sam_allcomp.paac_ct")
        df16 = pd.read_csv("sam_allcomp.apaac_ct")
        df17 = pd.read_csv("sam_allcomp.qso_ct")
        df18 = pd.read_csv("sam_allcomp.soc_ct")
        df19 = pd.concat([df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17,df18],axis=1)
        df19.to_csv(result_filename, index=None)
    if nrest != 0 or crest != 0:
        df1 = pd.read_csv("sam_allcomp.aac_rt")
        df2 = pd.read_csv("sam_allcomp.dpc_rt")
        df3 = pd.read_csv("sam_allcomp.tpc_rt")
        df4 = pd.read_csv("sam_allcomp.atc_rt")
        df5 = pd.read_csv("sam_allcomp.btc_rt")
        df6 = pd.read_csv("sam_allcomp.pcp_rt")
        df7 = pd.read_csv("sam_allcomp.rri_rt")
        df8 = pd.read_csv("sam_allcomp.pri_rt")
        df9 = pd.read_csv("sam_allcomp.ddr_rt")
        df10 = pd.read_csv("sam_allcomp.sep_rt")
        df11 = pd.read_csv("sam_allcomp.ser_rt")
        df12 = pd.read_csv("sam_allcomp.spc_rt")
        df13 = pd.read_csv("sam_allcomp.ctc_rt")
        df14 = pd.read_csv("sam_allcomp.ctd_rt")
        df15 = pd.read_csv("sam_allcomp.paac_rt")
        df16 = pd.read_csv("sam_allcomp.apaac_rt")
        df17 = pd.read_csv("sam_allcomp.qso_rt")
        df18 = pd.read_csv("sam_allcomp.soc_rt")
        df19 = pd.concat([df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17,df18],axis=1)
        df19.to_csv(result_filename, index=None)
    if ncrest != 0:
        df1 = pd.read_csv("sam_allcomp.aac_rnc")
        df2 = pd.read_csv("sam_allcomp.dpc_rnc")
        df3 = pd.read_csv("sam_allcomp.tpc_rnc")
        df4 = pd.read_csv("sam_allcomp.atc_rnc")
        df5 = pd.read_csv("sam_allcomp.btc_rnc")
        df6 = pd.read_csv("sam_allcomp.pcp_rnc")
        df7 = pd.read_csv("sam_allcomp.rri_rnc")
        df8 = pd.read_csv("sam_allcomp.pri_rnc")
        df9 = pd.read_csv("sam_allcomp.ddr_rnc")
        df10 = pd.read_csv("sam_allcomp.sep_rnc")
        df11 = pd.read_csv("sam_allcomp.ser_rnc")
        df12 = pd.read_csv("sam_allcomp.spc_rnc")
        df13 = pd.read_csv("sam_allcomp.ctc_rnc")
        df14 = pd.read_csv("sam_allcomp.ctd_rnc")
        df15 = pd.read_csv("sam_allcomp.paac_rnc")
        df16 = pd.read_csv("sam_allcomp.apaac_rnc")
        df17 = pd.read_csv("sam_allcomp.qso_rnc")
        df18 = pd.read_csv("sam_allcomp.soc_rnc")
        df19 = pd.concat([df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17,df18],axis=1)
        df19.to_csv(result_filename, index=None)
    if ncter != 0:
        df1 = pd.read_csv("sam_allcomp.aac_nct")
        df2 = pd.read_csv("sam_allcomp.dpc_nct")
        df3 = pd.read_csv("sam_allcomp.tpc_nct")
        df4 = pd.read_csv("sam_allcomp.atc_nct")
        df5 = pd.read_csv("sam_allcomp.btc_nct")
        df6 = pd.read_csv("sam_allcomp.pcp_nct")
        df7 = pd.read_csv("sam_allcomp.rri_nct")
        df8 = pd.read_csv("sam_allcomp.pri_nct")
        df9 = pd.read_csv("sam_allcomp.ddr_nct")
        df10 = pd.read_csv("sam_allcomp.sep_nct")
        df11 = pd.read_csv("sam_allcomp.ser_nct")
        df12 = pd.read_csv("sam_allcomp.spc_nct")
        df13 = pd.read_csv("sam_allcomp.ctc_nct")
        df14 = pd.read_csv("sam_allcomp.ctd_nct")
        df15 = pd.read_csv("sam_allcomp.paac_nct")
        df16 = pd.read_csv("sam_allcomp.apaac_nct")
        df17 = pd.read_csv("sam_allcomp.qso_nct")
        df18 = pd.read_csv("sam_allcomp.soc_nct")
        df19 = pd.concat([df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17,df18],axis=1)
        df19.to_csv(result_filename, index=None)
    if sp != 0:
        df1 = pd.read_csv("sam_allcomp.aac_st")
        df2 = pd.read_csv("sam_allcomp.dpc_st")
        df3 = pd.read_csv("sam_allcomp.tpc_st")
        df4 = pd.read_csv("sam_allcomp.atc_st")
        df5 = pd.read_csv("sam_allcomp.btc_st")
        df6 = pd.read_csv("sam_allcomp.pcp_st")
        df7 = pd.read_csv("sam_allcomp.rri_st")
        df8 = pd.read_csv("sam_allcomp.pri_st")
        df9 = pd.read_csv("sam_allcomp.ddr_st")
        df10 = pd.read_csv("sam_allcomp.sep_st")
        df11 = pd.read_csv("sam_allcomp.ser_st")
        df12 = pd.read_csv("sam_allcomp.spc_st")
        df13 = pd.read_csv("sam_allcomp.ctc_st")
        df14 = pd.read_csv("sam_allcomp.ctd_st")
        df15 = pd.read_csv("sam_allcomp.paac_st")
        df16 = pd.read_csv("sam_allcomp.apaac_st")
        df17 = pd.read_csv("sam_allcomp.qso_st")
        df18 = pd.read_csv("sam_allcomp.soc_st")
        df19 = pd.concat([df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17,df18],axis=1)
        df19.to_csv(result_filename, index=None)
filelist=glob.glob("sam_allcomp*")
for file_2 in filelist:
    os.remove(file_2)
os.remove(file_output)
