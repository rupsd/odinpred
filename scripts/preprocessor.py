#!/usr/bin/python
from __future__ import division
import sys 
import re


line=str(sys.argv[1])
fp2=open('input_seq2.txt','w')
fp3=open('pred_output.txt','w')
all_residues=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X']

line=line.replace(' ','')   ###Remove space
line=line.replace('\n','')  
line=''.join([i for i in line if i.isalpha()])   ## Remove everything other than alphabets

line3=line.upper()  ## change things to upper
line3=[x if x in all_residues else 'A' for x in line3]  ## change all alphabets other than amino acids to A
line3=''.join(line3) # Join everything to a string



fp2.write("%s\n" %(str(line3)))
for k in line3:
    fp3.write("%s\n" %(str(k)))
fp3.close()
fp2.close()
