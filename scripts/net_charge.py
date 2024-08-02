#!/usr/bin/python
import sys 
import os
import subprocess
from subprocess import call
import numpy as np
import math
from math import log
n={'K','R'}
p={"D", "E"}
win=9
win2=int((win-1)/2)
text_file2= open("net_charge.txt", "w")
name2=str(sys.argv[1])
fp=open(name2,'r')
dicto=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X']
  
for line in fp:
      line=line.replace('\n','')
      line=line.replace(' ','')
      line0=''.join(['X'*win2,line,'X'*win2])
      f=  [line0[i:i+win] for i in range(len(line0)-win+1)]  
      for t in f:   
               number=0  
               fraction1=0 
               for s in t :
                if s in dicto:
                    if s in p: 
                      number +=1 
                      fraction1 = float(number)/(len(t))
                    if s in n:
                      number -=1 
                      fraction1 = float(number)/(len(t))
               text_file2.write("%s \n" %(fraction1))
      text_file2.close()              
