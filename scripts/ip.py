#!/usr/bin/python
import sys 
import os
import numpy as np
import math
from math import log
dicto2=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X']
dicto={'A':6.01,'C':5.07,'D':2.77,'E':3.22,'F':5.48,'G':5.97,'H':7.59,'I':6.02,'K':9.74,'L':5.98,'M':5.47,'N':5.41,'P':6.48,'Q':5.65,'R':10.76,'S':5.68,'T':5.87,'V':5.97,'W':5.89,'Y':5.67,'X':0 }


text_file2= open("ip.txt", "w")
name2=str(sys.argv[1])
fp=open(name2,'r')
def ip(seq):
     win=7
     win2=int((win-1)/2)
     ip_list=[]
     line0=''.join(['X'*win2,line,'X'*win2])
     f=  [line0[i:i+win] for i in range(len(line0)-win+1)]  
     for t in f:    
               o=list(set(t)) 
               pI=0     
               for s in o :
                 if s in dicto:
                      number = t.count(s) 
                      pI += number*dicto[s]
               pI2= float(pI)/len(t)
               ip_list.append(pI2)
                    
     return ip_list
for line in fp:
           line=line.replace('\n','')
           line=line.replace(' ','')
           final=ip(line)
           for i in final:
               text_file2.write("%s \n" %(i))
text_file2.close()               
