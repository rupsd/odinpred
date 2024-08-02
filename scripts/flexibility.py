#!/usr/bin/python

import sys
import math
from math import log


dicto={'A':0.984,'C':0.906,'D':1.068,'E':1.094,'F':0.915,'G':1.031,'H':0.950,'I':0.927,'K':1.102,'L':0.935,'M':0.952,'N':1.048,'P':1.049,'Q':1.037,'R':1.008,'S':1.046,'T':0.997,'V':0.931,'W':0.904,'Y':0.929,'x':0.9,'X':0.9}

win=7
win2=int((win-1)/2)
ct=0
text_file2= open("flexibility.txt", "w")
text_file2= open("flexibility.txt", "r+")
name2=str(sys.argv[1])
fp=open(name2,'r')

for line in fp:
      line=line.replace('\n','')
      line=line.replace(' ','')
      line0=''.join(['X'*win2,line,'X'*win2])
      f=  [line0[i:i+win] for i in range(len(line0)-win+1)]  
      for t in f: 
               F1=0 
               for o in range(3,len(t)-3):
                   
                  
                   
                   F1+=((dicto[t[o]]+0.75*(dicto[t[o-1]]+dicto[t[o+1]])+0.5*(dicto[t[o-2]]+dicto[t[o+2]])+0.25*(dicto[t[o-3]]+dicto[t[o+3]]))/4)
                   text_file2.write("%s \n" %(F1))
                   
text_file2.close()               
   
