#!/usr/bin/python
from __future__ import division
import sys 
import math
from math import log
import CFdef0 as CF
import subprocess
from subprocess import call



text_file= open("Chou.txt", "w")
ct=0
##dict0=['A','B','T']   ####### turn,beta,helix
dict0=['T','B','A']
win=15
win2=int((win-1)/2)
seq=str(sys.argv[1])
for s in dict0:
 ct=ct+1
 fp1=open(seq,'r')
 for line in fp1:
     line=line.replace('\n','')
     line=line.replace(' ','')
     text_file2= open("CF.txt", "w")
     y=CF.ChouFasman(line)    #############yttttttbbbbaaaa
     y=str(y)
     y=y.replace(' ','X')
     line0=''.join(['X'*win2,y,'X'*win2])
     y1=  [line0[i:i+win] for i in range(len(line0)-win+1)]
     for t in y1:
            if s in t:
                    number = t.count(s)/len(t) 
                    text_file2.write("%s \n" %(number))
            else:
                     number=0         
                     text_file2.write("%s \n" %(number))
 text_file2.close()
 subprocess.call("paste Chou.txt CF.txt > test.txt && mv test.txt Chou.txt", shell=True)   
text_file.close()
