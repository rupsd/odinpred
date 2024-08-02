#!/usr/bin/python

#from __future__ import division
import sys 
import math
from math import log


win=51
win2=int((win-1)/2)
text_file2= open("entropy.txt", "w")
name2=str(sys.argv[1])
fp=open(name2,'r')

  
for line in fp:
      line=line.replace('\n','')
      line=line.replace(' ','')
      line0=''.join(['X'*win2,line,'X'*win2])
      f=  [line0[i:i+win] for i in range(len(line0)-win+1)]  
      for t in f:  
               entropy=0  
               o=list(set(t))          
               for s in o :
                      number = t.count(s) 
                      fraction1 = float(number)/(len(t)-1)
                      entropy+=(fraction1*log(fraction1,2))      
               entropy=-entropy
               text_file2.write("%s \n" %(entropy))
text_file2.close()              
