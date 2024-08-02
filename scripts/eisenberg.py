#!/usr/bin/python
from __future__ import division
import sys 
import os
import math
text_file= open("eis_alpha.txt", "w")
text_file2= open("eis_beta.txt", "w")

win=9
win2=int((win-1)/2)
hydroatph7 = {'F': 0.8111, 'I': 1.0, 'W': 0.4, 'L': 0.9222, 'V': 0.9667, 'M': 0.7111, 'Y': 0.3556, 'C': 0.7778, 'A': 0.7000, 'T':0.4222, 'H': 0.1444, 'G':0.4556, 'S':0.4111, 'Q':0.1111, 'R':0, 'L':0.9222, 'N': 0.1111,'E':0.1111, 'P':0.3222,'D':0.1111, 'K': 0.0667}
name2=str(sys.argv[1])
fp=open(name2,'r')


for line in fp:
    line=line.replace('\n','')
    line=line.replace(' ','')
    line0=''.join(['X'*win2,line,'X'*win2])
    f=  [line0[i:i+win] for i in range(len(line0)-win+1)]
    for t in f:
            hydrophobicity=0
            hydro_moment_ww=0
            hydro_moment_ww_b=0
            si=0
            co=0
            sib=0
            cob=0
            i=0
            for s in t:                               
              if s in hydroatph7:
                      hydrophobicity += hydroatph7[s]
                      i=i+1
                      co +=hydroatph7[s]*(math.cos(math.radians(i*100)))
                      si +=hydroatph7[s]*(math.sin(math.radians(i*100)))
                      cob +=hydroatph7[s]*(math.cos(math.radians(i*170)))
                      sib +=hydroatph7[s]*(math.sin(math.radians(i*170)))
                      
            hydro_moment_ww += abs(math.sqrt(co**2 + si**2))
            hydro_moment_ww2=hydro_moment_ww /len(t)
            text_file.write("%s\n" %(hydro_moment_ww2))
            hydro_moment_ww_b += abs(math.sqrt(cob**2 + sib**2))
            hydro_moment_ww2_b=hydro_moment_ww_b /len(t)
            text_file2.write("%s\n" %(hydro_moment_ww2_b))


            
text_file.close()
text_file2.close()
        
