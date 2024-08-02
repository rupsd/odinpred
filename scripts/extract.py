#!/usr/bin/python

import numpy as np
from numpy import *
from numpy import isnan
import math
from math import log
import sys 
import os


i=str(sys.argv[1])
mnop=sys.argv[2]
win=11
for m in range(0,9):
 text_file= open('evolv_featuresunity_training'+str(m)+'.txt', "w")
 with open(mnop) as fp:   ### i is
             r=[]
             r2=[]
             y=[]
             for line in fp:
                   
                  r.append(line.split()[m])
             
             r = [0 if x=='nan' else x for x in r]
             r2.extend(r[0:((win-1)/2)])
             r2.extend(r)
             r2.extend(r[len(r)-((win-1)/2): len(r)])
             r3=list(r2)
             f=[r3[k:k+win] for k in xrange(len(r3)-win+1)]
           
             for t in f:
                    s=np.array(t).astype(np.float)
                    average= np.mean(s)
                    text_file.write("%s \n" %(average))
 text_file.close()
 
win=17
for m in range(9,29):
 text_file= open('evolv_featuresunity_training'+str(m)+'.txt', "w")
 with open(mnop) as fp:
             r=[]
             r2=[]
             y=[]
             for line in fp:
                  r.append(line.split()[m])
             r = [0 if x=='nan' else x for x in r]
             r2.extend(r[0:((win-1)/2)])
             r2.extend(r)
             r2.extend(r[len(r)-((win-1)/2): len(r)])
             r3=list(r2)
             f=[r3[k:k+win] for k in xrange(len(r3)-win+1)]
             #print len(f)
             for t in f:
                    s=np.array(t).astype(np.float)
                    average= np.mean(s)
                    text_file.write("%s \n" %(average))


 text_file.close()                   
                  
            
   
