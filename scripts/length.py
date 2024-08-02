#!/usr/bin/python
import sys 

import subprocess
import math
from math import log

text_file2= open("length.txt", "w")
name2=str(sys.argv[1])
fp=open(name2,'r')

  
for line in fp:
      u=len(line)-1
      G=math.log10(min(500,u))
      for i in range(u):
              ##G=math.log10(u)
              text_file2.write("%s \n" %(G))

text_file2.close()              
