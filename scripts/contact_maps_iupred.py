#!/usr/bin/python
from __future__ import division
import sys
import numpy as np
import math
from math import log
import os
import pandas as pd

text_file2 = open("energy_25_again.txt", "w")
win = 25
win2 = int((win - 1) / 2)
order = {'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17, 'V': 18, 'W': 19, 'Y': 20, 'X': 21}
script_dir = os.path.dirname(os.path.abspath(__file__))

name2 = str(sys.argv[1])
fp = open(name2, 'r')
print(name2)
x_file = open(os.path.join(script_dir, "matrixP.txt"), 'r')
lines = x_file.readlines()

for line in fp:
    W = []
    line = line.replace('\n', '')
    line = line.replace(' ', '')
    line0 = ''.join(['X' * win2, line, 'X' * win2])
    f = [line0[i:i + win] for i in range(len(line0) - win + 1)]
    
    for t in f:
        if 'X' in t:
            t = t.replace('X', '')
        mat = []
        for y in range(len(t)):
            lst = []
            for x in range(len(t)):
                i = order[t[x]]
                j = order[t[y]]
                f1 = lines[j - 1]
                aa2 = f1.split()
                val = [eval(s) for s in aa2]
                l = val[i - 1]
                m = math.exp(-l / 4.26)
                lst.append(m)
            mat.append(lst)
        M = np.array(mat)
        H = np.sum(np.log(M) * M)
        text_file2.write("%6.4f\n" % (H))

text_file2.close()
