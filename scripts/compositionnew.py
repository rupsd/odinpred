#!/usr/bin/python
# https://cran.r-project.org/web/packages/protr/vignettes/protr.html

from __future__ import division
import sys
import re
import math
import os
import subprocess
from numpy import array, log, cumsum
from subprocess import call

name2 = str(sys.argv[1])
fp = open(name2, 'r')

Hydro_polar = {'R', 'K', 'E', 'D', 'Q', 'N'}
Hydro_neutral = {'G', 'A', 'S', 'T', 'P', 'H', 'Y'}
hydro = {'C', 'L', 'V', 'I', 'M', 'F', 'W'}

charge1 = {'K', 'R'}
charge2 = {'A', 'N', 'C', 'Q', 'G', 'H', 'I', 'L', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}
charge3 = {'D', 'E'}

dipole1 = {'A', 'G', 'V'}
dipole2 = {'I', 'L', 'F', 'P'}
dipole3 = {'T', 'Y', 'M', 'S'}
dipole4 = {'H', 'N', 'Q', 'W'}
dipole5 = {'R', 'L'}
dipole6 = {'D', 'E'}
dipole7 = {'C'}

polarity1 = {'L', 'I', 'F', 'W', 'C', 'M', 'V', 'Y'}
polarity2 = {'P', 'A', 'T', 'G', 'S'}
polarity3 = {'H', 'Q', 'R', 'K', 'N', 'E', 'D'}

polarizibility1 = {'G', 'A', 'S', 'D', 'T'}
polarizibility2 = {'C', 'P', 'N', 'V', 'E', 'Q', 'I', 'L'}
polarizibility3 = {'M', 'H', 'K', 'F', 'R', 'Y', 'W'}

helix = {'E', 'A', 'L', 'M', 'Q', 'K', 'R', 'H'}
strand = {'V', 'I', 'Y', 'C', 'W', 'F', 'T'}
coil = {'G', 'N', 'P', 'S', 'D'}

solvent_buried = {'A', 'L', 'F', 'C', 'G', 'I', 'V', 'W'}
solvent_exposed = {'R', 'K', 'Q', 'E', 'N', 'D'}
solvent_intermediate = {'M', 'S', 'P', 'T', 'H', 'Y'}

# alphasyn

def seq2state(seq, groups):
    st = []
    for s in seq:
        isin = False
        for i, gri in enumerate(groups):
            if s in gri:
                st.append(i)
                isin = True
                # break
        if not isin:
            print('warning: not in any group', s)
    return st

def st2trans(st):
    conv = {(0, 0): 0, (1, 1): 0, (2, 2): 0, (0, 1): 1, (1, 0): 1, (0, 2): 2, (2, 0): 2, (2, 1): 3, (1, 2): 3}
    tr = []
    for i in range(len(st) - 1):
        pair = st[i], st[i + 1]
        tr.append(conv[pair])
    tr.append(0)  # to make length identical to seq len
    return tr

def calc_fraction(st, win, num):
    allcounts = []
    for i in range(num):
        ari = array(st, dtype=int) == i
        si = cumsum(ari)
        counts = [si[k + win] - si[k] for k in range(len(st) - win)]
        allcounts.append(array([counts[0]] * int(win / 2 + 1) + counts + [counts[-1]] * int(win - win / 2)) * 1.0 / win)
    return allcounts

def computeCompTrans(seq, group, win=25):
    st = seq2state(seq, group)
    print(''.join([str(x) for x in st]))
    tr = st2trans(st)
    print(''.join([str(x) for x in tr]))
    comp = calc_fraction(st, win, 3)
    trans4 = calc_fraction(tr, win, 4)
    T = array(trans4)
    entr = []
    for i in range(len(st)):
        fri = T[:, i]
        prod = fri * log(fri)
        use = prod < 999
        dont = array(1 - use, dtype=bool)
        prod[dont] = 0
        entr.append(sum(prod))
    return comp, trans4[1:], entr

vw1 = {'G', 'A', 'S', 'T', 'P', 'D', 'C'}
vw2 = {'N', 'V', 'E', 'Q', 'I', 'L'}
vw3 = {'M', 'H', 'K', 'F', 'R', 'Y', 'W'}

outfile = open('comptransnew.txt', 'w')

for line1 in fp:
    line1 = line1.replace('\n', '')
    line1 = line1.replace(' ', '')
    line1 = line1.replace('X', 'A')
    GR1 = (Hydro_polar, Hydro_neutral, hydro)
    GR2 = (charge1, charge2, charge3)
    GR3 = (polarity1, polarity2, polarity3)
    GR4 = (polarizibility1, polarizibility2, polarizibility3)
    GR5 = (helix, strand, coil)
    GR6 = (solvent_buried, solvent_exposed, solvent_intermediate)
    GR7 = (vw1, vw2, vw3)
    allgroups = (GR1, GR2, GR3, GR4, GR5, GR6, GR7)
    features = []
    seq = line1
    for gri in allgroups:
        comp, trans, entr = computeCompTrans(seq, gri, win=25)
        features.append((comp, trans, entr))
    for i in range(len(seq)):
        for ig in range(len(allgroups)):
            for n in range(3):
                fign = features[ig][n]
                if n < 2:
                    values = [fx[i] for fx in fign]
                    outfile.write('%6.4f %6.4f %6.4f ' % tuple(values))
                else:
                    outfile.write('%7.4f ' % fign[0])
        outfile.write('\n')

outfile.close()
