#!/usr/bin/python
#########################################################
# This script calculates number of repeats in a protein (doublets/triplets etc.)
from __future__ import division
import sys
import subprocess
import string
from numpy import *
import numpy as np
from math import factorial

win = 9
win2 = int((win - 1) / 2)
extension = [0] * win2

text_file1 = open("repeatslist20.txt", "w")
text_file2 = open("repeatslist201.txt", "w")
text_file3 = open("repeatslist202.txt", "w")
text_file4 = open("repeatslist203.txt", "w")
text_file5 = open("repeatslist204.txt", "w")

aa13dict = {'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY', 'F': 'PHE', 'I': 'ILE', 'H': 'HIS', 'K': 'LYS', 
            'M': 'MET', 'L': 'LEU', 'N': 'ASN', 'Q': 'GLN', 'P': 'PRO', 'S': 'SER', 'R': 'ARG', 'T': 'THR', 'W': 'TRP', 
            'V': 'VAL', 'Y': 'TYR', 'X': 'NaN'}

aa31dict = {'CYS': 'C', 'GLN': 'Q', 'ILE': 'I', 'SER': 'S', 'VAL': 'V', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'LYS': 'K', 
            'THR': 'T', 'PHE': 'F', 'ALA': 'A', 'HIS': 'H', 'GLY': 'G', 'ASP': 'D', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
            'GLU': 'E', 'TYR': 'Y', 'NaN': 'X'}

# Convert dictionary keys to list and sort
aa1s = list(aa13dict.keys())
aa1s.sort()

# from 9-gram
# HighlyConserved=['C'] #0
# Proline=['P'] #1
# PolarAcids=['N','Q'] #2
# polar=['D','E'] #3
# Alcohols=['S','T'] #4
# Gly=['G'] #5
# Aliphatic=['A','I','V','L','M'] #6
# Aromatic=['F','Y','W'] #7
# Bases=['H','K','R'] #8

dct9 = {'C': '0', 'M': '6', 'N': '2', 'Q': '2', 'D': '3', 'E': '3', 'S': '4', 'T': '4', 'P': '1', 'A': '6', 'G': '5', 
        'I': '6', 'V': '6', 'L': '6', 'F': '7', 'Y': '7', 'W': '7', 'H': '8', 'K': '8', 'R': '8', '-': '-', 'X': '-'}

def expprob(N, l, p, k):
    h = (N - l + 1) * p
    return exp(-h) * pow(h, k) / factorial(k)

name2 = str(sys.argv[1])
fp = open(name2, 'r')

for line in fp:
    seq = ''
    probcut = 0.1
    line = line.replace('\n', '')
    line = line.replace(' ', '')
    if len(line) > 0 and line[0] != ">":
        seq += line
    aacount = [seq.count(aa) for aa in aa1s]
    fracs = array(aacount) * 1.0 / len(seq)
    # print fracs

    N = len(seq)
    nums = zeros(len(seq), dtype=int)
    mult = zeros(len(seq), dtype=int)
    probs = ones(len(seq))
    dists = zeros(len(seq))
    symms = zeros(len(seq))
    groups = [None for _ in range(N)]
    
    def getrep(seq, delta, pos):
        # print '--------------------',delta,pos
        if len(pos) == 0:
            return None
        pairs = [seq[i:i+delta] for i in pos]
        sp = set(pairs)
        lp = array(list(sp))
        # print len(sp)
        # print sp
        cp = [pairs.count(spi) for spi in sp]
        rep = array(cp) > 1
        reps = lp[rep]
        # print len(reps)
        # print reps
        ap = array(pairs)
        repeats = [[r2i, pos[ap == r2i], []] for r2i in reps]
        for tup in repeats:
            segm, posnew, children = tup
            L = len(segm)
            M = len(posnew)
            child = getrep(seq, delta + 1, posnew)
            if len(child) == 0:
                for pn in posnew:
                    for i in range(pn, pn + delta):
                        comb = prod([fracs[aa1s.index(sm)] for sm in segm])
                        prob = expprob(N, L, comb, M)
                        if prob < probs[i]:
                            nums[i] = L
                            mult[i] = M
                            # prob = expprob(N, L, 0.05**L, M)
                            probs[i] = prob
                            groups[i] = segm
                            devs = [posnew[n + 1] - posnew[n] for n in range(len(posnew) - 1)]
                            dists[i] = 1.0 / (max(1, average(devs) - L))
                            if len(posnew) > 2:
                                symms[i] = 1.0 / (max(1, std(devs)))
                children.append(child)
        return repeats

    getrep(seq, 2, arange(len(seq) - 1))
    alpha = '-.23456789abcdefghijklmnopqrstuvxyz'
    deg = ''.join([alpha[n] for n in nums])
    # print string.join([alpha[n] for n in mult], '')

    removed = probs > probcut
    dists[removed] = 0.0
    symms[removed] = 0.0

    nums = np.append(extension, nums)
    nums = np.append(nums, extension)
    ct = 0
    fnums = [nums[i:i + win] for i in range(len(nums) - win + 1)]
    for t in fnums:
        x = sum(t) / len(t)
        # kprobs=np.append(kprobs,t/len(t))
        text_file4.write("%s \n" % (x))

    mult = np.append(extension, mult)
    mult = np.append(mult, extension)
    ct = 0
    fmult = [mult[i:i + win] for i in range(len(mult) - win + 1)]
    for t in fmult:
        x = sum(t) / len(t)
        # kprobs=np.append(kprobs,t/len(t))
        text_file5.write("%s \n" % (x))

    probs = np.append(extension, probs)
    probs = np.append(probs, extension)
    ct = 0
    fprobs = [probs[i:i + win] for i in range(len(probs) - win + 1)]
    for t in fprobs:
        x = sum(t) / len(t)
        # kprobs=np.append(kprobs,t/len(t))
        text_file1.write("%s \n" % (-log(x)))
    # 3
    print(ct)
    dists = np.append(extension, dists)
    dists = np.append(dists, extension)
    fdists = [dists[i:i + win] for i in range(len(dists) - win + 1)]
    kdists = []
    for t in fdists:
        y = sum(t) / len(t)
        text_file2.write("%s \n" % (y))

    symms = np.append(extension, symms)
    symms = np.append(symms, extension)
    fsymms = [symms[i:i + win] for i in range(len(symms) - win + 1)]
    ksymms = []
    for t in fsymms:
        z = sum(t) / len(t)
        text_file3.write("%s \n" % (z))

text_file1.close()
text_file2.close()
text_file3.close()
text_file4.close()
text_file5.close()
