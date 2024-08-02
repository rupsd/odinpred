import argparse
import string
from scipy.optimize import curve_fit
from scipy.special import erfc
import numpy as np
##import matplotlib.pyplot as plt
from numpy import zeros
from itertools import compress, count


def k(ion, T, eps):
    return np.sqrt((8.0 * np.pi * ion * 0.2003457539870666) / (eps * T * 0.0019872041))

def W1(r, ion, T, eps):
    kappa = k(ion, T, eps)
    x = kappa * r / np.sqrt(6.0)
    return 332.286 * np.sqrt(6.0 / np.pi) * (1 - np.sqrt(np.pi) * x * np.exp(x ** 2) * erfc(x)) / (eps * r)   ### main formula

def id2seq(id,pref='../evolution/seqlib'):
    f=open('%s/%s.fasta'%(pref,id))
    buf=[lin for lin in f]
    seq=''
    for i in range(1,len(buf)):
      if i<len(buf)-1:seq+=buf[i][:-1]
      else:seq+=buf[i]
    return seq


def electrostatic_potential(seq, doplot=False, ion=0.0001, temp=298, eps=83.83, gca=5.0, gcb=7.5):

  q0 = {"n": 1.0, "D": -1.0, "E": -1.0, "H": 1.0, "K": 1.0, "R": 1.0, "c": -1.0}

  ##pos = np.array([i for i in xrange(len(seq)) if seq[i] in q0.keys()])  #all indices with acidic or basic residues.
  pos = np.array([i for i in range(len(seq)) if seq[i] in q0.keys()])  #all indices with acidic or basic residues.

  charges=zeros(len(seq))
  for i in pos:
    charges[i]=q0[seq[i]]

  electrostats=[]

  for position,residue in enumerate(seq):
                l1 = np.array(abs(position-pos))    ##distances
               
                q=np.array(charges[pos])

                indices = [i for i, x in enumerate(l1) if x==0]
                l= [s for j, s in enumerate(l1) if j not in indices]
                #l=l1
                #q1=q
                q1=[i for n, i in enumerate(q) if n not in indices]
                d = (gca + np.sqrt(l) * gcb)
                tmp = (sum(W1(d, ion, 298, eps)/q1))
                electrostats.append(tmp)
                ##filename1.write("%s \n" %(tmp)) 
               
  if doplot:
    plt.plot(electrostats)
    plt.axis([1, len(seq), 3,-3.5])

    plt.plot(charges,'ro')
    plt.show()        

  return electrostats 

