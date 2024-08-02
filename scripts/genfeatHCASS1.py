##from pylab import *
import os
import sys
from numpy import *

accdelta=4
ss3delta=7

from elecpot4 import electrostatic_potential
##from genfeatureshcamsa2 import fracdom


def winave(s,win):
  N=len(s)
  if win>1:
    t=sum([s[delta:N-win+1+delta] for delta in range(win)],axis=0)/win
  else:t=s
  half=(win-1)/2
  half=int(half)
  return hstack((ones(half)*t[0],t,ones(half)*t[-1]))

def gettotcharge(segm):
  cK,cR,cD,cE=segm.count('K'),segm.count('R'),segm.count('D'),segm.count('E')
  return abs(cK+cR-cD-cE),(0.0+cK+cR+cD+cE)/len(segm)

def main(bmrid,numseqs=1,evo=False,doplot=False):
  os.system('mv /home/server/programs/RaptorX_Property_Fast-master/tempout/%s .'%bmrid)

  f=open('%s/%s.all'%(bmrid,bmrid),'r')
  bufa=[lin[:-1] for lin in f]
  seq=bufa[1]
  ss3=bufa[2]
  ss8=bufa[3]
  acc=bufa[4]

  os.system('grep domain %s.hca | grep -v "#">tempdom'%bmrid)
  os.system('grep cluster %s.hca | grep -v "#">tempclu'%bmrid)
  os.system('grep ">" %s.hca | grep -v "#">tempsum'%bmrid)
  doms=[]
  f=open('tempdom','r')
  for lin in f:
    tup=lin.split()
    doms.append(tup[1:])
  clus=[]
  f=open('tempclu','r')
  for lin in f:
    tup=lin.split()
    clus.append(tup[1:])
  f=open('tempsum','r')
  for lin in f:
    tup=lin.split()
    seqlen=int(tup[1])
  print(bmrid,seqlen,len(doms),len(clus),doms,clus)

  nan=0.04
  F=zeros((seqlen,10))
  for dom in doms:
    first=int(dom[0])-1
    final=int(dom[1])
    domlen=final-first
    pval=eval(dom[2])
    ##scor=eval(dom[3])
    scor=eval(dom[3])+2.5
    segm=seq[first:final]
    totcharge,chargedens=gettotcharge(segm)
    print(domlen,pval,scor,totcharge/len(segm),chargedens,segm)
    F[first:final,0]=log(domlen)
    F[first:final,1]=log(0.05/pval)*2.0
    F[first:final,2]=scor
    F[first:final,5]=len(segm)/10.0/(1+totcharge)
    F[first:final,6]=len(segm)/10.0/(1+chargedens)

  cludat=[]
  for clu in clus:
    first=int(clu[0])-1
    final=int(clu[1])##-1
    clulen=final-first
    if clulen>1:
     cludat.append((first,final))
     binclu=clu[2]
     numhph=binclu.count('1')
     segm=seq[first:final]
     totcharge,chargedens=gettotcharge(segm)
     print('cluster:',clulen,numhph,totcharge/len(segm),chargedens,segm,binclu)
     F[first:final,3]=clulen
     F[first:final,4]=numhph*10.0/clulen
     F[first:final,7]=len(segm)/3.0/(1+totcharge)
     F[first:final,8]=len(segm)/3.0/(1+chargedens)
  print(cludat)
  cludat=[(-99,-99)]+cludat+[(999999,999999)]
  dists=[(cludat[i][0]-cludat[i-1][1],-cludat[i][1]+cludat[i+1][0]) for i in range(1,len(cludat)-1)]
  mids=[(cludat[i][1]+cludat[i][0])/2 for i in range(1,len(cludat)-1)]
  mids=[int(mx) for mx in mids]
  print(mids)
  dfull=ones(len(seq))*30
  for i,dat in enumerate(cludat[1:len(cludat)-2]):
    if i>0:
      dfull[mids[i-1]:mids[i]]  =dists[i][0]
    if i<len(cludat)-2:
      dfull[mids[i]:  mids[i+1]]=dists[i][1]
  wdf=winave(dfull,13)
  F[:,9]=50.0/wdf

  #RaptorX sec struct preds

  f=open('%s/%s.ss8'%(bmrid,bmrid),'r')
  ##probabilities are in the order of H G I E B T S L(loops), the 8 secondary structure types used in DSSP 
  buf=[lin.split() for lin in f]
  pss8=[[eval(x) for x in tup[3:]] for tup in buf[2:]]
  P8=array(pss8) 
  E8=log(P8)*P8
  E8[P8<0.000001]=0
  E=sum(E8,axis=1).reshape((seqlen,1)) #entropy in p-column
  print(E.shape)

  ##probabilities are in the order of H E C, the 3 secondary structure types used in DSSP 
  f=open('%s/%s.ss3'%(bmrid,bmrid),'r')
  buf=[lin.split() for lin in f]
  pss3=[[eval(x) for x in tup[3:]] for tup in buf[2:]]
  P3=array(pss3)
  P3Ws=transpose(array([winave(P3[:,I],ss3delta*2+1) for I in range(3)]))
  pc=P3[:,2]
  pcw=winave(pc,ss3delta*2+1)

  C=zeros((len(seq),3))
  start=0
  for i in range(len(seq)):
    if ss8[i] in 'LS':
      C[start:i,0]=log(i-start)*2
      start=i
  start=0
  for i in range(len(seq)):
    if ss8[i] not in 'LS':
      C[start:i,1]=6-log(i-start)*2
      start=i
  start=0
  for i in range(len(seq)):
    if acc[i] in 'BM':
      C[start:i,2]=6-log(i-start)*2
      start=i

  accseq='E'*accdelta+acc+'E'*accdelta

  accentr=[]
  for i in range(len(seq)):
    segm=accseq[i-accdelta:i+accdelta+1]
    cnts=array([segm.count(x) for x in 'EBM'])
    frc=cnts[cnts>0]/(1.0+accdelta*2)
    accentr.append(exp(-sum(log(frc)*frc)))
  EA=array(accentr).reshape((seqlen,1))
  PC=array(pcw).reshape((seqlen,1))

  f=open('%s/%s.acc'%(bmrid,bmrid),'r')
  #probabilities are in the order of B (Bury, pACC: 0-10), M (Medium, pACC: 11-40) and E (Exposed, pACC: 41-100), 
  #where pACC is the relative solvent accessibility value calculated by DSSP 
  buf=[lin.split() for lin in f]
  pacc=[[eval(x) for x in tup[3:]] for tup in buf[3:]]
  PA=array(pacc)
  PAWs=transpose(array([winave(PA[:,I],accdelta*2+1) for I in range(3)]))
 
  ##F=hstack((F,(1-P8[:,7:])*10,(1-PC)*10,exp(-E)*2,EA*2,C)) #adding prob of 'L'
  G=hstack((F,P8*10,10*(1-P3Ws),exp(-E)*2,EA*2,C,PAWs*10))

  epot=electrostatic_potential(seq)
  print(len(seq),len(epot))
  ##print(epot)
  if evo:FD=fracdom(bmrid,seqlen=len(seq),numalns=numseqs)
  else:FD=zeros((seqlen,0))

  if doplot:
   subplot(411)
   title(bmrid)
   imshow(rot90(G),origin='lower',interpolation='none',cmap=cm.hot_r,vmin=0,vmax=13)
   subplot(412)
   ##title(bmrid)
   for i in range(len(G[0])):
     plot(range(len(G)),G[:,i])
   axis((0,len(G),0,20))
   subplot(414)
   plot(range(len(seq)),abs(array(epot)))
   axis((0,len(G),0,4))
   subplot(413)
   for i in range(len(FD[0])):
     plot(range(len(FD)),FD[:,i])
   axis((0,len(G),0,1))
   show()

  EP=array(epot).reshape((seqlen,1))
  F=hstack((F,P8*10,10*(1-P3Ws),exp(-E)*2,EA*2,C,PAWs*10,EP,FD))
  print(F.shape)

  fout=open('featHCASS%s.txt'%bmrid,'w')
  for j in range(len(F)):
    fmt='%8.5f '*len(F[0])
    fmt=fmt[:-1]+'\n'
    fout.write(fmt%tuple(F[j]))
  fout.close()
  if not doplot:return
  else:show()

def writefastasingle(seq,id,libname='seqlib'):
    singfasta=open('%s/%s.fasta'%(libname,id),'w')
    singfasta.write("> %s\n"%id)
    for i in range((len(seq)-1)/80+1):
        singfasta.write(seq[i*80:min((i+1)*80,len(seq))]+'\n')

def runHCAandSS(bmrid,runpath):
        print('running pyHCA for id',bmrid)
        os.system('/home/server/programs/pyHCA/pyHCA/bin/hcatk segment -i %s.fasta -m domain -o %s.hca'%(bmrid,bmrid))
        print('running RaptorX for id',bmrid)
        os.system('/home/server/programs/odinpred3/scripts/runRaptorX.sh %s %s'%(bmrid,runpath))
        os.system('mv /home/server/programs/RaptorX_Property_Fast-master/tempout/%s .'%bmrid)

def copyfasta(entry):
    ##os.system('cp /home/au122487/seqlib/%s.fasta .'%entry)
    os.system('cp seqlib/%s.fasta .'%entry)

def fasta2seq(entry):
    fil=open('%s.fasta'%entry,'r')
    seq=''
    for lin in fil:
        if lin[0]!='>':
            seq+=lin.strip('\n')
    return seq

def runHCAmulti(bmrid):
        buf=[lin for lin in open('filedump/msa%s.txt'%bmrid,'r')]
        numseqs=eval(buf[0].split()[0])
        fout=open('mfas%s.txt'%bmrid,'w')
        for i in range(numseqs):
            mtup=buf[1+i].split()
            seqi,dbid,sco=mtup
            fout.write('>%s\n'%dbid)
            fout.write('%s\n'%seqi)
        fout.close()
        os.system('/home/server/programs/pyHCA/pyHCA/bin/hcatk segment -i mfas%s.txt -m domain -o %sm.hca'%(bmrid,bmrid))
        os.system('grep domain %sm.hca | grep -v "#" > msadoms%s.txt'%(bmrid,bmrid))
        return numseqs



def fracdom(bmrid,seqlen=None,numalns=None,doplot=False):
  ##f=open('tempdom','r')
  f=open('msadoms%s.txt'%bmrid,'r')
  doms=[]
  for lin in f:
    tup=lin.split()
    doms.append(tup[1:])
  ##print(bmrid,len(doms))
  nan=0.04
  firsts=[int(dom[0])-1 for dom in doms]
  finals=[int(dom[1])   for dom in doms]
  pvals =[eval(dom[2])  for dom in doms]
  scors =[eval(dom[3])  for dom in doms]
  if seqlen==None:seqlen=max(finals)
  
  F=zeros((seqlen,2))
  for i,pval in enumerate(pvals):
    first=firsts[i]
    final=finals[i]
    scor=scors[i]
    domlen=final-first
    ##print(domlen,pval,scor)
    if True:
      F[first:final,0]+=1
    if pval<0.03:
      F[first:final,1]+=1

  F*=1.0
  if numalns==None:numalns=F.max() 
  F/=numalns

  if doplot:
   for i in range(len(F[0])):
    plot(range(len(F)),F[:,i])
   axis((0,len(F),0,1))
   show()
  return F

def cleanup(entry):
    os.system('rm -r %s'%entry)
    os.system('rm %s.*'%entry)
    os.system('rm %sm.hca'%entry)
    os.system('rm m*%s.txt'%entry)

def run(evo=True):
   numseqs=1
   entry=sys.argv[1]
   runpath=sys.argv[2]
   if len(sys.argv)>3:
       evoflag=sys.argv[3]
       if evoflag in ("n","no","-n"):evo=False
   ##copyfasta(entry)
   ##seq=fasta2seq(entry)
   runHCAandSS(entry,runpath)
   if evo:numseqs=runHCAmulti(entry)
   else:numseqs=1
   main(entry,numseqs,evo=evo,doplot=False)
   cleanup(entry)

run()
