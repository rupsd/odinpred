from numpy import *
import os
import sys
import time
import numpy
import math
##import matplotlib.mlab as mlab

aa13dict={'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY', 'F': 'PHE', 'I': 'ILE', 'H': 'HIS', 'K': 'LYS', 'M': 'MET', 'L': 'LEU', 'N': 'ASN', 'Q': 'GLN', 'P': 'PRO', 'S': 'SER', 'R': 'ARG', 'T': 'THR', 'W': 'TRP', 'V': 'VAL', 'Y': 'TYR'}
aa31dict={'CYS': 'C', 'GLN': 'Q', 'ILE': 'I', 'SER': 'S', 'VAL': 'V', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'LYS': 'K', 'THR': 'T', 'PHE': 'F', 'ALA': 'A', 'HIS': 'H', 'GLY': 'G', 'ASP': 'D', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'GLU': 'E', 'TYR': 'Y'}
aa3s=aa31dict.keys();aa3s.sort()#introduce ordering
aa1s=aa13dict.keys();aa1s.sort()
aa3s=['ALA','ARG' ,'ASP' ,'ASN' ,'CYS' ,'GLU' ,'GLN' ,'GLY' ,'HIS' ,'ILE' ,'LEU' ,'LYS' ,'MET' ,'PHE' ,'PRO' ,'SER' ,'THR' ,'TRP' ,'TYR' ,'VAL']
aa1s3=[aa31dict[k] for k in aa3s]

#from 9-gram
HighlyConserved=['C'] #0
Hydrophobic=['M'] #1
PolarAcids=['N','Q'] #2
polar=['D','E'] #3
Alcohols=['S','T'] #4
Aliphatic=['P','A','G'] #5
Aliphatic_big=['I','V','L'] #6
Aromatic=['F','Y','W'] #7
Bases=['H','K','R'] #8

dct9={'C':'0', 'M':'1', 'N':'2', 'Q':'2', 'D':'3', 'E':'3', 'S':'4', 'T':'4', 'P':'5', 'A':'5', 'G':'5', 'I':'6', 'V':'6', 'L':'6', 'F':'7', 'Y':'7', 'W':'7', 'H':'8', 'K':'9', 'R':'9', '-':'-', 'X':'-'}

refdat='''Ala  8.25
Arg  5.53
Asn  4.06
Asp  5.45
Cys  1.37
Gln  3.93
Glu  6.75
Gly  7.07
His  2.27
Ile  5.96
Leu  9.66
Lys  5.84
Met  2.42
Phe  3.86
Pro  4.70
Ser  6.56
Thr  5.34
Trp  1.08
Tyr  2.92
Val  6.87'''

def initfil2(filename):
    file=open(filename,'r')
    buffer=file.readlines()
    file.close()
    for i in range(len(buffer)):
       ##buffer[i]=string.split(buffer[i])
       buffer[i]=buffer[i].split()
    return buffer

def initfil(filename):
    file=open(filename,'r')
    buffer=file.readlines()
    file.close()
    return buffer

def get_refcomposition():
    ##buf=initfil2('refcompositionProtscale.txt')
    data=string.split(refdat,'\n')
    buf=[string.split(rd) for rd in data]
    dct={}
    for lin in buf:
	aa1=aa31dict[string.upper(lin[0])]
	dct[aa1]=eval(lin[1])
    return dct

blosum62={
('W', 'F') : 1, ('L', 'R') : -2, ('S', 'P') : -1, ('V', 'T') : 0, 
('Q', 'Q') : 5, ('N', 'A') : -2, ('Z', 'Y') : -2, ('W', 'R') : -3, 
('Q', 'A') : -1, ('S', 'D') : 0, ('H', 'H') : 8, ('S', 'H') : -1, 
('H', 'D') : -1, ('L', 'N') : -3, ('W', 'A') : -3, ('Y', 'M') : -1, 
('G', 'R') : -2, ('Y', 'I') : -1, ('Y', 'E') : -2, ('B', 'Y') : -3, 
('Y', 'A') : -2, ('V', 'D') : -3, ('B', 'S') : 0, ('Y', 'Y') : 7, 
('G', 'N') : 0, ('E', 'C') : -4, ('Y', 'Q') : -1, ('Z', 'Z') : 4, 
('V', 'A') : 0, ('C', 'C') : 9, ('M', 'R') : -1, ('V', 'E') : -2, 
('T', 'N') : 0, ('P', 'P') : 7, ('V', 'I') : 3, ('V', 'S') : -2, 
('Z', 'P') : -1, ('V', 'M') : 1, ('T', 'F') : -2, ('V', 'Q') : -2, 
('K', 'K') : 5, ('P', 'D') : -1, ('I', 'H') : -3, ('I', 'D') : -3, 
('T', 'R') : -1, ('P', 'L') : -3, ('K', 'G') : -2, ('M', 'N') : -2, 
('P', 'H') : -2, ('F', 'Q') : -3, ('Z', 'G') : -2, ('X', 'L') : -1, 
('T', 'M') : -1, ('Z', 'C') : -3, ('X', 'H') : -1, ('D', 'R') : -2, 
('B', 'W') : -4, ('X', 'D') : -1, ('Z', 'K') : 1, ('F', 'A') : -2, 
('Z', 'W') : -3, ('F', 'E') : -3, ('D', 'N') : 1, ('B', 'K') : 0, 
('X', 'X') : -1, ('F', 'I') : 0, ('B', 'G') : -1, ('X', 'T') : 0, 
('F', 'M') : 0, ('B', 'C') : -3, ('Z', 'I') : -3, ('Z', 'V') : -2, 
('S', 'S') : 4, ('L', 'Q') : -2, ('W', 'E') : -3, ('Q', 'R') : 1, 
('N', 'N') : 6, ('W', 'M') : -1, ('Q', 'C') : -3, ('W', 'I') : -3, 
('S', 'C') : -1, ('L', 'A') : -1, ('S', 'G') : 0, ('L', 'E') : -3, 
('W', 'Q') : -2, ('H', 'G') : -2, ('S', 'K') : 0, ('Q', 'N') : 0, 
('N', 'R') : 0, ('H', 'C') : -3, ('Y', 'N') : -2, ('G', 'Q') : -2, 
('Y', 'F') : 3, ('C', 'A') : 0, ('V', 'L') : 1, ('G', 'E') : -2, 
('G', 'A') : 0, ('K', 'R') : 2, ('E', 'D') : 2, ('Y', 'R') : -2, 
('M', 'Q') : 0, ('T', 'I') : -1, ('C', 'D') : -3, ('V', 'F') : -1, 
('T', 'A') : 0, ('T', 'P') : -1, ('B', 'P') : -2, ('T', 'E') : -1, 
('V', 'N') : -3, ('P', 'G') : -2, ('M', 'A') : -1, ('K', 'H') : -1, 
('V', 'R') : -3, ('P', 'C') : -3, ('M', 'E') : -2, ('K', 'L') : -2, 
('V', 'V') : 4, ('M', 'I') : 1, ('T', 'Q') : -1, ('I', 'G') : -4, 
('P', 'K') : -1, ('M', 'M') : 5, ('K', 'D') : -1, ('I', 'C') : -1, 
('Z', 'D') : 1, ('F', 'R') : -3, ('X', 'K') : -1, ('Q', 'D') : 0, 
('X', 'G') : -1, ('Z', 'L') : -3, ('X', 'C') : -2, ('Z', 'H') : 0, 
('B', 'L') : -4, ('B', 'H') : 0, ('F', 'F') : 6, ('X', 'W') : -2, 
('B', 'D') : 4, ('D', 'A') : -2, ('S', 'L') : -2, ('X', 'S') : 0, 
('F', 'N') : -3, ('S', 'R') : -1, ('W', 'D') : -4, ('V', 'Y') : -1, 
('W', 'L') : -2, ('H', 'R') : 0, ('W', 'H') : -2, ('H', 'N') : 1, 
('W', 'T') : -2, ('T', 'T') : 5, ('S', 'F') : -2, ('W', 'P') : -4, 
('L', 'D') : -4, ('B', 'I') : -3, ('L', 'H') : -3, ('S', 'N') : 1, 
('B', 'T') : -1, ('L', 'L') : 4, ('Y', 'K') : -2, ('E', 'Q') : 2, 
('Y', 'G') : -3, ('Z', 'S') : 0, ('Y', 'C') : -2, ('G', 'D') : -1, 
('B', 'V') : -3, ('E', 'A') : -1, ('Y', 'W') : 2, ('E', 'E') : 5, 
('Y', 'S') : -2, ('C', 'N') : -3, ('V', 'C') : -1, ('T', 'H') : -2, 
('P', 'R') : -2, ('V', 'G') : -3, ('T', 'L') : -1, ('V', 'K') : -2, 
('K', 'Q') : 1, ('R', 'A') : -1, ('I', 'R') : -3, ('T', 'D') : -1, 
('P', 'F') : -4, ('I', 'N') : -3, ('K', 'I') : -3, ('M', 'D') : -3, 
('V', 'W') : -3, ('W', 'W') : 11, ('M', 'H') : -2, ('P', 'N') : -2, 
('K', 'A') : -1, ('M', 'L') : 2, ('K', 'E') : 1, ('Z', 'E') : 4, 
('X', 'N') : -1, ('Z', 'A') : -1, ('Z', 'M') : -1, ('X', 'F') : -1, 
('K', 'C') : -3, ('B', 'Q') : 0, ('X', 'B') : -1, ('B', 'M') : -3, 
('F', 'C') : -2, ('Z', 'Q') : 3, ('X', 'Z') : -1, ('F', 'G') : -3, 
('B', 'E') : 1, ('X', 'V') : -1, ('F', 'K') : -3, ('B', 'A') : -2, 
('X', 'R') : -1, ('D', 'D') : 6, ('W', 'G') : -2, ('Z', 'F') : -3, 
('S', 'Q') : 0, ('W', 'C') : -2, ('W', 'K') : -3, ('H', 'Q') : 0, 
('L', 'C') : -1, ('W', 'N') : -4, ('S', 'A') : 1, ('L', 'G') : -4, 
('W', 'S') : -3, ('S', 'E') : 0, ('H', 'E') : 0, ('S', 'I') : -2, 
('H', 'A') : -2, ('S', 'M') : -1, ('Y', 'L') : -1, ('Y', 'H') : 2, 
('Y', 'D') : -3, ('E', 'R') : 0, ('X', 'P') : -2, ('G', 'G') : 6, 
('G', 'C') : -3, ('E', 'N') : 0, ('Y', 'T') : -2, ('Y', 'P') : -3, 
('T', 'K') : -1, ('A', 'A') : 4, ('P', 'Q') : -1, ('T', 'C') : -1, 
('V', 'H') : -3, ('T', 'G') : -2, ('I', 'Q') : -3, ('Z', 'T') : -1, 
('C', 'R') : -3, ('V', 'P') : -2, ('P', 'E') : -1, ('M', 'C') : -1, 
('K', 'N') : 0, ('I', 'I') : 4, ('P', 'A') : -1, ('M', 'G') : -3, 
('T', 'S') : 1, ('I', 'E') : -3, ('P', 'M') : -2, ('M', 'K') : -1, 
('I', 'A') : -1, ('P', 'I') : -3, ('R', 'R') : 5, ('X', 'M') : -1, 
('L', 'I') : 2, ('X', 'I') : -1, ('Z', 'B') : 1, ('X', 'E') : -1, 
('Z', 'N') : 0, ('X', 'A') : 0, ('B', 'R') : -1, ('B', 'N') : 3, 
('F', 'D') : -3, ('X', 'Y') : -1, ('Z', 'R') : 0, ('F', 'H') : -1, 
('B', 'F') : -3, ('F', 'L') : 0, ('X', 'Q') : -1, ('B', 'B') : 4,
('X', 'X') : -1, ('-','-') : 1, ('X','-') : -1, ('-','X') : -1 
}


for aa1 in aa1s:
    blosum62[(aa1,'-')]=-5
    blosum62[(aa1,aa1)]=5
for pair in blosum62.keys():blosum62[(pair[1],pair[0])]=blosum62[pair]

def calc_jsd(p,q):
    #calculate the Jensen-Shannon Divergence for p vs. q
    mix=(p+q)*0.5
    return 0.5*sum(p*log(p/mix))+0.5*sum(q*log(q/mix))

def calcfeatures(bmrID,seq,trimmedlist,doplot=False,dowrite=True,writefull=True,wmode='Henikoff',returnarc=False):
    rcdct=get_refcomposition()
    buf=initfil2(runpath+'filedump/msa%s.txt'%bmrID)
    ##alns=[s[:-1] for s in buf]
    alns=[]
    print buf[0]
    Na=string.atoi(buf[0][0])
    Nr=string.atoi(buf[0][1])
    Ni=string.atoi(buf[0][2])#number of insersions
    dave=eval(buf[0][3])#average depth of phylogenetivc tree
    lnf =eval(buf[0][4])#tree diversity parameter...
    seqbuf=[(buf[1+i][0],buf[1+i][1],eval(buf[1+i][2])) for i in range(Na)]
    seqbuf.sort(lambda x,y: cmp(x[-1],y[-1]))#the reference sequence is the first
    phydists=[]##tup[-1] for tup in seqbuf]
    ##phyweights=1.0/(array(phydists)+0.05)
    print seqbuf[0]
    print seqbuf[-1]
    frins=zeros(Nr)
    for lin in buf[Na+2:]:
	fri=eval(lin[1])#fraction of insertions in alignment at this pos
        if fri<0.00001:
	  frent=(1-fri)*log(1-fri)
        elif fri>0.99999:
	  frent=fri*log(fri)
        else:
	  frent=fri*log(fri)+(1-fri)*log(1-fri)#entropy of 2-state prob
	frins[string.atoi(lin[0])]=-frent#was fri
    seqdct={}
    skipempty=0
    triminds=[]
    for i,tup in enumerate(seqbuf):
       s=tup[0]
       seqid=tup[1]
       treedist=tup[2]
       aln=''
       for x in s:
	if x in aa1s or x=='-':aln+=x
	else:aln+='X'
       if len(set(list(aln))-set(['-']))==0:
	print 'warning: empty sequence!:',bmrID,i
	skipempty+=1
       else:
        alns.append(aln)
        seqdct[seqid]=aln
	phydists.append(treedist)
        if trimmedlist=='all' or seqid in trimmedlist:triminds.append(i-skipempty)
    print len(alns),len(alns[0])
    Na=len(alns)
    phyweights=1.0/(array(phydists)+0.05)
    ##Nr=len(alns[0])
    cols=[[alns[i][j] for i in range(Na)] for j in range(Nr)]
    triminds=array(triminds)

    if False:
     #now calculate the ancestral sequence evolution........
     ancseqs=initfil2(runpath+'filedump/arcseqs%s.txt'%bmrID)[0]
     print ancseqs
     ancscos=zeros(Nr)
     for i in range(len(ancseqs)-1):
	s0=seqdct[ancseqs[i]]
	s1=seqdct[ancseqs[i+1]]
	print s0,s0[50],s1[50],blosum62[(s0[50],s1[50])]
        ancscos+=array([blosum62[(s0[j],s1[j])] for j in range(Nr)])
     ancscos/=(len(ancseqs)-1)
     ancscos/=5.0
     ##print ancscos

    #now calculate the Henikoff weights..................
    mat=[]
    aass=[]
    wHen=ones(Na)#used if wmode=='unity'
    if True:##wmode=='Henikoff': (#we need the aass)
     for j in range(Nr):
      si=set(cols[j])
      aasi=list(si)##;aasi.sort()
      aass.append(aasi)
      cnti=[cols[j].count(aai) for aai in aasi]
      dcti=dict(zip(aasi,cnti))
      wcols=[1.0/dcti[aa] for aa in cols[j]]
      ##print 'tjekv',j,string.join(cols[j],''),wcols
      mat.append(wcols)
     M=array(mat)
     sums=sum(M,axis=0)
     wHen=sums/sum(sums)#henikoff_weights=sums/sum(sums)
     wHenord=argsort(wHen)
    if wmode=='Phylo':
	wHen=phyweights
	print 'using weights from phylogenetic distances',len(wHen),Na
    elif wmode=='unity':wHen=ones(Na)
    ##wHenord=argsort(phyweights)#better?
    ##print wHen
    print len(wHen),Na

    #Stats for the aligned sequences (alns).............
    if True:
     allfracs=[];allaas=[]
     refset=set(rcdct.keys())
     for i in range(Na):
      ##aasi=list(set(list(alns[i])));aasi.sort()
      si=set(list(alns[i]))&refset
      aasi=list(si);aasi.sort()
      cnti=array([alns[i].count(aai) for aai in aasi])# if aai in aa1s])
      fracsi=cnti*1.0/sum(cnti)
      allfracs.append(fracsi)
      allaas.append(aasi)
      #cntgap=alns[i].count('-')
      if i==0:
       fracs0=fracsi;aas0=aasi
       reffracs=array([rcdct[aai] for aai in aasi])/100.0## if aai in aa1s])/100.0
       mix=(fracsi+reffracs)/2
       jsd=0.5*sum(fracsi*log(fracsi/mix))+0.5*sum(reffracs*log(reffracs/mix))#Jensen-Shannon Divergence
       print 'sqJSD0',bmrID,sqrt(jsd),seq
     alljsds=[]
     for i in range(1,Na):
      aasi=allaas[i]
      isin0=array([aa in aas0 for aa in aasi])#Boolean array
      isini=array([aa in aasi for aa in aas0])
      ##print i,aasi,isin0,alns[i]
      fracsi=allfracs[i][isin0]
      fracs0i=fracs0[isini]
      jsdi=calc_jsd(fracsi,fracs0i)
      ##print i,jsdi
      alljsds.append(jsdi)
     print 'averageJSD',bmrID,average(sqrt(alljsds),weights=wHen[1:]),Na,Nr
     if False:
       plot(range(1,Na),alljsds)
       title(bmrID)
       axis([0,2500,0.0,0.1])
       return

    #Derive the features from the columns......................
    ##print string.join(cols[2],'')
    ents=[];jsds=[];fgaps=[];entsstd=[];ents9=[]
    for j in range(Nr):
      ##g9j=[dct9[cj] for cj in cols[j]]
      ##g9j=[dct9[cj] for cj in array(cols[j])[wHenord][:500]]
      g9j=[dct9[cj] for cj in array(cols[j])[triminds]]
      g9js=list(set(g9j))
      ##print(str.join('',g9j))
      fracsj=[average(array(cols[j])==aai,weights=wHen) for aai in aass[j]]## if aai in aa1s]
      ##fracs9=[average(array(g9j)==g9i) for g9i in g9js if g9i!='-']#weights removed in average
      fracs9=[average(array(g9j)==g9i) for g9i in g9js]
      ##print(('fracs9',j,len(fracs9),fracs9))
      fracsjstd=[average(array(cols[j])==aai,weights=wHen) for aai in aass[j] if aai in aa1s]
      entrop=sum(fracsj*log(fracsj))#Shannon Entropy
      entrop9=sum(fracs9*log(fracs9))
      ##print(('etrop9',j,len(fracs9),entrop9))
      entropstd=sum(fracsjstd*log(fracsjstd))
      ents.append(entrop)
      ents9.append(entrop9)
      entsstd.append(entropstd)
      reffracs=array([rcdct[aai] for aai in aass[j] if aai in aa1s])/100.0
      jsd=calc_jsd(fracsjstd,reffracs)
      ##mix=(fracsj+reffracs)/2
      ##jsd=0.5*sum(fracsj*log(fracsj/mix))+0.5*sum(reffracs*log(reffracs/mix))#Jensen-Shannon Divergence
      jsds.append(sqrt(jsd))
      fgaps.append(average(array(cols[j])=='-',weights=wHen))
    cons=[average(array(cols[j])==seq[j],weights=wHen) for j in range(Nr)]
    scos=[average([blosum62[(cols[j][i],seq[j])] for i in range(Na)],weights=wHen) for j in range(Nr)]
    scos=array(scos)/5.0
    E=array(ents)
    E9=array(ents9)
    Estd=array(entsstd)
    ##print cons

    if False:
      plot(range(Nr),scos)
      plot(range(Nr),exp(ents))
      plot(range(Nr),ancscos)
      title(bmrID)
      show();1/0

    #Features for correlated mutations
    T_0=time.time()
    #since list comprehension can produce memory error for size Na*Nr^2 we no longer do this:
    ##pcols=[[alns[i][j]+alns[i][k] for i in range(Na)] for j in range(Nr) for k in range(Nr)]
    pents=[]
    ##for k in range(Nr**2):
    if Na>=500:
        print 'using approx. truncating corr. evol. to 500 trimmed seqs',bmrID,Na
        ##search_is=wHenord[:500]
        search_is=list(triminds)
    else:search_is=range(Na)
    for k in range(Nr):
      for j in range(Nr):
       if j<=k:pents.append(-999)
       else:
        ##pcols=[alns[i][j]+alns[i][k] for i in range(Na)]
        pcols=[dct9[alns[i][j]]+dct9[alns[i][k]] for i in search_is]##range(Na)]
	sk=set(pcols)
	pairs=list(sk)
        pairsnog=[pn for pn in pairs if not '-' in pn]
        ##fracsk=[average(array(pcols)==pk,weights=wHen) for pk in pairs]
        ##fracsk=array([average(array(pcols)==pk,weights=wHen) for pk in pairsnog]) #dont weight the average...
        ##fracsk=array([average(array(pcols)==pk) for pk in pairsnog])
        fracsk=array([average(array(pcols)==pk) for pk in pairs])
        fracsk/=sum(fracsk)
	entrop=sum(fracsk*log(fracsk))
	pents.append(entrop)
        ##print((k,j,len(pairs),pairs,entrop))
        #print((k,j,len(pairs),len(pairsnog),pairsnog,entrop))
    P=array(pents).reshape((Nr,Nr))
    for k in range(Nr):
      for j in range(k):P[k,j]=P[j,k]
    ##C=transpose(P-E)-E
    ##C=transpose(P-Estd)-Estd
    C=transpose(P-E9)-E9
    C[range(Nr),range(Nr)]=0.0
    cents=sum(C,axis=0)/Nr
    print(('TIME:',time.time()-T_0)),sum(cents)
    if False:
      ##imshow(P,origin='lower',interpolation='none',cmap=cm.jet)
      ##imshow(C,origin='lower',interpolation='none',cmap=cm.jet)
      imshow(C,origin='lower',interpolation='none',cmap=cm.RdBu_r)
      ##imshow(exp(-P),origin='lower',interpolation='none',cmap=cm.jet)
      colorbar()
      show();1/0
      
    if doplot:
      fga=array(fgaps)
      zfga=fga==0
      efga=log(fga)*fga
      efga[zfga]=0.0
      plot(range(Nr),cons)
      plot(range(Nr),scos)
      plot(range(Nr),exp(ents))
      plot(range(Nr),cents)
      ##plot(range(Nr),jsds,'k-')
      ##plot(range(Nr),fgaps,'k-')
      plot(range(Nr),efga,'k-')
      title(bmrID)
      show();1/0
    print 'averages',bmrID,average(cons),average(scos),average(exp(ents)),average(jsds),average(cents),dave,lnf,Na,Nr
    consseq=''
    if dowrite:outfile=open(runpath+'filedump/evolustats'+wmode+bmrID+'.txt','w')
    for j in range(Nr):
      if cons[j]>0.8:consseq+=seq[j]
      else:consseq+=string.lower(seq[j])
      if dowrite:
	##outfile.write('%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f'%(cons[j],scos[j],exp(ents[j]),cents[j],jsds[j],fgaps[j]))
	outfile.write('%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f  %6.4f %6.4f %6.4f'%
                      (cons[j],scos[j],exp(ents[j]),cents[j],jsds[j],fgaps[j],frins[j],dave,lnf))
        if writefull:
         fracsj=[average(array(cols[j])==aai,weights=wHen) for aai in aa1s]
	 outfile.write((' %6.3f'*20+'\n')%tuple(fracsj))
        else:
	 if cons[j]<0.5:
	  outfile.write(' 0'*20+'\n')
	 else:
	  ##ci=aa1s.index(seq[j])
	  ci=aa1s.index(seq[j])+1
	  outfile.write(' 0'*(ci-1))#was BUG: ...this gives one too many if ci==0 (aa1=='A')
	  outfile.write(' %7.4f'%cons[j])
	  outfile.write(' 0'*(20-ci)+'\n')
    print 'consseq:',bmrID,consseq
    ##xticks(range(Nr),seq)

def get_refcomp_ss():
    buf=initfil2('refcompositionSS_RefDB_DSSP809.txt')
    fhelix=[eval(lin[1]) for lin in buf]
    fsheet=[eval(lin[2]) for lin in buf]
    return array(fhelix),array(fsheet)#frequencies ordered as in aa1s (and evolstats output final 20 cols)
   
from Bio.SeqIO.FastaIO import *

import sys
import string
import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Align import AlignInfo
from numpy import *

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

##from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import ClustalOmegaCommandline

from Bio.Align.AlignInfo import SummaryInfo

from Bio import AlignIO

from Bio.Phylo.Applications import PhymlCommandline

from Bio import Phylo

from Bio.Blast.Applications import NcbiblastpCommandline

from Bio.SubsMat import FreqTable


def get_aligninfo(hsp,parentseq):
  aln=''
  if hsp.query_start>1:aln+=((hsp.query_start-1)*'-')
  for first in range(len(hsp.query)):
    if hsp.query[first]!='-':break
  firstgaps=len(aln)
  gapcount=0
  for i in range(first,len(hsp.query)):
    if hsp.query[i]=='-':gapcount+=1
    else:
      aln+=hsp.sbjct[i]
      if hsp.sbjct[i]=='-':gapcount+=1#OK??
  endgaps=0
  if len(aln)<len(parentseq):
    print('adding gaps to end of alignment')
    endgaps=len(parentseq)-len(aln)
    aln+=('-'*(len(parentseq)-len(aln)))
  ##return aln
  return aln,firstgaps,gapcount,endgaps
      
def getstats(alns,parentseq):
  print('--------------------------------------------------')
  cons=''
  numcons=0
  allfris=[]
  for i in range(len(alns[0])):
    res=[alns[j][i] for j in range(len(alns))]
    ##fri=res.count(res[0])/len(res)
    aai=parentseq[i]
    fri=res.count(aai)/len(res)
    si=set(res)
    if len(si)==1:
      cons+=aai
      numcons+=1
      allfris.append(1.0)
    else:
      print('nonconserved:',i,fri,si)
      cons+=aai.lower()
      allfris.append(fri)
  fris=array(allfris)
  print(cons)
  avfri=average(fris)
  avfri2=average(log(1-fris+0.01))
  print('consfrac:',bmrID,numcons/len(alns[0]),avfri,avfri2)
    
def print_info_content(summary_info, ftab, fout=None, rep_record=0):
    """3 column output: position, aa in representative sequence, ic_vector value."""
    fout = fout or sys.stdout
    ics=[];ics2=[]
    if not summary_info.ic_vector:
        summary_info.information_content(e_freq_table=ftab,chars_to_ignore=['-','X','Z','U','J','B','O'])
    rep_sequence = summary_info.alignment[rep_record].seq
    for pos, ic in enumerate(summary_info.ic_vector):
        ##fout.write("%d %s %.3f\n" % (pos, rep_sequence[pos], ic))
        ics.append(ic)
        if rep_sequence[pos]!='-':ics2.append(ic)
    return ics,ics2
    
##E_VALUE_THRESH = 0.00001
E_VALUE_THRESH = 0.001

def count_accepted(alignments, threshold=E_VALUE_THRESH):
    best=None
    bestexpect=999
    cnt=0
    seqs=[]
    for aln in alignments:
        for i,hsp in enumerate(aln.hsps):
            if hsp.expect < threshold:
               ##cnt+=1
               if i==0:
                 cnt+=1
                 seqs.append(hsp.sbjct)
               ##print(i,aln.accession,cnt,hsp.expect,hsp.sbjct)
            if hsp.expect < bestexpect:
               best=aln.accession
               bestexpect=hsp.expect
    numdistinct=len(set(seqs))    
    print('number of distinct sequences',numdistinct)
    ##return cnt,best
    return numdistinct,best

def get_seqrecs(alignments, parent, threshold=E_VALUE_THRESH, useonly='all'):
    dct={}
    for aln in alignments:
        ##for hsp in aln.hsps:
        for i,hsp in enumerate(aln.hsps):
            if hsp.expect < threshold:
              if useonly=='all' or aln.accession in useonly:
                ##yield SeqRecord(Seq(hsp.sbjct), id=aln.accession)
                sr=SeqRecord(Seq(hsp.sbjct), id=aln.accession)
                ##dct[hsp.sbjct]=sr
                dct[hsp.sbjct]=hsp.expect,aln.accession,sr
                ##yield SeqRecord(Seq(hsp.sbjct), id=aln.accession+'_'+str(i))
                break #only use the very first hsp
    print('numuniq:',len(dct.values()))
    srp=SeqRecord(Seq(parent), id='parent')
    dct[parent]=0.0,'parent',srp
    ##newbest=min(dct.values())[2].id
    ##return list(dct.values()) #only use entries with unique sbjct sequences!
    return [tup[2] for tup in dct.values()]##,newbest #only use entries with unique sbjct sequences!

def align_best_vs_parrent(alignments,parent,best,bmrID):
    for aln in alignments:
        for hsp in aln.hsps:break
        break
    print(aln.accession,best)
    if aln.accession!=best:return None,0,0,0,None
    print('getting gapped info...')
    bestaln,firstgaps,gapcount,endgaps=get_aligninfo(hsp,parent)
    print(bestaln)
    gapfrac=(firstgaps+gapcount+endgaps)/len(parent)
    print('gapinfo:',bmrID,firstgaps,gapcount,endgaps,gapfrac)
    return bestaln,firstgaps,gapcount,endgaps,gapfrac


def get_blasts(dbname,parent,savename,proc='local'):
 try:b_handle = open(savename, 'r')
 except IOError:
  if proc=='www':
   b_results = NCBIWWW.qblast('blastp',dbname,parent,hitlist_size=2500)#was 50
   bres=b_results.read()
   save_file=open(savename,'w')
   save_file.write(bres)
   save_file.close()
  elif proc=='local':
   cline = NcbiblastpCommandline(cmd="/home/server/programs/ncbi-blast-2.7.1+/bin/blastp",
            query=bmrID+".fasta", db="/home/server/programs/blastdb/swissprot/swissprot",
            outfmt=5, out=savename)
   print(cline)
   cline()
  b_handle = open(savename, 'r')
 b_it = NCBIXML.parse(b_handle)
 for blast in b_it:return blast#return the first (and only)

def get_inserts(alignment,pssm,firstgaps,endgaps,best,bestaln,parent):
 print(alignment)
 lenal=len(alignment)
 print('length of alignment:',lenal)#i.e. the number of sequences (rows in pssm)
 print(alignment[:,17]) 
 ##print(pssm)
 print(bestaln)
 print('gaps:',firstgaps,endgaps)
 for n in range(len(alignment)):
   if alignment[n].id==best:
        refalign=alignment[n]
        nref=n
        print('reference(best) alignment at index:',n)
        print(refalign)
        print(refalign.__class__)
        ##print(dir(refalign))
        print(refalign.seq)
        ##print(refalign[3])
        print(len(refalign)) 
 cnt=0
 gapopen=False
 ins=[]
 for i in range(len(refalign)):
   if refalign[i]=='-':
     if not gapopen:
       go=i
     gapopen=True
   else:
     if gapopen:
       gc=i
       ins.append((go,gc))
     gapopen=False
 finalpos=i
 if gapopen:ins.append((go,-1))
 print(ins)
 if len(ins)==0:
   print('no insertions',bmrID,bestaln)
   return [('-'*firstgaps)+str(alseq.seq)+('-'*endgaps) for alseq in alignment],{}
 ##if len(ins)>0:print(pssm[ins[0][1]])
 if ins[0][0]==0:#alignment starts with an insertion
   start=ins[0][1]-firstgaps
   #can be negative
 ##else:start=0
 else:start=-firstgaps
 if ins[-1][1]==-1:endins=finalpos-ins[-1][0]+1#trailing number of insertions
 else:endins=0
 print('starting from',start,endins)
 pins=[]#position of insersion in refseq
 ##for n in list(range(2))+[nref]:
 allseqs=[]
 for n in range(lenal):
   if ins[0][0]>0:
     if start>=0:seqn=alignment[n][start:ins[0][0]]
     else:seqn=('-'*abs(start))+alignment[n][:ins[0][0]]
   else:
     seqn=('-'*(ins[0][1]-start))
   for gi in range(len(ins)-1):
     if n==nref:pins.append(len(seqn))
     seqn+=alignment[n][ins[gi][1]:ins[gi+1][0]]
   ##if n==nref:pins.append(len(seqn))
   if n==nref:pins.append(len(seqn)-1)
   if ins[-1][1]!=-1:#this means endins==0
    seqn+=alignment[n][ins[-1][1]:]
    if endgaps>0:
       seqn+=('-'*endgaps)
   else:
    if endgaps>0:
     if endgaps<=endins: 
       seqn+=alignment[n][ins[-1][0]:ins[-1][0]+endgaps]
     else:
       print('adding end gaps...',endgaps,endins)
       seqn+=alignment[n][ins[-1][0]:ins[-1][0]+endins]
       seqn+=('-'*(endgaps-endins))
   print(seqn.seq)
   allseqs.append(seqn.seq)
   if n==nref:refseq=seqn
 frdct={}
 print(pins,len(refseq))
 for ig,gi in enumerate(ins):
   if gi[1]==-1:searchrange=range(gi[0],finalpos)
   else:searchrange=range(gi[0],gi[1])
   fracins=1.0-average([alignment[:,j].count('-') for j in searchrange])/lenal
   print(ig,gi,pins[ig])
   print(pins[ig],refseq[pins[ig]-1],refseq[pins[ig]],gi,fracins)
   frdct[pins[ig]]=fracins
 if bestaln!=refseq.seq:
   print('warning differing sequences',bmrID,bestaln,refseq.seq)
 if len(parent)!=len(refseq.seq):
   print('warning differing sequence lengths',bmrID,len(parent),len(refseq.seq))
 return allseqs,frdct
   
def writeresults(bmrID,allseqs,ids,dists,dave,lnf,frdct):
  file=open(runpath+'filedump/msa'+bmrID+'.txt','w')
  file.write('%3d %4d %2d %6.4f %7.4f\n'%(len(allseqs),len(allseqs[0]),len(frdct),dave,lnf))
  for n,seq in enumerate(allseqs):
    file.write('%s %s %6.4f\n'%(seq,ids[n],dists[n]))
  sins=[' ']*len(seq)
  for i in frdct:sins[i]='|'
  file.write(str.join('',sins)+'\n')
  for i in frdct: file.write('%3d %6.4f\n'%(i,frdct[i]))
  file.close()

def findnearest(my_tree,cl,depth,best):
 if cl.name!=None:return depth,0,cl.name
 ##else:for cli in cl:return findnearest(cli) this is depth-first
 else:
   depth+=1
   results=[]
   for i in range(len(cl)):
     cliname=cl[i].name
     if cliname!=None:results.append((depth,my_tree.distance(cliname,best),cliname))
   if len(results)==0:return min([findnearest(my_tree,cli,depth,best) for cli in cl])
   return min(results)

def preparealign(bmrID,parent,runpath):
 print('preparealign for',parent)
 try:
   file=open(runpath+'filedump/msa'+bmrID+'.txt','r')
   return
 except IOError:pass
 phyname=runpath+'filedump/phy'+bmrID+'.phy'
 try:
  open(phyname,'r')
  doalign=False
 except IOError:
  doalign=True
 if True:##doalign:
  savename=runpath+'filedump/alnssp'+bmrID+'.out'
  blast=get_blasts('swissprot',parent,savename)
  numacc,best=count_accepted(blast.alignments, threshold=E_VALUE_THRESH)
  numaln = len(blast.alignments)
  print('number of alignments(swissprot)',bmrID,numaln,numacc)
  if numacc<10:
   print('getting alignment using full nr db',bmrID,numaln,numacc)
   savename=runpath+'filedump/alnsnr'+bmrID+'.out'
   blast=get_blasts('nr',parent,savename,proc='www')
   numacc,best=count_accepted(blast.alignments, threshold=E_VALUE_THRESH)
   numaln = len(blast.alignments)
   print('number of alignments(nonredundant)',bmrID,numaln,numacc)

  if doalign:
   best_seqs = get_seqrecs(blast.alignments,parent)
   best=parent
   SeqIO.write(best_seqs, runpath+'filedump/seqs'+bmrID+'.fasta', 'fasta')
   clustalo_exe="/home/server/programs/clustalo-1.2.4-Ubuntu-x86_64"
   clustalomega_cline = ClustalOmegaCommandline(cmd=clustalo_exe,infile=runpath+"filedump/seqs"+bmrID+".fasta",
       outfile=runpath+"filedump/aligned"+bmrID+".aln", outfmt="clu", force="force")#clustal format
   clustalomega_cline()
   AlignIO.convert(runpath+"filedump/aligned"+bmrID+".aln", "clustal", runpath+"filedump/phy"+bmrID+".phy", "phylip-relaxed")

  ##bestaln,firstgaps,internalgaps,endgaps,gapfrac=align_best_vs_parrent(blast.alignments,parent,best,bmrID)
  bestaln,firstgaps,internalgaps,endgaps,gapfrac='parent',0,0,0,0

  if bestaln==None:
   savename=runpath+'filedump/alnsnr'+bmrID+'.out'
   blast=get_blasts('nr',parent,savename,proc='www')
   bestaln,firstgaps,internalgaps,endgaps,gapfrac=align_best_vs_parrent(blast.alignments,parent,best,bmrID)
  elif gapfrac>0.2 or internalgaps>0:
   savename=runpath+'filedump/alnsnr'+bmrID+'.out'
   blast=get_blasts('nr',parent,savename,proc='www')
   print('redoing alignmentdata',bmrID,gapfrac,internalgaps)
   if doalign:
     best_seqs= get_seqrecs(blast.alignments,parent)
     best=parent
     SeqIO.write(best_seqs, runpath+'filedump/seqs'+bmrID+'.fasta', 'fasta')
     clustalomega_cline = ClustalOmegaCommandline(cmd=clustalo_exe,infile=runpath+"filedump/seqs"+bmrID+".fasta",
       ##outfile="filedump/phy"+bmrID+".phy", outfmt="phy", force="force")
       outfile=runpath+"filedump/aligned"+bmrID+".aln", outfmt="clu", force="force")#clustal format
     clustalomega_cline()
     AlignIO.convert(runpath+"filedump/aligned"+bmrID+".aln", "clustal", runpath+"filedump/phy"+bmrID+".phy", "phylip-relaxed")

 #end-except IOError:
 ##return
 os.system("/home/server/programs/FastTree -out "+runpath+"filedump/phy%s.phy_fasttree.txt %sfiledump/phy%s.phy"%(bmrID,runpath,bmrID))
 test_tree = Phylo.read(runpath+"filedump/phy"+bmrID+".phy_fasttree.txt", "newick")
 #-------------trim tree using Treemmer...--------
 if numacc>=500:
   stepsize=int(numacc/10+(numacc**2)/2500/2500*100)
   print(('trimming tree',bmrID,numacc,stepsize))
   os.system('python2.7 /home/server/programs/odinpred_backupagain/scripts/Treemmer_v0.1_betaJTN.py -X 500 -r %d -np -c 1 -lp 1 %sfiledump/phy%s.phy_fasttree.txt'%(stepsize,runpath,bmrID))
   trimname=runpath+'filedump/phy%s.phy_fasttree.txt_trimmed_list_X_500'%bmrID
   trimmedlist=[lin.strip('\n') for lin in initfil(trimname)]
 else:trimmedlist='all'
 print trimmedlist

 if True:
  Phylo.draw_ascii(test_tree)
  from Bio.Phylo import PhyloXML

  # Promote the basic tree to PhyloXML
  test_phy = test_tree.as_phyloxml()

  # Make a lookup table for sequences
  best_seqs= get_seqrecs(blast.alignments,parent)
  best='parent'
  lookup = dict((rec.id, str(rec.seq)) for rec in best_seqs)
  print('numseqrecs:',len(lookup))

  for clade in test_phy.get_terminals():
    key = clade.name
    accession = PhyloXML.Accession(key, 'NCBI')
    mol_seq = PhyloXML.MolSeq(lookup[key], is_aligned=True)
    sequence = PhyloXML.Sequence(type='aa', accession=accession, mol_seq=mol_seq)
    clade.sequences.append(sequence)

  # Save the annotated phyloXML file
  ##Phylo.write(test_phy,annotname, 'phyloxml')
  ##my_tree=Phylo.read(annotname, 'phyloxml')
  #Phylo.draw_ascii(my_tree)
 ##my_tree=Phylo.read(annotname, 'phyloxml')
 my_tree=test_phy
 totnum=my_tree.count_terminals()
 nodes=my_tree.get_nonterminals()
 nodefrac=len(nodes)*1.0/totnum
 print('number of clades in tree',bmrID,totnum)
 print('number of nodes in tree',bmrID,len(nodes))
 depths=my_tree.depths()
 ##print(depths)
 depthsu= my_tree.depths(unit_branch_lengths=True)
 ##print(depths.values())
 dave=average(list(depths.values()))
 ##print(depthsu.values())
 daveu=average(list(depthsu.values()))
 print('most similar sequence',bmrID,best)
 pclades=my_tree.get_path(best)
 print(len(pclades))
 print(pclades)
 print([cl.name for cl in pclades])
 alldists=[my_tree.distance(clade.name,best) for clade in my_tree.get_terminals()]
 ##print(alldists)
 aveall=average(alldists)

 alignment=AlignIO.read(runpath+"filedump/aligned"+bmrID+".aln", "clustal")
 idsaln=[aln.id for aln in alignment]
 alldistsaln=[my_tree.distance(alnid,best) for alnid in idsaln]
 print(average(alldistsaln),aveall,dave)
 print(len(alldistsaln))
 sinf=SummaryInfo(alignment)
 print(sinf)
 ##pssm=sinf.pos_specific_score_matrix(chars_to_ignore=['X','Z','U','J','B','O'])
 pssm=None
 allseqs,frdct=get_inserts(alignment,pssm,firstgaps,endgaps,best,bestaln,parent)

 cons=sinf.dumb_consensus()
 print(cons)
 gapc=sinf.gap_consensus()
 print(gapc)
 ##AlignInfo.print_info_content(sinf)

 bave=my_tree.total_branch_length()/totnum
 tdist=my_tree.distance(best)#length of path to root
 if False: #not used
   rcdct=get_refcomposition()
   ftab=FreqTable.FreqTable(rcdct,FreqTable.FREQ)
   ics,ics2=print_info_content(sinf,ftab)
   avic=average(ics)
   avic2=average(ics2)#excluding gapped
   ##print('sumstat',bmrID,totnum,dave,daveu,dave/totnum,daveu/totnum,aveall,average(ics))
   print('sumstat %5s %4d %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6.3f %7.4f %7.4f %7.4f'%
          (bmrID,totnum,dave,daveu,dave/totnum,daveu/totnum,aveall,avic,avic2,nodefrac,log(1-nodefrac),bave,tdist))
 writeresults(bmrID,allseqs,idsaln,alldistsaln,dave,log(1-nodefrac),frdct)
 return trimmedlist


bmrID=sys.argv[1]
runpath=sys.argv[2]+'/'
##parent=sys.argv[2]
fbuf=initfil(runpath+bmrID+'.fasta')
parent=string.join([lin[:-1] for lin in fbuf[1:]],'')
trimmedlist=preparealign(bmrID,parent,runpath)
calcfeatures(bmrID,parent,trimmedlist,wmode='unity')


