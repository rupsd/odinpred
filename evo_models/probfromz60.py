from numpy import *
import time
import string
from scipy.stats import norm

from scipy.integrate import quad
import scipy.interpolate as interpolate
##from numpy.polynomial.hermite import hermgauss

def read_dba(filnavn):
  from cPickle import load
  fil=open(filnavn,'r')
  dba=load(fil)
  return dba

std2zerrdat='''0.0  1.87365 4.71915
0.1  2.76697 1.83914
0.2  2.66592 1.67254
0.3  2.80460 1.62927
0.4  2.84893 1.55950
0.5  3.01659 1.65633
0.6  3.19727 1.87820
0.7  3.58060 2.19015
0.8  3.99582 2.47217
0.9  4.41776 2.67658
1.0  4.83769 3.02099
1.1  5.23488 3.25355
1.2  5.73597 3.51153
1.3  6.04533 3.73392
1.4  6.28619 4.03466
1.5  6.76000 4.26975
1.6  7.34251 4.63529
1.7  7.26797 5.09968
1.8  7.66240 5.45547
1.9  8.12756 5.75006
2.0  8.66108 5.92507
2.1  8.00264 6.33943
2.2  8.36033 6.82380
2.3  9.57945 7.22904
2.4  10.0182 9.55935'''
##std2zerrdatas=string.split(std2zerrdat,'\n')
evoZe=[1.190337119, 1.09722654, 1.5184969, 1.55149632, 1.61350234, 1.62655598, 1.62872112, 1.65231031, 1.76850308, 1.8980599, 2.16951512,
         2.4577851, 2.80050940, 3.1211849, 3.41051604, 3.50174896, 3.65288783, 4.11117191, 4.67562047, 6.26930760, 5.0117271, 5.00937022, 8.09654326, 8.04418484]


def bisqnorm(z,mode=1,fracidp=0.33333):
    loc1,scal1,skewa1,wei1,loc2,scal2,skewa2=-3.291881970407184, 7.7508542282050312, 6.6246663141506641, 0.20223826159617481, 15.284259283569316, 1.7168652246554941, -3.6335101511642236
    wa=(1.89811878217,2.10188121783)
    ##fr0=wa[0]/(wa[0]+wa[1])
    ##print 'calculating',z,mode
    if mode==1:#disorder
      if z<-5:return bisqnorm(-5,mode,fracidp)
      elif z>12:return bisqnorm(12,mode,fracidp)
      zstn1=(z-loc1)/scal1
      p1=norm.pdf(zstn1)*norm.cdf(skewa1*zstn1)*wei1
      return p1/wa[0]*fracidp
    elif mode==2:
      if z<5:return bisqnorm(5,mode,fracidp)
      elif z>16:return bisqnorm(16,mode,fracidp)
      zstn2=(z-loc2)/scal2
      p2=norm.pdf(zstn2)*norm.cdf(skewa2*zstn2)##*wei2
      return p2/wa[1]*(1-fracidp)

def bisqnorm_exp(z,z0,zerr=4.5,mode=None):
    return bisqnorm(z,mode)*norm.pdf(z,z0,zerr)

def integrated_bisqnorm(z,zerr=4.5):
    pDis=quad(bisqnorm_exp,-5,16.2,args=(z,zerr,1),limit=20,epsrel=0.00001)[0]
    pOrd=quad(bisqnorm_exp,-5,16.2,args=(z,zerr,2),limit=20,epsrel=0.00001)[0]
    return pDis/(pDis+pOrd)

def genzmat():
    std2zerrdatas=string.split(std2zerrdat,'\n')
    return [[eval(x) for x in string.split(lin)] for lin in std2zerrdatas]
matrixZe=genzmat()

def getZerror(zstd,useEvo=True):
    if useEvo:
        lzstd=log(zstd)
	if lzstd<-1.4:return 1.6
	elif zstd>2.7:return 6.5
	else:
	  zind=int((lzstd+2.0)/(3.7/24))
	  return evoZe[zind]
    else:
	if zstd<0.6:return 2.7
	elif zstd>1.5:return min(9,zstd*2.0+4.0)#linear model
	else:
	  zind=int(zstd*10)
	  return matrixZe[zind][1]


if False:
 data=[]
 for i,zerr in enumerate(zerrs):
  probs=[integrated_bisqnorm(z,zerr) for z in zs]
  data.append(probs)
  print 'done',zerr,average(probs)
 print 'TIME:',time.time()-T_0
 save_dba2(data,'probfromzdata33.pck')

fullzs=arange(-5.5,17.1,0.1)


def getprobability(Zave,Zstd,useEvo=True,runpath=''):
    intpz=read_dba(runpath+'probfromzinterpq.pck')
    Zerr=getZerror(Zstd,useEvo)
    return intpz(Zave,Zerr)

if __name__ == '__main__':
  from pylab import *
  p05 = getprobability(fullzs,0.5)
  p01 = getprobability(fullzs,1.0)
  p15 = getprobability(fullzs,2.0)
  p05ne = getprobability(fullzs,0.5,False)
  p01ne = getprobability(fullzs,1.0,False)
  p15ne = getprobability(fullzs,1.9,False)
  plot(fullzs,p05,'b')
  plot(fullzs,p05ne,'r')
  plot(fullzs,p15,'k')
  plot(fullzs,p15ne,'m')
  plot(fullzs,p01ne,'g')
  plot(fullzs,p01,'c')
  show()
