from numpy import *
import time
import string
from scipy.stats import norm
##from jakob_util import *

from scipy.integrate import quad
import scipy.interpolate as interpolate
##from numpy.polynomial.hermite import hermgauss

def read_dba(filnavn):
  from cPickle import load
  fil=open(filnavn,'r')
  dba=load(fil)
  return dba

std2zerrdat='''0.0  1.87365 0.0
0.1  2.76697 2.26713
0.2  2.66592 2.30841
0.3  2.80460 2.18114
0.4  2.84893 2.13970
0.5  3.01659 2.22641
0.6  3.19727 2.51455
0.7  3.58060 2.96549
0.8  3.99582 3.50521
0.9  4.41776 3.88958
1.0  4.83769 4.33756
1.1  5.23488 4.72794
1.2  5.73597 5.15091
1.3  6.04533 5.46166
1.4  6.28619 5.73674
1.5  6.76000 5.98994
1.6  7.34251 6.65581
1.7  7.26797 6.79843
1.8  7.66240 7.17653
1.9  8.12756 7.92629
2.0  8.66108 8.52250
2.1  8.00264 8.81416
2.2  8.36033 8.74740
2.3  9.57945 8.58681
2.4  10.0182 9.94675'''
std2zerrdatas=string.split(std2zerrdat,'\n')

evoZe=[1.190337119, 1.09722654, 1.5184969, 1.55149632, 1.61350234, 1.62655598, 1.62872112, 1.65231031, 1.76850308, 1.8980599, 2.16951512,
         2.4577851, 2.80050940, 3.1211849, 3.41051604, 3.50174896, 3.65288783, 4.11117191, 4.67562047, 6.26930760, 5.0117271, 5.00937022, 8.09654326, 8.04418484]

def bisqnorm(z,mode=1,fracidp=0.33333):
    ##loc1,scal1,skewa1,wei1,loc2,scal2,skewa2=-3.291881970407184, 7.7508542282050312, 6.6246663141506641, 0.20223826159617481, 15.284259283569316, 1.7168652246554941, -3.6335101511642236
    ##wa=(1.89811878217,2.10188121783)
    loc1,scal1,skewa1,wei1,loc2,scal2,skewa2=-2.5735848312114489, 6.8936499377226141, 7.0651939067336889, 0.18898238375553189, 15.237901315036559, 2.2769128730764292, -4.9021600247562018
    wa=(1.45046161584, 2.54953838416)
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
 T_0=time.time()
 zs=hstack(([-5.5,-4],arange(-3,14,0.25),[15,16,17]))
 zerrs=hstack((arange(2.0,7.1,0.25),[7.5,8.0,9.0]))
 data=[]
 for i,zerr in enumerate(zerrs):
  probs=[integrated_bisqnorm(z,zerr) for z in zs]
  data.append(probs)
  print 'done',zerr,average(probs)
 save_dba2(data,'probfromzdata33_1325.pck')
 print 'TIME:',time.time()-T_0
 A=array(data)
 intpz = interpolate.interp2d(zs,zerrs, A, kind='cubic')
 intpz = interpolate.interp2d(zs,zerrs, A, kind='quintic')
 save_dba2(intpz,'probfromzinterpq_1325.pck')
 save_dba2(intpz,'probfromzinterpc_1325.pck')

fullzs=arange(-5.5,17.1,0.1)

def getprobability(Zave,Zstd,useEvo=True,runpath=''):
    intpz=read_dba(runpath+'probfromzinterpq_1325.pck')
    Zerr=getZerror(Zstd,useEvo)
    return intpz(Zave,Zerr)

def getprobability2(Zave,useEvo=True):
    if useEvo:Zerr=4.16
    else:     Zerr=4.46
    return intpz(Zave,Zerr)

if __name__ == '__main__':
  from pylab import *
  p05 = getprobability(fullzs,0.5)
  p01 = getprobability(fullzs,1.0)
  p15 = getprobability(fullzs,2.0)
  p05ne = getprobability(fullzs,0.5,False)
  p01ne = getprobability(fullzs,1.0,False)
  p15ne = getprobability(fullzs,1.9,False)
  title('Disorder prob vs. Z for different Zstds (1325)')
  plot(fullzs,p05,'b')
  plot(fullzs,p05ne,'r')
  plot(fullzs,p15,'k')
  plot(fullzs,p15ne,'m')
  plot(fullzs,p01ne,'g')
  plot(fullzs,p01,'c')
  show()

#this is 1325 data
