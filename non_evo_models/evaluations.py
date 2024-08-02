
from __future__ import division
from tensorflow import keras
import tensorflow as tf
##from tensorflow.python import keras
##import keras
from keras.models import load_model
##tf.keras.backend.clear_session()
##import keras.models
from keras.initializers import glorot_uniform
from keras.utils import CustomObjectScope
import numpy as np
import sys
from numpy import genfromtxt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import os
sys.path.append(os.path.join(os.path.dirname(sys.argv[3]), os.path.dirname(sys.argv[4])))
##from probfromz5 import getprobability
from probfromz6 import getprobability

 
f=sys.argv[3]
a=str(sys.argv[1])
b=sys.argv[2]
text_file2= open('final_predictions.txt', "w")
ps=[]
pred=np.array([])
prob_predictions=np.array([])
tf.keras.backend.clear_session()
with CustomObjectScope({'GlorotUniform': glorot_uniform()}):
 for d in range(1,11):
 ##for d in range(10):
   x=np.genfromtxt(b)
   z=np.array(x)
   ##print(z.shape)
   x_val2=np.nan_to_num(z)
   ##print(x_val2.shape)
   meaners2=np.loadtxt(f+'/'+'meaners'+str(d)+'.txt', delimiter=',' ,unpack=True)
   ##print(meaners2.shape)
   stder2=np.loadtxt(f+'/'+'stders'+str(d)+'.txt', delimiter=',' ,unpack=True)
   ##new_model = load_model(f+'/'+'my_model_adam'+str(d)+'.h1')
   x_val2= (x_val2 - meaners2) / stder2
   ##new_model = load_model(f+'/'+'JT_best_model'+str(d)+'nn128reg0.01nn180nn220nn315nn410.h5')
   ##new_model = keras.models.load_model(f+'/'+'JT_best_model'+str(d)+'nn128reg0.01nn180nn220nn315nn410.h5')
   new_model = keras.models.load_model(f+'/'+'model_nonevo'+str(d)+'.h5')
   classes6 = new_model.predict(x_val2)
   classes6[classes6>1116.5]=1116.5
   classes6[classes6<-5]=-5
   pred=np.hstack([pred, classes6]) if pred.size else classes6
Zave=np.mean(pred,axis=1)
Zave=Zave.reshape(-1, 1)
Zstd=np.std(pred,axis=1)
Zstd=Zstd.reshape(-1, 1)
for h,k in zip(Zave,Zstd):
 prob_predictions=getprobability(h,k,useEvo=True,runpath=f+'/')
 ps.append(prob_predictions)
for l in range(len(ps)):
      c=l+1
      text_file2.write("%3s %7.4f  %6.4f\n" %(c,Zave[l],ps[l]))
aps=np.array(ps)
dis=aps>0.5
disofrac=np.average(dis)*100.0
print('Disordered residues: %6.1f'%disofrac+"%") 
text_file2.close()
fig=plt.figure()
plt.plot(ps)
plt.xlim((1, x_val2.shape[0]))
plt.xlabel('Residue number')
plt.ylim((0, 1))
plt.ylabel('Probability of disorder')
plt.title('Disorder predictions for (no evolution)'+' ' +str(a))
plt.savefig('DisorderPredictions'+str(a)+'.pdf')
