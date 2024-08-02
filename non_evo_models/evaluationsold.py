#/usr/bin/python
from __future__ import division
from tensorflow import keras
import numpy as np
from numpy import genfromtxt
#import matplotlib.pyplot as plt
import subprocess

def Average(lst): 
    return mean(lst) 


pred=np.array([])
prob_predictions=np.array([])
for d in range(1,11):
   x=np.genfromtxt('total.txt', delimiter='\t' ,unpack=True)
   z=np.array(x)
   x_val2=np.nan_to_num(z)
   print (x_val2.shape)
   meaners2=np.loadtxt('meaners'+str(d)+'.txt', delimiter=',' ,unpack=True)
   meaners2=np.reshape(meaners2,(len(meaners2),1))
   stder2=np.loadtxt('stders'+str(d)+'.txt', delimiter=',' ,unpack=True)
   stder2=np.reshape(stder2,(len(meaners2),1))
   x_val2= (x_val2 - meaners2) / stder2
   x_val2=np.transpose(x_val2)
   new_model = keras.models.load_model('my_model_adam'+str(d)+'.h1')
   classes6 = new_model.predict(x_val2)
   classes6[classes6>15]=15
   classes6[classes6<-5]=5
   #plt.plot(classes6)
   pred=np.hstack([pred, classes6]) if pred.size else classes6
final_predictions=np.mean(pred,axis=1)
final_predictions=final_predictions.reshape((-1, 1))
Zaverage=Average(final_predictions)
output=np.hstack([final_predictions,prob_predictions]) 
print output
np.savetxt('final_predictions.txt',output, delimiter=',')
#plt.plot(final_predictions, 'r-')
#plt.xlim((1, x_val2.shape[0]))
#plt.xlabel('residue no.')
#plt.ylim((-5, 15))
#plt.ylabel('predicted zscores')
#plt.title('ODiNPred')
#plt.savefig("plot.pdf")
#plt.show(block=False)
subprocess.call("paste pred_output.txt final_predictions.txt > total2.txt && rm pred_output.txt", shell=True)




