#!/usr/bin/python
import sys 
import subprocess

#!/usr/bin/python
import sys
import subprocess
from subprocess import call
ct=0
dicto=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X']
win=25
win2=int((win-1)/2)
seq=str(sys.argv[1])
text_file=open('total.txt','w')
for s in dicto:
  ct=ct+1
  fp1=open(seq,'r')
  for line1 in fp1:
      line1=line1.replace('\n','')
      line1=line1.replace(' ','')
      text_file2= open("AA.txt", "w")
      line0=''.join(['X'*win2,line1,'X'*win2])
      f=  [line0[i:i+win] for i in range(len(line0)-win+1)]
      for t in f:
               if s in t :
                      number = t.count(s)
                      text_file2.write("%s \n" %(number))

               else:
                      number=0
                      text_file2.write("%s \n" %(number))

  text_file2.close()
  if ct==1:
    subprocess.call("mv AA.txt total.txt", shell=True)
  else:
    subprocess.call("paste total.txt AA.txt > total2.txt && mv total2.txt total.txt", shell=True)
text_file.close()
subprocess.call(" mv total.txt totalAA.txt", shell=True)

