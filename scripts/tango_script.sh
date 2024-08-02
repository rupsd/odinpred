#!/bin/bash
count=1
file="$1"


for i in `cat $1` ; 

do
mod=${i//[X]/A}

tango_x86_64_release $count'tango' ct="N" nt="N" ph="7" te="298.15" io="0.1000" seq=$mod
count=`expr $count + 1`
done
echo $count
sed -i '1,1d' *tango.txt

count=`expr $count - 1`

for i in `seq 1 $count`

do
  
  awk '{print $3}' $i'tango.txt'> $i'beta.txt'
  awk '{print $4}' $i'tango.txt' > $i'turn.txt'
  awk '{print $5}' $i'tango.txt' > $i'helix.txt'
  awk '{print $6}' $i'tango.txt' > $i'polyproline.txt'
  i=`expr $i + 1`
done

printf "%s\0" 1beta.txt | sort -zn | xargs -0 cat > beta_oneline.txt 
printf "%s\0" 1turn.txt | sort -zn | xargs -0 cat > turn_oneline.txt 
printf "%s\0" 1helix.txt | sort -zn | xargs -0 cat > helix_oneline.txt
printf "%s\0" 1polyproline.txt | sort -zn | xargs -0 cat > polyproline_oneline.txt 
paste beta_oneline.txt helix_oneline.txt turn_oneline.txt  polyproline_oneline.txt > totaltango.txt
rm beta_oneline.txt helix_oneline.txt turn_oneline.txt polyproline_oneline.txt 1beta.txt 1turn.txt 1helix.txt 1polyproline.txt 1tango.txt
