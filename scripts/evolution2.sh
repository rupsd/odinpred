#!/bin/bash
>totalxyz.txt
for (( i=0; i <= 28; ++i ))
do
 paste 'totalxyz.txt' 'evolv_featuresunity_training'$i'.txt' > 'total2.txt'
 mv 'total2.txt' 'totalxyz.txt'
 rm 'evolv_featuresunity_training'$i'.txt'
done


#for (( i=4; i <= 28; ++i ))
#do
# paste 'totalxyz.txt' 'evolv_featuresunity_training'$i'.txt' > 'total2.txt#'
 #mv 'total2.txt' 'totalxyz.txt'
 #rm 'evolv_featuresunity_training'$i'.txt'
#d#one
mv totalxyz.txt evolution_matrix.txt

