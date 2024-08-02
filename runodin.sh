#!/bin/bash

#######################################
# Usage
#######################################
usage() {
   echo "#######################"
   echo "Usage: ./runodin.sh arg1 arg2 arg3"
   echo "arg1: a text file containing sequences: if the file contains both fasta and simple sequences , leave a space in between"
   echo "arg2: a folder name in which you need results"
   echo "arg3: -y if you want evolution parameter."
   echo "eg. ./runodin.sh input_seq.txt result_folder -y"
   echo "#####################"
   exit 0
}


########################################
# Setup
########################################
setup_variables() {
   # Odin
   odin_dir=/scripts
   scripts_dir=/scripts
   
   # Evo models
   non_evo_model_dir=$odin_dir/scripts/non_evo_models
   evo_model_dir=$odin_dir/scripts/evo_models
   
   # Current/scratch/results directories
   cwd=$PWD
   scratch=$cwd/temp_$1/
   results=$cwd/$2/

   if [ "$3" == "-y" ]; then
      evo=true
   else
      evo=false
   fi
}

########################################
# Run OdinPred on input
########################################
run_odin_pred() {
   mypython=python3
   ppython=python3

   idx=$2
   protein=$3

  	$mypython $scripts_dir./preprocessor.py $protein
  	mv 'input_seq2.txt' 'tk'$1
  	filename='tk'$1
  	$mypython $scripts_dir/AA.py  $filename
   $mypython $scripts_dir motifs.py $filename
   $mypython $scripts_dir./net_charge.py $filename
   $mypython $scripts_dir./length.py $filename
   $mypython $scripts_dir./repeats_improved.py $filename >>/dev/null
   $scripts_dir./tango_script.sh $filename $scripts_dir >>/dev/null
  	$mypython $scripts_dir./eisenberg.py $filename
   $mypython $scripts_dir./flexibility.py $filename
   $mypython $scripts_dir./entropy.py $filename >>/dev/null
   $mypython $scripts_dir./ip.py $filename 
   $mypython $scripts_dir./compositionnew.py $filename >>/dev/null
   $mypython $scripts_dir./chou_fasmann3_features.py $filename >>/dev/null
  	$mypython $scripts_dir./contact_maps_iupred.py $filename $scratch $scripts_dir >>/dev/null
  	seqx=$(head -n 1 $filename)
  	echo '>'$idx>> $idx.fasta 
  	echo $seqx>> $idx.fasta
  	echo $'\n' >> $idx.fasta
   if [ "$evo" = true ]; then
      echo "Entering evolution"
  	   mkdir filedump
  	   $mypython $scripts_dir./evolall3.py $idx $scratch >& log.txt
      if [ ! -e filedump/evolustatsunity$idx.txt ];then
         echo 'No alignments found. Predicting without evolution' > $results/message.txt
  	      ##paste totalAA.txt motifs*.txt net_charge.txt length.txt repeatslist20.txt repeatslist20?.txt totaltango.txt eis_*.txt flexibility.txt entropy.txt ip.txt comptransnew.txt Chou.txt energy_25_again.txt featHCASS$idx.txt > total$1
         $mypython $scripts_dir./genfeatHCASS1.py $idx $scratch -n > loghcass.txt
  	      paste totalAA.txt energy_25_again.txt entropy.txt comptransnew.txt eis_*.txt repeatslist20.txt repeatslist20?.txt flexibility.txt ip.txt length.txt motifs???.txt Chou.txt net_charge.txt featHCASS$idx.txt totaltango.txt > total$1
        ##$mypython $non_evo_model_dir/./evaluations.py $idx total$1 $non_evo_model_dir $scratch
        $ppython $non_evo_model_dir/./evaluations.py $idx total$1 $non_evo_model_dir $scratch
      else
         echo 'Alignments found. Predicting with evolution' > $results/message.txt
  	      ##cp filedump/evolustatsunity$idx.txt $results
         $mypython $scripts_dir./genfeatHCASS1.py $idx $scratch> loghcass.txt
  	      cd filedump
  	      $mypython $scripts_dir./extract.py $idx evolustatsunity$idx.txt 
  	      $scripts_dir./evolution2.sh
  	      mv evolution_matrix.txt ../
  	      rm *
  	      cd $scratch
  	      paste totalAA.txt energy_25_again.txt entropy.txt comptransnew.txt eis_*.txt repeatslist20.txt repeatslist20?.txt flexibility.txt ip.txt length.txt motifs???.txt Chou.txt net_charge.txt featHCASS$idx.txt totaltango.txt evolution_matrix.txt > total$1
         ##$mypython $evo_model_dir/./evaluations.py $idx total$1 $evo_model_dir $scratch
         $ppython $evo_model_dir/./evaluations.py $idx total$1 $evo_model_dir $scratch
      fi
  	else
         $mypython $scripts_dir./genfeatHCASS1.py $idx $scratch -n > loghcass.txt
  	      paste totalAA.txt energy_25_again.txt entropy.txt comptransnew.txt eis_*.txt repeatslist20.txt repeatslist20?.txt flexibility.txt ip.txt length.txt motifs???.txt Chou.txt net_charge.txt featHCASS$idx.txt totaltango.txt > total$1

  	   ##$mypython $non_evo_model_dir/./evaluations.py $idx total$1 $non_evo_model_dir $scratch
  	   $ppython $non_evo_model_dir/./evaluations.py $idx total$1 $non_evo_model_dir $scratch
  	fi

  	paste pred_output.txt final_predictions.txt > DisorderPredictions$idx.txt
  	sed -i '1i\'"Residue  No. Zscore Disorder-Probability" DisorderPredictions$idx.txt
  	mv DisorderPredictions$idx.txt $results/
  	mv DisorderPredictions$idx.pdf $results/
  	##cp total$1 $results/allfeatures$idx.txt
}

########################################
# Read input file
########################################
read_input_file() {
   i=0 # for naming in non-fasta format
   already_read=false
   while $already_read == true || IFS='' read -r line || [[ -n "$line" ]]; do
      already_read=false
      
      # Fix-up line, removing newline and trimming spaces.
      line=`echo $line | tr '\n' ' '`
      line=`echo $line | tr -d '[:space:]'`
      
      if [[ -z "$line" ]] || [[ ${line:0:1} == "#" ]]; then
         # Empty line or comment
         continue
      elif [[ ${line:0:1} == ">" ]]; then
         # Fasta format
         id=${line:1}
         protein=""

         # If fasta, we read until ">" or blank line
         while $already_read == true || IFS='' read -r line || [[ -n "$line" ]]; do
            if [[ ${line:0:1} == ">" ]]; then
               already_read=true
               break
            elif [[ -z "$line" ]] || [[ ${line:0:1} == "#" ]]; then
               break
            fi

            protein=$protein$line
         done
      else
         # If not fasta, each line is a protein
         i=$(($i + 1))
         id=$i
         protein=$line
      fi

      # Fix ID
      id=`echo $id | sed 's/\s.*$//'`  # Remove everything after first space
      id=`echo $id | sed -e 's/|/_/g'` # Change '|' to '_'

      # Fix Protein
      protein=`echo $protein | sed -e 's/ //g'`
      
      # Print what we found
      echo "ID:      "$id
      echo "PROTEIN: "$protein

      # Run Odin Pred on input
      ##if [ -e "$id.done" ]; then
      ##   echo "ID '"$id"' is already in use."
      ##else
         echo "Running ID '"$id"'."
         run_odin_pred $1 $id $protein
         touch $id.done
      ##fi
      
      echo ""
   done < $1 
}

  
########################################
# Main
########################################
# Check commandline input
if [ $# -eq 0 ]; then
   echo "No arguments supplied."
   echo "Type -h for help"
   echo ""
   usage
fi

if [ "$1" == "-h" ] ; then
   usage
fi

if [ -z "$2" ]; then
   second=result
else
   second=$2
fi

# Setup
echo "Setting up directories"
setup_variables $1 $second $3

rm -rf $scratch
mkdir -p $scratch
mkdir -p $results
cp $1 $scratch
cd $scratch

# Read input and run odin
read_input_file $1
