cp $2/$1.fasta /home/server/programs/RaptorX_Property_Fast-master/
cd /home/server/programs/RaptorX_Property_Fast-master/
./oneline_command.sh $1.fasta tempout 1 0
rm $1.fasta
