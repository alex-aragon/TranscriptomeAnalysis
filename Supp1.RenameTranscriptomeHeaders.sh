#FINAL VERSION
#Changing sequence IDs

#PBS -q default
#PBS -l nodes=1:ppn=16,mem=32Gb,vmem=32Gb,walltime=400:00:00
#PBS -N TranscriptomeRename
#PBS -V


cd $PBS_O_WORKDIR



module load seqkit/2018



awk '{print $1}' CeratopterisAssembly_FilteredList.txt | sed 's/_i.*//g' | sort -u  > temporalIndex.txt

seqkit split CeratopterisAssembly_FilteredSequences.fasta -s 10000

cp temporalIndex.txt IDIndex_Trinity-to-Custom.txt

tnum=1

for j in CeratopterisAssembly_FilteredSequences.fasta.split/*.fasta
do
  while grep -e "^>TRINITY" -q in $j && read i
  do
      sed -i 's/^'"$i"'/'"$i"' Crichardii-Hnn_ARJA-v1.0_transcript'"$tnum"'/g' IDIndex_Trinity-to-Custom.txt
      sed -i 's/^>'"$i"'_/>Crichardii-Hnn_ARJA-v1.0_transcript'"$tnum"'./g' $j
      tnum=$(($tnum+1))
  done < temporalIndex.txt
  awk '/Crichardii/ {} !/Crichardii/ {print $1}' IDIndex_Trinity-to-Custom.txt > temporalIndex.txt
done

rm temporalIndex.txt

cat CeratopterisAssembly_FilteredSequences.fasta.split/*.fasta > CeratopterisTranscriptome_ARJA-v1.0.fasta

sed -i 's/path=.*//g' CeratopterisTranscriptome_ARJA-v1.0.fasta
