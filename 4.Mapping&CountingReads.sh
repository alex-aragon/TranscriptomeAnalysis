#PBS -q ensam
#PBS -l nodes=1:ppn=10,vmem=60Gb,walltime=100:00:00
#PBS -N TranscriptomeMpSamples
#PBS -V

cd $PBS_O_WORKDIR

module load bowtie2/2.3.4.2
module load samtools/1.9
module load express/1.5.1


#Mapping and counting reads from samples to final version of the transcriptome

bowtie2-build CeratopterisTranscriptome_ARJA-v1.0.fasta IndexCeratopterisTranscriptome

while read i;
do

  bowtie2 -p 10 -q -a -X 600 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 --no-discordant --no-mixed \
  -x IndexCeratopterisTranscriptome \
  -1 ../../Processing/PairedSamples/${i}_FW_PR.fastq.gz \
  -2 ../../Processing/PairedSamples/${i}_RV_PR.fastq.gz \
  -S ${i}_CriTranscriptome-v1.0_bwOutput \
  2>${i}_CriTranscriptome-v1.0_align_stats.txt

  samtools view -S -b ${i}_CriTranscriptome-v1.0_bwOutput > ${i}_CriTranscriptome-v1.0_bwOutput.bam
  samtools sort -n ${i}_CriTranscriptome-v1.0_bwOutput.bam ${i}_CriTranscriptome-v1.0_bwOutput.bam.sorted

  express --calc-covar --no-bias-correct -o ${i}_CriTranscriptome-v1.0_eXpressOutput \
          CeratopterisTranscriptome_ARJA-v1.0.fasta \
          ${i}_CriTranscriptome-v1.0_bwOutput.bam.sorted

done < samples.txt

rm *_CriTranscriptome-v1.0_bwOutput





#Mapping and counting reads from other samples to final version of the transcriptome

while read j;
do
  while read i;
  do

    bowtie2 -p 10 -q -a -X 600 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 --no-discordant --no-mixed \
    -x IndexCeratopterisTranscriptome \
    -1 ../../RawData/${j}/${i}_1_PR.fastq.gz \
    -2 ../../RawData/${j}/${i}_2_PR.fastq.gz \
    -S ${j}.${i}_CriTranscriptome-v1.0_bwOutput \
    2>${j}.${i}_CriTranscriptome-v1.0_align_stats.txt

    samtools view -S -b ${j}.${i}_CriTranscriptome-v1.0_bwOutput > ${j}.${i}_CriTranscriptome-v1.0_bwOutput.bam
    samtools sort -n ${j}.${i}_CriTranscriptome-v1.0_bwOutput.bam ${j}.${i}_CriTranscriptome-v1.0_bwOutput.bam.sorted

    express --calc-covar --no-bias-correct -o ${i}_CriTranscriptome-v1.0_eXpressOutput \
            CeratopterisTranscriptome_ARJA-v1.0.fasta \
            ${j}.${i}_CriTranscriptome-v1.0_bwOutput.bam.sorted


  done < samples_otherdata.txt
done < directories_otherdata.txt

rm *_CriTranscriptome-v1.0_bwOutput
