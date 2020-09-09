#PBS -q default
#PBS -l nodes=1:ppn=16,mem=16Gb,vmem=32Gb,walltime=100:00:00
#PBS -N TrimmingRNASeqSamples
#PBS -V

cd $PBS_O_WORKDIR

module load java/1.8
module load Trimmomatic/0.32
module load FastQC/0.11.2
module load MultiQC

for (( i = 2; i < 10; ++i ));
do
  java -jar /LUSTRE/storage/data/software/Trimmomatic-0.32/trimmomatic-0.32.jar PE SAMPLE${i}_FW.fastq SAMPLE${i}_RV.fastq SAMPLE${i}_FW_PR.fastq SAMPLE${i}_FW_UP.fastq SAMPLE${i}_RV_PR.fastq SAMPLE${i}_RV_UP.fastq ILLUMINACLIP:TruSeq2_PE_Alex.fa:2:30:7 SLIDINGWINDOW:4:15 MINLEN:60
  fastqc SAMPLE${i}_*_*.fastq
  multiqc SAMPLE${i}_*_*.zip -o SAMPLE${i}
done
