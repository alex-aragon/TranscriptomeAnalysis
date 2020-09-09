#Transcriptome clusterization to reduce the number of sequences for further analysis.
#Two different software were used sequencially: (1) CD-HIT, based on sequence similarity; (2) Compacta, based on sharing same reads mapped to different sequences.

#PBS -q default
#PBS -l nodes=1:ppn=16,mem=32Gb,vmem=32Gb,walltime=100:00:00
#PBS -N Clustering
#PBS -V

cd $PBS_O_WORKDIR

module load cd-hit/4.6
module load bowtie2/2.3.4.2
module load samtools/1.3.1
module load Compacta/1.01





#First clustering step with CD-HIT

cd-hit-est -i CeratopterisRawAssembly.fasta -o CeratopterisAssembly_ClusteringCDHIT.fasta -c 0.95 -n 10 -r 1 -g 1 -M 200000 -T 16 -d 40





#Mapping to generate Compacta inputs

bowtie2-build CeratopterisAssembly_ClusteringCDHIT.fasta IndexCeratopteris_T25bw2

while read i;
do
  bowtie2 -p 12 -q --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant -k 200 \
  -x IndexCeratopteris_T25bw2 \
  -1 ../../Processing/PairedSamples/${i}_FW_PR.fastq.gz \
  -2 ../../Processing/PairedSamples/${i}_RV_PR.fastq.gz \
  -S ${i}_CriT25_bw2Output \
  2>${i}_bw2CriT25_postCDHIT_align_stats.txt
  samtools view -S -b ${i}_CriT25_bw2Output > ${i}_CriT25_bw2Output.bam
done < samples.txt

rm *_CriT25_bw2Output





#Second clustering step with Compacta

Compacta -b SAMPLE1_CriT25_bw2Output.bam,\
SAMPLE2_CriT25_bw2Output.bam,\
SAMPLE3_CriT25_bw2Output.bam,\
SAMPLE4_CriT25_bw2Output.bam,\
SAMPLE5_CriT25_bw2Output.bam,\
SAMPLE6_CriT25_bw2Output.bam,\
SAMPLE7_CriT25_bw2Output.bam,\
SAMPLE8_CriT25_bw2Output.bam,\
SAMPLE9_CriT25_bw2Output.bam \
-n 9 -s SAMPLE1,SAMPLE2,SAMPLE3,SAMPLE4,SAMPLE5,SAMPLE6,SAMPLE7,SAMPLE8,SAMPLE9 \
-g Leaf,Leaf,Leaf,Root,Root,Root,RootTip,RootTip,RootTip \
-t 16 \
-o CeratopterisAssembly_CompactaOutput

rm *_CriT25_bw2Output.bam





#Filtering unique IDs per Compacta clusters

module load seqkit/2018

seqkit sort -2 --quiet -j 4 CeratopterisAssembly_ClusteringCDHIT.fasta > CeratopterisAssembly_ClusteringCDHIT_sorted.fasta

awk '{if ($3 ~ /1/) {print $1}}' CeratopterisAssembly_Compacta2Output_clusters.txt | sort -n > CeratopterisAssembly_FilteredList.txt





#Getting unique transcripts from CDHIT output using IDS from Compacta list

module load seqtk/1.0

seqtk subseq CeratopterisAssembly_ClusteringCDHIT_sorted.fasta CeratopterisAssembly_FilteredList.txt > CeratopterisAssembly_FilteredSequences.fasta
