#Assembly with trinity/2.4.0 with Trinity default settings
#PBS -q ensam
#PBS -l nodes=1:ppn=20,mem=500gb,vmem=750gb,walltime=2000:00:00
#PBS -N Transcriptome_CriSamples_Assembly2
#PBS -V

cd $PBS_O_WORKDIR

module load java/1.8
module load samtools/1.9
module load jellyfish/2.2.10
module load Salmon/0.11.3
module load bowtie2/2.3.4.2
module load trinity/2.4.0

Trinity --seqType fq --samples_file samples.txt --CPU 20 --max_memory 750G --no_normalize_reads --output CeratopterisAssembly_trinity2.4.0 --full_cleanup
