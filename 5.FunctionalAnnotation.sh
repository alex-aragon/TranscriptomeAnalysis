#Perfoming BUSCO for different versions of the transcriptome clustering
#PBS -q default
#PBS -l nodes=1:ppn=16,mem=16Gb,vmem=16Gb,walltime=400:00:00
#PBS -N BUSCO-TranscriptomeVersions
#PBS -V

cd $PBS_O_WORKDIR

module load ncbi-blast+/2.2.31
module load hmmer/3.1b2
module load BUSCO/2.0.1

while read i
do
  python3  /data/software/busco/BUSCO.py -i ../../1.Clustering/CeratopterisRawAssembly.fasta -o RawAssembly-${i}_BUSCO -l ../../BUSCOAnalysis/${i} -m tran -c 16 -sp arabidopsis

  python3  /data/software/busco/BUSCO.py -i ../../1.Clustering/CeratopterisAssembly_ClusteringCDHIT.fasta -o CDHIT-Filtering-${i}_BUSCO -l ../../BUSCOAnalysis/${i} -m tran -c 16 -sp arabidopsis

  python3  /data/software/busco/BUSCO.py -i ../../1.Clustering/CeratopterisAssembly_FilteredSequences.fasta -o Compacta-Filtering-${i}_BUSCO -l ../../BUSCOAnalysis/${i} -m tran -c 16 -sp arabidopsis
done < ../../BUSCOAnalysis/database.txt











#PBS -q default
#PBS -l nodes=1:ppn=16,mem=16Gb,vmem=16Gb,walltime=400:00:00
#PBS -N CDSPrediction-TranscriptomeVersions
#PBS -V

cd $PBS_O_WORKDIR

module load ncbi-blast+/2.2.31
module load TransDecoder/5.3.0


TransDecoder.LongOrfs -t ../../1.Clustering/CeratopterisRawAssembly.fasta
TransDecoder.Predict -t ../../1.Clustering/CeratopterisRawAssembly.fasta

TransDecoder.LongOrfs -t ../../1.Clustering/CeratopterisAssembly_ClusteringCDHIT.fasta
TransDecoder.Predict -t ../../1.Clustering/CeratopterisAssembly_ClusteringCDHIT.fasta

TransDecoder.LongOrfs -t ../../1.Clustering/CeratopterisTranscriptome_ARJA-v1.0.fasta
TransDecoder.Predict -t ../../1.Clustering/CeratopterisTranscriptome_ARJA-v1.0.fasta

mv CeratopterisTranscriptome_ARJA-v1.0.fasta.transdecoder.cds CeratopterisTranscriptome_CDS_ARJA-v1.0.fasta
mv CeratopterisTranscriptome_ARJA-v1.0.fasta.transdecoder.pep CeratopterisTranscriptome_Prot_ARJA-v1.0.fasta









#Annotation with different tools and putting all in a Trinotate file.
#PBS -q default
#PBS -l nodes=1:ppn=16,mem=16Gb,vmem=16Gb,walltime=400:00:00
#PBS -N PreparingTrinotateInputs
#PBS -V



cd $PBS_O_WORKDIR




module load trinity/2.4.0
module load Trinotate/3.2.1
module load seqkit/2018
module load ncbi-blast+/2.2.31
module load hmmer/3.1b2
module load SignalP/4.1
module load tmhmm/2.0c





#Copy and split peptides prediction file
mkdir ProtPrediction_Split
cp CeratopterisTranscriptome_Prot_ARJA-v1.0.fasta ProtPrediction_Split/
seqkit split ProtPrediction_Split/CeratopterisTranscriptome_Prot_ARJA-v1.0.fasta -s 1000

#Copy and split CDS prediction file
mkdir Transcriptome_Split
cp CeratopterisTranscriptome_ARJA-v1.0.fasta Transcriptome_Split/
seqkit split Transcriptome_Split/CeratopterisTranscriptome_ARJA-v1.0.fasta -s 1000







gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
export PFAMDB="Pfam-A.hmm"
export Viridiplantae="3.4.AnnotationDBs/Viridiplantae_AllReviewed_SwissProt.fasta"
export Arabidopsis_RefProteome="3.4.AnnotationDBs/Arabidopsis_ReferenceProteome_Uniprot.fasta"
export Arabidopsis_Araport11="3.4.AnnotationDBs/Arabidopsis-thaliana_Proteome.fasta"
export Azolla_FernBase="3.4.AnnotationDBs/Azolla-filiculoides_Proteome.fasta"
export Salvinia_FernBase="3.4.AnnotationDBs/Salvinia-cucullata_Proteome.fasta"

#Run BLASTp for peptides prediction file
for i in Proteome_Split/*.fasta
do
  blastp\
    -db $Viridiplantae\
    -query $i\
    -num_threads 16\
    -max_target_seqs 1\
    -outfmt 6\
    -evalue 1e-5\
    -out ${i}_Viridiplantae.out
done
cat Proteome_Split/*_Viridiplantae.out > CeratopterisProteome_BLASTptoUniprot.out

#Run BLASTx for transcriptome file
for i in Transcriptome_Split/*.fasta
do
  blastx\
    -db $Viridiplantae\
    -query $i\
    -num_threads 16\
    -max_target_seqs 1\
    -outfmt 6\
    -evalue 1e-5\
    -out ${i}_Viridiplantae.out
done
cat Transcriptome_Split/*_Viridiplantae.out > CeratopterisTranscriptome_BLASTxtoUniprot.out


#Run HMMER for peptides prediction file
hmmscan --cpu 16 --domtblout CeratopterisProt_HMMERScantoPfam.out  $PFAMDB CeratopterisTranscriptome_PredictedProteome_ARJA-v1.0.fasta > pfam.log



#Run of other tools for annotation
signalp -f short -n signalp.out CeratopterisTranscriptome_PredictedProteome_ARJA-v1.0.fasta

tmhmm --short < CeratopterisTranscriptome_PredictedProteome_ARJA-v1.0.fasta > tmhmm.out

../../../../../storage/data/software/Trinotate-3.0.2/util/rnammer_support/RnammerTranscriptome.pl --transcriptome ../1.Clustering/CeratopterisTranscriptome_ARJA-v1.0.fasta --path_to_rnammer /LUSTRE/usuario/jaragon/CeratopterisAnalysis/Post-Assembly/3.Annotation/rnammer

for i in Transcriptome_Split/*.fasta
do
  ../../../../../storage/data/software/Trinotate-3.0.2/util/rnammer_support/RnammerTranscriptome.pl --transcriptome $i --path_to_rnammer /LUSTRE/usuario/jaragon/CeratopterisAnalysis/Post-Assembly/3.Annotation/rnammer
done
cat *.fasta.rnammer.gff > CeratopterisTranscriptome.rnammer.gff



#Run of extra databases for annotation
tRNAscan-SE -E --max CeratopterisTranscriptome_ARJA-v1.0.fasta \
-o CeratopterisTranscriptome_tRNAPrediction \
-f CeratopterisTranscriptome_tRNA2ndStructures \
-b CeratopterisTranscriptome_tRNAOutput --thread 16

awk '/^Crichardii/ {print $1,"tRNA"$2"-"$5$6,"100","100","0","0",$3,$4,"1","100","0.0",$9}'  CeratopterisTranscriptome_tRNAPrediction > CeratopterisTranscriptome_tRNAs

sed -i 's/ /\t/g' CeratopterisTranscriptome_tRNAs



for i in Proteome_Split/*.fasta
do
  blastp\
    -db $Arabidopsis_Araport11\
    -query $i\
    -num_threads 16\
    -max_target_seqs 1\
    -outfmt 6\
    -evalue 1e-5\
    -out ${i}_ArabidopsisAraport11.out
done

cat Proteome_Split/*_ArabidopsisAraport11.out > CeratopterisProteome_BLASTptoArabidopsis.out



for i in Proteome_Split/*.fasta
do
  blastp\
    -db $Azolla_FernBase\
    -query $i\
    -num_threads 16\
    -max_target_seqs 1\
    -outfmt 6\
    -evalue 1e-5\
    -out ${i}_AzollaFernBase.out
done

cat Proteome_Split/*_AzollaFernBase.out > CeratopterisProteome_BLASTptoAzolla.out



for i in Proteome_Split/*.fasta
do
  blastp\
    -db $Salvinia_FernBase\
    -query $i\
    -num_threads 16\
    -max_target_seqs 1\
    -outfmt 6\
    -evalue 1e-5\
    -out ${i}_SalviniaFernBase.out
done

cat Proteome_Split/*_SalviniaFernBase.out > CeratopterisProteome_BLASTptoSalvinia.out





#Constructing the gene and transcripts matrix
awk '/^>Crichardii/ {print $1,$1}' CeratopterisTranscriptome_ARJA-v1.0.fasta | sed 's/_i.* /\t/g'  > CeratopterisTranscriptome_ARJA-v1.0_gene_trans_map



#Initial Run of Trinotate
/data/software/Trinotate-3.0.2/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate

#Run of Trinotate - BasicData
Trinotate Trinotate.sqlite init --gene_trans_map CeratopterisTranscriptome_gene_trans_map \
--transcript_fasta ../CeratopterisTranscriptome_ARJA-v1.0.fasta \
--transdecoder_pep ../CeratopterisTranscriptome_PredictedProteome_ARJA-v1.0.fasta



#Run of Trinotate - Part 2 (Uniprot & Pfam data)
Trinotate Trinotate.sqlite LOAD_swissprot_blastp CeratopterisProteome_BLASTptoUniprot.out

Trinotate Trinotate.sqlite LOAD_swissprot_blastx CeratopterisTranscriptome_BLASTxtoUniprot.out



Trinotate Trinotate.sqlite LOAD_pfam CeratopterisProt_HMMERScantoPfam.out

Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out

Trinotate Trinotate.sqlite LOAD_signalp signalp.out

Trinotate Trinotate.sqlite LOAD_rnammer CeratopterisTranscriptome.rnammer.gff



#Run of Trinotate - Extra data

Trinotate Trinotate.sqlite LOAD_custom_blast \
--outfmt6 CeratopterisProteome_BLASTptoArabidopsis.out \
--prog blastp --dbtype Arabidopsis_Araport11

Trinotate Trinotate.sqlite LOAD_custom_blast \
--outfmt6 CeratopterisProteome_BLASTptoAzolla.out \
--prog blastp --dbtype Azolla_FernBase

Trinotate Trinotate.sqlite LOAD_custom_blast \
--outfmt6 CeratopterisProteome_BLASTptoSalvinia.out \
--prog blastp --dbtype Salvinia_FernBase

Trinotate Trinotate.sqlite LOAD_custom_blast \
--outfmt6 CeratopterisTranscriptome_tRNAs \
--prog blastx --dbtype tRNAsPrediction



#Generate Trinotate report
Trinotate Trinotate.sqlite report > CeratopterisTranscriptomeAnnotation.xls
