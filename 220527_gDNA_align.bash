#!/bin/bash 
#SBATCH --mail-type=end  
#SBATCH --job-name=gDNA
#SBATCH --mail-user=ewitt@rockefeller.edu

 


bwa mem /ru-auth/local/home/lzhao/Data_scratch/witt/witt/Genomes/dmel_r6.38_FB2021_01/fasta/dmel-all-chromosome-r6.38.fasta /ru-auth/local/home/lzhao/DataforEvan/DNAseq_forscRNA/raw_data/Old3DNA/Old3DNA_CSFP220022597-1a_HN2MFDSX3_L1_1.fq.gz /ru-auth/local/home/lzhao/DataforEvan/DNAseq_forscRNA/raw_data/Old3DNA/Old3DNA_CSFP220022597-1a_HN2MFDSX3_L1_2.fq.gz | samtools sort -O BAM -o Old3a_gDNA.bam


bwa mem /ru-auth/local/home/lzhao/Data_scratch/witt/witt/Genomes/dmel_r6.38_FB2021_01/fasta/dmel-all-chromosome-r6.38.fasta /ru-auth/local/home/lzhao/DataforEvan/DNAseq_forscRNA/raw_data/Old4DNA/Old4DNA_CSFP220022598-1a_HN2MFDSX3_L1_1.fq.gz /ru-auth/local/home/lzhao/DataforEvan/DNAseq_forscRNA/raw_data/Old4DNA/Old4DNA_CSFP220022598-1a_HN2MFDSX3_L1_2.fq.gz | samtools sort -O BAM -o Old3b_gDNA.bam


bwa mem /ru-auth/local/home/lzhao/Data_scratch/witt/witt/Genomes/dmel_r6.38_FB2021_01/fasta/dmel-all-chromosome-r6.38.fasta /ru-auth/local/home/lzhao/DataforEvan/DNAseq_forscRNA/raw_data/Young3DNA/Young3DNA_CSFP220022596-1a_HN2MFDSX3_L1_1.fq.gz /ru-auth/local/home/lzhao/DataforEvan/DNAseq_forscRNA/raw_data/Young3DNA/Young3DNA_CSFP220022596-1a_HN2MFDSX3_L1_2.fq.gz | samtools sort -O BAM -o Young3_gDNA.bam





bcftools mpileup  -Q 25 -f /ru-auth/local/home/lzhao/Data_scratch/witt/witt/Genomes/dmel_r6.38_FB2021_01/fasta/dmel-all-chromosome-r6.38.fasta Old3a_gDNA.bam  | bcftools call -mv -Ob -o Old3a_gDNA.bcf

bcftools mpileup  -Q 25 -f /ru-auth/local/home/lzhao/Data_scratch/witt/witt/Genomes/dmel_r6.38_FB2021_01/fasta/dmel-all-chromosome-r6.38.fasta Old3b_gDNA.bam  | bcftools call -mv -Ob -o Old3b_gDNA.bcf

bcftools mpileup  -Q 25 -f /ru-auth/local/home/lzhao/Data_scratch/witt/witt/Genomes/dmel_r6.38_FB2021_01/fasta/dmel-all-chromosome-r6.38.fasta Young3_gDNA.bam  | bcftools call -mv -Ob -o Young3_gDNA.bcf

