#!/bin/bash 
#SBATCH --mail-type=end  
#SBATCH --mail-user=ewitt1093@gmail.com

#get private SNPS

bcftools index Old3a.bcf
bcftools index Old3b.bcf
bcftools index YoungLib3.bcf
bcftools index Old3a_gDNA.bcf
bcftools index Old3b_gDNA.bcf
bcftools index Young3_gDNA.bcf


bcftools isec Old3a.bcf Old3a_gDNA.bcf -p isec
cp isec/0001.vcf Old3a.sc_private.vcf


bcftools isec Old3b.bcf Old3b_gDNA.bcf -p isec
cp isec/0001.vcf Old3b.sc_private.vcf

bcftools isec YoungLib3.bcf Young3_gDNA.bcf -p isec
cp isec/0001.vcf YoungLib3.sc_private.vcf

for i in `ls | grep private.vcf`
do
basename=$(echo $i| sed 's/\.sc_private.vcf//g')
rm $basename.read_names.txt
rm $basename.snps_ref_reads.txt
bcftools view -v snps $i|grep -v "^#" | grep -v Scaffold |grep -v 211000 | grep -v mitochondrion| grep -v rDNA| while read line ; do

run=$(echo ${basename}_run )
position=$(echo $line |awk -v OFS="" '{print $2}')
chr=$(echo $line |awk -v OFS="" '{print $1}')
title=$(echo $line |awk -v OFS="" '{print $1,":", $2,"-",$2}')
alt=$(echo $line |awk -v OFS="" '{print $5}')
ref=$(echo $line |awk -v OFS="" '{print $4}')
#echo $position
#echo $title
#echo $alt

read=$(samtools view -b $run/outs/possorted_genome_bam.bam $title | samtools fillmd -e - /ru-auth/local/home/lzhao/Data_scratch/witt/witt/Genomes/dmel_r6.38_FB2021_01/fasta/dmel-all-chromosome-r6.38.fasta | grep -v "^@"| awk -v pos=$position 'BEGIN {OFS = FS = "\t" } ; {n=split($10,a,"") ; if(a[(pos-$4)+1] != "=" ) print pos,(pos-$4)+1, a[(pos-$4)+1], $1, $4, $10 }'|awk -v alt="$alt" '{if($3==alt) print $4}')

#echo $read
for i in `echo $read`
do
echo $chr $position $ref $alt $i >>$basename.read_names.txt 
done


read2=$(samtools view -b  $run/outs/possorted_genome_bam.bam $title | samtools fillmd -e - /ru-auth/local/home/lzhao/Data_scratch/witt/witt/Genomes/dmel_r6.38_FB2021_01/fasta/dmel-all-chromosome-r6.38.fasta | grep -v "^@"| awk -v pos=$position 'BEGIN {OFS = FS = "\t" } ; {n=split($10,a,"") ; if(a[(pos-$4)+1] != "=" ) print pos,(pos-$4)+1, a[(pos-$4)+1], $1, $4, $10 }'|awk -v ref="$ref" '{if($3==ref) print $4}')
#echo $read2
nread=$(echo $read2 | wc -w  )

echo $chr $position $ref $alt $nread >> $basename.snps_ref_reads.txt

done


rm $basename.read_names_and_barcodes.txt

#samtools view $run/outs/possorted_genome_bam.bam | fgrep -w -f <(awk '{print $5}' $basename.read_names.txt ) | awk '{print $1, $(NF-3)}' >>$basename.read_names_and_barcodes.txt



done

samtools view Old3a_run/outs/possorted_genome_bam.bam | fgrep -w -f <(awk '{print $5}' Old3a.read_names.txt ) | awk '{print $1, $(NF-3)}' >read_annotation/Old3a.read_names_and_barcodes.txt
samtools view Old3b_run/outs/possorted_genome_bam.bam | fgrep -w -f <(awk '{print $5}' read_annotation/Old3b.read_names.txt ) | awk '{print $1, $(NF-3)}' >read_annotation/Oldb.read_names_and_barcodes.txt
samtools view YoungLib3_run/outs/possorted_genome_bam.bam | fgrep -w -f <(awk '{print $5}' read_annotation/YoungLib3.read_names.txt ) | awk '{print $1, $(NF-3)}' >read_annotation/YoungLib3.read_names_and_barcodes.txt




samtools depth Old3a_run/outs/possorted_genome_bam.bam >Old3a_coverage.txt
samtools depth Old3b_run/outs/possorted_genome_bam.bam >Old3b_coverage.txt
samtools depth YoungLib3_run/outs/possorted_genome_bam.bam >YoungLib3_coverage.txt

gzip Old3a_coverage.txt
gzip Old3b_coverage.txt
gzip YoungLib3_coverage.txt
