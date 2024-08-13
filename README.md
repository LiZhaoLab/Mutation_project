# Mutation_project
This repository supports our paper on the transcriptional and mutational consequences of aging in the male Drosophila germline.

[updated 2024 Because Dropbox folder rearragment] For Large files please see:
https://www.dropbox.com/scl/fo/o4bfvxq7e59vn7yovq8gu/AJsCJKZK452f9rvDr6Fm6q4?rlkey=d7chwem9pdxq9gfq7wj3s2v4c&st=v53pkwof&dl=0
220607_6reps.RDS: This is the RDS file. It is a very large file, so we cannot submit to figshare. You can easily generate from the raw data from NCBI. If you have difficulty accessing it, please contact lzhao@rockefeller.edu for a copy. 
210603_R517_Old_Young_SNP_database.csv: This is all the SNP information (not unique lines, all relevant read information )
210603_R517_Old_Young_SNP_database.simp.txt: A simplified file to show unique SNP information from all cells, totaling of ~5000 SNPs from hundreds of thousand germ cells

Raw sequence files for sc-RNA-seq and gDNA-seq are available on SRA with accession number PRJNA777411.

220527_gDNA_align.bash shows the steps we used to align gDNA to the reference genome for SNP calling.
220601_Cellranger.bash aligns the single cell data,
220614_read_annotation.bash gets the cell barcode information from mutated read,
220908_Witt_Old_Young_SC.rmd is our code to produce figures.

