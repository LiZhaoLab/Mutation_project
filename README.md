# Mutation_project
This repository will support our upcoming paper on the transcriptional and mutational consequences of aging in the male Drosophila germline.

The SNP database is located at https://www.dropbox.com/s/7t4b75h5uvh6lot/210603_R517_Old_Young_SNP_database.csv?dl=0 and seurat RDS for the main analysis is located at https://www.dropbox.com/s/pf1x781c7poubrp/220607_6reps.RDS?dl=0 .
Seurat RDS for the alternate TE/de novo gene reference is located at https://www.dropbox.com/s/algc8fvobs5oxhi/220607_6reps_te_denovo.RDS?dl=0

Raw sequence files for sc-RNA-seq and gDNA-seq will be released on SRA.

220527_gDNA_align.bash shows the steps we used to align gDNA to the reference genome for SNP calling.
220601_Cellranger.bash aligns the single cell data
220614_read_annotation.bash gets the cell barcode information from mutated read
220531_sc.R goes through the single cell analysis

210830_Aging_paper_code.Rmd is an older version of the analysis.
