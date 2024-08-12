# Mutation_project
This repository will support our upcoming paper on the transcriptional and mutational consequences of aging in the male Drosophila germline.

Theseurat RDS for the main analysis is located at [https://www.dropbox.com/s/pf1x781c7poubrp/220607_6reps.RDS?dl=0](https://www.dropbox.com/scl/fi/ommrr3nf1aat54vxz7ak9/220607_6reps.RDS?rlkey=31iyq6zfudbzkkb6m6doq762n&st=4bdh2qzd&dl=0). Seurat RDS for the alternate TE/de novo gene reference is located at [https://www.dropbox.com/scl/fi/ommrr3nf1aat54vxz7ak9/220607_6reps.RDS?rlkey=31iyq6zfudbzkkb6m6doq762n&dl=0](https://www.dropbox.com/scl/fi/ommrr3nf1aat54vxz7ak9/220607_6reps.RDS?rlkey=31iyq6zfudbzkkb6m6doq762n&st=4bdh2qzd&dl=0).
[update: If the links are not availble, please regenerate from the raw reads, or please contact lzhao@rockefeller.edu for a copy.]

Raw sequence files for sc-RNA-seq and gDNA-seq are available on SRA with accession number PRJNA777411.

220527_gDNA_align.bash shows the steps we used to align gDNA to the reference genome for SNP calling.
220601_Cellranger.bash aligns the single cell data,
220614_read_annotation.bash gets the cell barcode information from mutated read,
220908_Witt_Old_Young_SC.rmd is our code to produce figures.

