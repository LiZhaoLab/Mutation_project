#!/bin/bash
#SBATCH --mail-type=END
#SBATCH --mail-user=ewitt@rockefeller.edu
#this is happening in /ru-auth/local/home/ewitt/220516_novogene/usftp21.novogene.com/raw_data


#one for TE_denovo
/ru-auth/local/home/lzhao/Data_scratch/witt/witt/cellranger-4.0.0/cellranger count --id=Old3a_TE_denovo_run\
                   --transcriptome=/ru-auth/local/home/lzhao/Data_scratch/witt/witt/singlecell/Dmel_TE_denovo_genome/ \
                   --fastqs=Old3a \
                   --sample=Old3a_CKDL220010126-1a-SI_TT_B8_HN2LGDSX3

/ru-auth/local/home/lzhao/Data_scratch/witt/witt/cellranger-4.0.0/cellranger count --id=Old3b_TE_denovo_run\
                   --transcriptome=/ru-auth/local/home/lzhao/Data_scratch/witt/witt/singlecell/Dmel_TE_denovo_genome/ \
                   --fastqs=Old3b \
                   --sample=Old3b_CKDL220010127-1a-SI_TT_B9_HN2LGDSX3

/ru-auth/local/home/lzhao/Data_scratch/witt/witt/cellranger-4.0.0/cellranger count --id=YoungLib3_TE_denovo_run\
                   --transcriptome=/ru-auth/local/home/lzhao/Data_scratch/witt/witt/singlecell/Dmel_TE_denovo_genome \
                   --fastqs=YoungLib3 \
                   --sample=YoungLib3_CKDL220010125-1a-SI_TT_G8_HN2LTDSX3
                   
                   
#one for normal:

/ru-auth/local/home/lzhao/Data_scratch/witt/witt/cellranger-4.0.0/cellranger count --id=Old3a_run\
                   --transcriptome=/ru-auth/local/home/lzhao/Data_scratch/witt/witt/singlecell/Dmel_genome/ \
                   --fastqs=Old3a \
                   --sample=Old3a_CKDL220010126-1a-SI_TT_B8_HN2LGDSX3

/ru-auth/local/home/lzhao/Data_scratch/witt/witt/cellranger-4.0.0/cellranger count --id=Old3b_run\
                   --transcriptome=/ru-auth/local/home/lzhao/Data_scratch/witt/witt/singlecell/Dmel_genome \
                   --fastqs=Old3b \
                   --sample=Old3b_CKDL220010127-1a-SI_TT_B9_HN2LGDSX3

/ru-auth/local/home/lzhao/Data_scratch/witt/witt/cellranger-4.0.0/cellranger count --id=YoungLib3_run\
                   --transcriptome=/ru-auth/local/home/lzhao/Data_scratch/witt/witt/singlecell/Dmel_genome \
                   --fastqs=YoungLib3\
                   --sample=YoungLib3_CKDL220010125-1a-SI_TT_G8_HN2LTDSX3