#!/bin/bash --login

#SBATCH -J index_bams
#SBATCH -N 20
#SBATCH -e logs/bash/index_bams.%J.err
#SBATCH -o logs/bash/index_bams.%J.out
#SBATCH --mem=64G
#SBATCH --time=72:00:00

#activate conda enviroment
conda activate basictools-env

#make .bai
for file in $(ls mapped_data/*/*.bam); do 
samtools index -b $file -@ 20
done
