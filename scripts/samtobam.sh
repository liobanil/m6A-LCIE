#!/bin/bash --login

#SBATCH -J sam_to_bam
#SBATCH -N 4
#SBATCH -e sam_to_bam.%J.err
#SBATCH -o sam_to_bam.%J.out
#SBATCH --mem=32G
#SBATCH --time=36:00:00


conda activate basictools-env

for file in $(ls *.sam | sed "s,.sam,,"); do 
samtools view -bo ${file}.bam ${file}.sam
done


