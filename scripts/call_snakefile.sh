#!/bin/bash --login

#SBATCH -J snakefile
#SBATCH -N 20
#SBATCH -e logs/bash/snakefile.%J.err
#SBATCH -o logs/bash/snakefile.%J.out
#SBATCH --mem=64G
#SBATCH --time=96:00:00

#activate conda enviroment
conda activate snakemake

#snakemake call
snakemake --snakefile snakefile --cores 20 --use-conda
