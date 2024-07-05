#!/bin/bash
#SBATCH  -t 2000 
#SBATCH --ntasks-per-core 1 
#SBATCH -A XXXX  
#SBATCH --mem-per-cpu 10G 
#SBATCH -c 8
#SBATCH -e job-snakemake.err 
#SBATCH -o job-snakemake.out 
module load bioinfo-tools
module load pysam/0.17.0-python3.9.5
module load fgbio/2.1.0
module load snakemake
module load bwa-mem2/2.2.1-20211213-edc703f
module load samtools/1.17
module load BEDTools/2.29.2
module load VarDictJava/1.8.3
module load vep/111.0
module load bcftools/1.19

snakemake --verbose -d .  --cores 8 --rerun-incomplete

