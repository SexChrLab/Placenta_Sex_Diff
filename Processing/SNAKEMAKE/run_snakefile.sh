#!/bin/bash
#SBATCH -n 6
#SBATCH --job-name=GZ_fastq
#SBATCH -t 1-0:0
#SBATCH -o slurm.GZ_fastq.out
#SBATCH -e slurm.GZ_fastq.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kcolney@asu.edu
#--------------------------

snakemake --snakefile Snakefile -j 30 --nolock --latency-wait 15 --rerun-incomplete --cluster "sbatch -n 1 --nodes 1 -c 8 -t 04:00:00"


