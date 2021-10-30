#!/bin/bash

#SBATCH --job-name='chipimputation'
#SBATCH --cpus-per-task=8
#SBATCH --mem=50GB
#SBATCH --output=chipimput-%j-stdout.log
#SBATCH --error=chipimput-%j-stderr.log 
#SBATCH --time=10-00:00:00

#cd /scratch3/users/nanje/chipimputation

echo "Submitting SLURM job"

nextflow  -c /scratch3/users/nanje/chipimputation/chipimputation/test_imp5.config \
  run /scratch3/users/nanje/chipimputation/chipimputation/main-dsl2_imp5.nf -profile slurm,singularity -resume
