#!/bin/sh

## script to run 1st step of REGENIE

#SBATCH --job-name=step1
#SBATCH --partition=compute
#SBATCH --account=your-slurm-account
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 120
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm_logs/%x_%j.out

## ============================================================
## Configuration: set paths for your HPC environment
## ============================================================
export proj_dir="/path/to/your/gwas/directory"
export genotype_dir="/path/to/ukb/pruned_genotypes"      # pruned UKB genotypes for REGENIE step 1
export input_dir="${proj_dir}/input"
export regenie="/path/to/regenie"                          # REGENIE binary
export sample_inclusion="path/to/sample_inclusion_list"
## ============================================================

## change to project directory
cd ${proj_dir}

## run REGENIE
${regenie} \
--step 1 \
--bed ${genotype_dir}/ukb22418_allChrs.pruned \
--extract ${genotype_dir}/ukb22418_allChrs.prune.in \
--keep ${sample_inclusion} \
--phenoFile ${input_dir}/phenotypes.txt \
--covarFile ${input_dir}/covariates.txt \
--threads 110 \
--qt \
--bsize 1000 \
--lowmem \
--lowmem-prefix tmpdir/regenie_tmp_preds \
--out input/ukb_step1_pruned
