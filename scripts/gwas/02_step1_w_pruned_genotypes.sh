#!/bin/sh

## script to run 1st step of REGENIE

#SBATCH --job-name=step1
#SBATCH --partition=compute
#SBATCH --account=sc-users
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 120
#SBATCH --mail-type=FAIL
#SBATCH --output=/sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/04_gwas/slurm_logs/%x_%j.out

## export location of genotype files to be used
export dir=/sc-projects/sc-proj-computational-medicine/data/UK_biobank/genotypes/regenie_44448/pruned_genotypes
export input_dir=/sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/04_gwas/input

## change to relevant directory
cd /sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/04_gwas

## run REGENIE
/sc-projects/sc-proj-computational-medicine/programs/regenie_v2.2.4.gz_x86_64_Centos7_mkl \
--step 1 \
--bed ${dir}/ukb22418_allChrs.pruned \
--extract ${dir}/ukb22418_allChrs.prune.in \
--keep ${input_dir}/sample_inclusion/EUR_panukbb_regenie_format_w20230425.id \
--phenoFile ${input_dir}/phenotypes.txt \
--covarFile ${input_dir}/covariates.txt \
--threads 110 \
--qt \
--bsize 1000 \
--lowmem \
--lowmem-prefix tmpdir/regenie_tmp_preds \
--out input/ukb_step1_pruned
