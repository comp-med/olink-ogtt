#!/bin/sh

## script to run 2nd step of REGENIE

#SBATCH --job-name=step2
#SBATCH --partition=compute
#SBATCH --account=sc-users
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 120
#SBATCH --array=1-23%3
#SBATCH --mail-type=FAIL
#SBATCH --output=/sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/04_gwas/slurm_logs/slurm-%x-%j.out

## export location of files to be used
export dir=/sc-projects/sc-proj-computational-medicine/data/UK_biobank/genotypes/bgen_files_44448
export input_dir=/sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/04_gwas/input
export var_qc=/sc-projects/sc-proj-computational-medicine/data/UK_biobank/genotypes/variant_qc

## get the chromosome
echo "Job ID: $SLURM_ARRAY_TASK_ID"

## extract the chromosome to be run
export chr=$SLURM_ARRAY_TASK_ID

## change to X for 23
if [[ $chr -eq 23 ]]; then
  export chr="X"  
fi

## change directory
cd /sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/04_gwas

## create variant inclusion list for the respective chromosome - already done in prev step for chr13
#cat ${var_qc}/output/ukb_imp_chr${chr}_qced.txt | awk -v chr=${chr} '{if(NR != 1 && $3 == chr) print $2}' - > tmpdir/tmp.ex.chr${chr}.list
  
## run REGENIE
/sc-projects/sc-proj-computational-medicine/programs/regenie_v2.2.4.gz_x86_64_Centos7_mkl \
  --step 2 \
  --bgen ${dir}/ukb22828_c${chr}_b0_v3.bgen \
  --ref-first \
  --extract tmpdir/tmp.ex.chr${chr}.list \
  --sample ${dir}/ukb22828_c${chr}_b0_v3.sample \
  --keep ${input_dir}/sample_inclusion/EUR_panukbb_regenie_format_w20230425.id \
  --phenoFile ${input_dir}/phenotypes.txt \
  --covarFile ${input_dir}/covariates.txt \
  --threads 120 \
  --qt \
  --pred ${input_dir}/ukb_step1_pruned_pred.list \
  --bsize 1000 \
  --out output/ukb_step2_chr${chr}