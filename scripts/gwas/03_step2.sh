#!/bin/sh

## script to run 2nd step of REGENIE

#SBATCH --job-name=step2
#SBATCH --partition=compute
#SBATCH --account=your-slurm-account
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 120
#SBATCH --array=1-23%10
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm_logs/slurm-%x-%j.out
#SBATCH --error=slurm_logs/slurm-%x-%j.err

## ============================================================
## Configuration: set paths for your HPC environment
## ============================================================
export proj_dir="/path/to/your/gwas/directory"
export input="${proj_dir}/input"
export output="${proj_dir}"
export bgen_imp="/path/to/ukb/imputed/bgen_files"          # UKB imputed BGEN files
export snpstats="/path/to/ukb/imputed/variant_qc/EUR"      # UKB variant QC stats
export sample_qc="path/to/sample_inclusion_list"
export regenie="/path/to/regenie"                           # REGENIE binary
## ============================================================

## get the chromosome
echo "Job ID: $SLURM_ARRAY_TASK_ID"

## extract the chromosome to be run
export chr=$SLURM_ARRAY_TASK_ID

## change to X for 23
if [[ $chr -eq 23 ]]; then
  export chr="X"
fi

## change directory
cd ${output}

# Create variant list if it doesnt already exist
if [ ! -f "${output}/varqc_lists/variants_chr${chr}.list" ]; then
  echo "Generating variant list for Chromosome ${chr}..."
# Necessary because the structure of chromosomes 1-22 and X is different
  if [ $chr = "X" ]; then
      cat ${snpstats}/ukb_imp_EUR_chr${chr}_snpstat.out | \
      grep -v '#' | \
      awk -F '\t' -v chr=${chr} ' { if ( (NR != 1 ) &&
                                      ( $8 > 1e-15 ) &&
                                      ( $17 > 0.001 ) )
                                      print $2 }' > ${output}/varqc_lists/variants_chr${chr}.list

      echo "ChrX variant list created."
  else

    cat ${snpstats}/ukb_imp_EUR_chr${chr}_snpstat.out | \
    grep -v '#' | \
    awk -F '\t' -v chr=${chr} ' { if ( (NR != 1 ) &&
                                      ( $8 > 1e-15 ) &&
                                      ( $14 > 0.001 ) )
                                      print $2 }' > ${output}/varqc_lists/variants_chr${chr}.list
    echo "Chr${chr} variant list created."
  fi
else
    echo "Variant list for Chromosome ${chr} already exists. Skipping."
fi

## run REGENIE
${regenie} \
  --step 2 \
  --bgen ${bgen_imp}/ukb22828_c${chr}_b0_v3.bgen \
  --ref-first \
  --extract ${output}/varqc_lists/variants_chr${chr}.list \
  --sample ${bgen_imp}/ukb22828_c${chr}_b0_v3.sample \
  --keep ${sample_qc} \
  --phenoFile ${input}/phenotypes.txt \
  --covarFile ${input}/covariates.txt \
  --phenoCol anxa10_int \
  --threads 120 \
  --qt \
  --pred ${input}/ukb_step1_pruned_pred.list \
  --bsize 1000 \
  --out ${output}/gwas_output/gwas_chr${chr}
