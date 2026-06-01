#!/bin/sh

## script to create a LD-matrix using LDstore2

## ============================================================
## Configuration: set paths for your HPC environment
## ============================================================
export bgen_dir="/path/to/ukb/imputed/bgen_files"   # UKB imputed BGEN files
export bgenix="/path/to/bgenix"                      # bgenix binary
export ldstore="/path/to/ldstore_v2.0/ldstore"       # LDstore2 binary
## ============================================================

## get the chromosome
export chr=${1}
export lowpos=${2}
export uppos=${3}
export pheno=${4}

export tmpdir=tmpdir/sakue

echo "Chromosome ${chr} : Locus start ${lowpos} : Locus end ${uppos}"

if [ $chr -eq 23 ]; then

  ## create subset BGEN file (rsids of interest)
  ${bgenix} \
  -g ${bgen_dir}/ukb22828_cX_b0_v3.bgen \
  -incl-rsids ${tmpdir}/snplist.${pheno}.${chr}.${lowpos}.${uppos}.lst > ${tmpdir}/filtered.${pheno}.${chr}.${lowpos}.${uppos}.bgen

else

  ## create subset BGEN file (rsids of interest)
  ${bgenix} \
  -g ${bgen_dir}/ukb22828_c${chr}_b0_v3.bgen \
  -incl-rsids ${tmpdir}/snplist.${pheno}.${chr}.${lowpos}.${uppos}.lst > ${tmpdir}/filtered.${pheno}.${chr}.${lowpos}.${uppos}.bgen

fi

## create corresponding index file
${bgenix} \
-g ${tmpdir}/filtered.${pheno}.${chr}.${lowpos}.${uppos}.bgen \
-index


## run LD correlation
${ldstore} \
--in-files ${tmpdir}/master.${pheno}.${chr}.${lowpos}.${uppos}.z \
--write-text \
--n-threads 30 \
--read-only-bgen
