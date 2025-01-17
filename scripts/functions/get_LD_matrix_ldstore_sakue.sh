#!/bin/sh

## script to create a LD-matrix
## Maik Pietzner 02/03/2022

## export location of files
# export dir=/sc-projects/sc-proj-computational-medicine/data/UK_biobank/genotypes/bgen_files_44448
export dir=/sc-resources/ukb/data/bulk/genetic/imputed/bgen_files_44448/

## get the chromosome
export chr=${1}
export lowpos=${2}
export uppos=${3}
export pheno=${4}

export tmpdir=tmpdir/sakue

## cd in the correct directory
# cd /sc-projects/sc-proj-computational-medicine/people/Maik/28_GWAS_EHR_embeddings/02_fine_mapping/01_SuSiE


echo "Chromosome ${chr} : Locus start ${lowpos} : Locus end ${uppos}"

if [ $chr -eq 23 ]; then

  ## create subset BGEN file (rsids of interest)
  /sc-projects/sc-proj-computational-medicine/programs/bgen/build/apps/bgenix \
  -g ${dir}/ukb22828_cX_b0_v3.bgen \
  -incl-rsids ${tmpdir}/snplist.${pheno}.${chr}.${lowpos}.${uppos}.lst > ${tmpdir}/filtered.${pheno}.${chr}.${lowpos}.${uppos}.bgen

else

  ## create subset BGEN file (rsids of interest)
  /sc-projects/sc-proj-computational-medicine/programs/bgen/build/apps/bgenix \
  -g ${dir}/ukb22828_c${chr}_b0_v3.bgen \
  -incl-rsids ${tmpdir}/snplist.${pheno}.${chr}.${lowpos}.${uppos}.lst > ${tmpdir}/filtered.${pheno}.${chr}.${lowpos}.${uppos}.bgen
  
fi

## create corresponding index file
/sc-projects/sc-proj-computational-medicine/programs/bgen/build/apps/bgenix \
-g ${tmpdir}/filtered.${pheno}.${chr}.${lowpos}.${uppos}.bgen \
-index


## run LD correlation
/sc-projects/sc-proj-computational-medicine/programs/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 \
--in-files ${tmpdir}/master.${pheno}.${chr}.${lowpos}.${uppos}.z \
--write-text \
--n-threads 30 \
--read-only-bgen

