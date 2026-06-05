# Quantifying the multi-tissue response to glucose ingestion in humans by plasma proteomics

This repository contains the analysis code accompanying the paper published in Diabetologia (2026).

> Uluvar B, Williamson A, Kolnes KJ, Jeppesen PB, Kolnes AJ, Koprulu M, Zoodsma M, Beuchel C, Bambal Y, Reines M, Yasmeen S, Maj C, Schumacher J, O’Rahilly S, van Heel DA, Bartfeld S, Carrasco-Zanini J, Jensen J, Pietzner M, Langenberg C. 
> *Quantifying the multi-tissue response to glucose ingestion in humans by plasma proteomics.* Diabetologia (2026).

---

#### Note
The provided scripts are not designed to work out of the box, but illustrate the main analytical steps used to generate the results reported in the manuscript.
File paths have been replaced with path/to/file.

---

### Project structure

| Script | Description |
|--------|-------------|
| **scripts** | |
| [`scripts/01_proteomics_analysis.R`](scripts/01_proteomics_analysis.R) | Import QC'ed proteomics dataset; run linear models on protein trajectories (time-varying during each OGTT, differential across both OGTTs, sex-differential); sensitivity analyses; plot linear model results; plot individual protein trajectories; tissue enrichment |
| [`scripts/02_clinical_biomarkers.R`](scripts/02_clinical_biomarkers.R) | Import clinical biomarker dataset; run linear models; plot glucose and insulin trajectories; compute insulin indices and plot |
| [`scripts/03_stomach_scrna.R`](scripts/03_stomach_scrna.R) | Import gastrointestinal scRNAseq dataset; plot UMAP of cell types and ANXA10 expression; compute correlations of each gene expression with ANXA10; dotplots of correlated genes |
| [`scripts/04_disease_coloc.R`](scripts/04_disease_coloc.R) | Read ANXA10 summary statistics (generated with scripts in the in gwas directory); generate LD-matrix; join summary statistics for gastric cancer, gastric polyp, peptic ulcer, selected proteins from UKB-PPP; run coloc for binary and continuous outcomes separately |
| **gwas** | |
| [`scripts/gwas/01_prepare_input.R`](scripts/gwas/01_prepare_input.R) | Prepare phenotype and covariate files for REGENIE |
| [`scripts/gwas/02_step1_w_pruned_genotypes.sh`](scripts/gwas/02_step1_w_pruned_genotypes.sh) | Run REGENIE Step1 with LD-pruned high-quality genotyped SNPs |
| [`scripts/gwas/03_step2.sh`](scripts/gwas/03_step2.sh) | Run REGENIE Step2 for ANXA10 levels after some variant filtering |
| [`scripts/gwas/04_collate_output.sh`](scripts/gwas/04_collate_output.sh) | Collate GWAS output into one file and get regional lead variants |
| [`scripts/gwas/05_manhattan_plot.R`](scripts/gwas/05_manhattan_plot.R) | Create manhattan plot of ANXA10 levels using plot_manhattan function |
| **functions** | |
| [`scripts/functions/get_LD_matrix_ldstore.sh`](scripts/functions/get_LD_matrix_ldstore.sh) | Function to generate LD matrix based on UKB imputed BGEN files using LDstore2 |
| [`scripts/functions/mixed_effect_regression.R`](scripts/functions/mixed_effect_regression.R) | Function to run linear model with mixed effects |
| [`scripts/functions/mixed_effect_regression_interaction.R`](scripts/functions/mixed_effect_regression_interaction.R) | Function to run linear model with mixed effects for differential analysis |
| [`scripts/functions/plot_manhattan.R`](scripts/functions/plot_manhattan.R) | Function to make manhattan plot with ggplot |

---

### Software

| Software | Version | Purpose |
|----------|---------|---------|
| R | 4.3.x | All statistical analyses |
| REGENIE | 3.2.5 | Genome-wide and exome-wide association testing |
| LDStore | 2.0 | Generating LD matrix |
| coloc | 5.x | Bayesian colocalisation |

---

### Citation

Uluvar B, Williamson A, Kolnes KJ, Jeppesen PB, Kolnes AJ, Koprulu M, Zoodsma M, Beuchel C, Bambal Y, Reines M, Yasmeen S, Maj C, Schumacher J, O’Rahilly S, van Heel DA, Bartfeld S, Carrasco-Zanini J, Jensen J, Pietzner M, Langenberg C. Quantifying the multi-tissue response to glucose ingestion in humans by plasma proteomics. Diabetologia (2026)