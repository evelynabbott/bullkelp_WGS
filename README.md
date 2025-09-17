# Details of data analysis pipeline for Puget Sound Bull kelp genomics project (2025).
Evelyn Abbott, Filipe Alberto

## Main directory:
**data_processing.sh**: creating VCF file, quality filtering, identity by state, and admixture analysis

**linkage_disequilibrium.sh**: estimating LD within each genetic cluster using plink  

The outputs from these scripts are plotted using **pop_structure.R**

## Subdirectories:  
**RDA_forest**: calculating importances of environmental variables and plotting adaptive neighborhoods 

**Environmental_data**: sourcing environmental data, imputing rasters, extracting data from the Salish Sea Model (Yang and Khangaonkar, 2010), and the final tables used in downstream analyses. 

**Functional_analyses**: identifying enriched gene ontology (GO) terms from environmentally-associated SNPs  

**Snp_environment_association**: calculating genotype x environment correlations using latent factor mixed models (LFMM)


