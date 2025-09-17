# SnpEff: variants which impact gene function
This repository details data analysis following the identification of variants using SnpEff ([Cingolani, 2012](https://www.tandfonline.com/doi/full/10.4161/fly.19695)). Information on upstream processes can be found [here](https://pcingola.github.io/SnpEff/snpeff/introduction/).  

## Scripts:  
**impact_by_cluster.R**: initial exploration + plotting variant frequencies by cluster  
**snpeff_analyses.R**: frequencies of variants by population and cluster; additional statistics  
**variant_effects.R**: create tables listing the specific effects each variant has on gene function  
**variant2protein_scan.R**: match variants with protein IDs (based on location)  

## Associated files:
**rawfiles.tar.gz**: result files returned by SnpEff  
**snpeff_results.Rdata**: SnpEff results as a list of dataframes in R  
**Nerluet1_FilteredModels1_2024-03-15.gff3.gz**: annotation file  
**env_65.Rdata**: environmental and sample metadata  
**impact.prot.Rdata**: dataframe of all SnpEff results along with protein IDs

