# Functional enrichment analysis using the R package topGO

### Scripts:
**lfmm2proteinID.R**: create input files for topGO using calibrated p-values from LFMM as input  
**top_go.R**: functional enrichment analysis and plotting  
**haplotype_frequencies.R**: compare functionally significant haplotypes from each genetic cluster and plot as pie charts

## 
### Associated files:
**qvals4go.RData.gz**: calibrated p-values from LFMM  
**sizes.genome**: scaffold sizes  
**Nerluet1_FilteredModels1_2024-03-15.gff3.gz**: annotation file containing GO terms, positions, and protein IDs  
GO.DF* : files created by **lfmm2proteinID.R** and used as input in **top_go.R**  
