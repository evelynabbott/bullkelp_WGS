# Sourcing environmental data and creating files for downstream analyses
### Most environmental data are sourced from NOAA coastwatch links provided in get_rasters.R.  
### Salinity raster is sourced from the Salish Sea Model  (SSM, [Yang and Khangaonkar, 2010](https://link.springer.com/article/10.1007/s10236-010-0348-5)).
## 

### Scripts:
get_rasters.R: links for downloading environmental data  
get_salinity_nodes.R: create list of nodes to extract from SSM  
salinity_ssm.py: extract and average salinity data across one year  
raster_imputation.R: impute missing data and standardize raster cell size  
Maps_with_rasters.R plots maps displaying environmental data
## 
### Associated files:
rasters4maps.Rdata.gz: dataframe containing all imputed/standardized rasters  
buoy_data.tsv: environmental data for each site (buoy)  
sample_coordinates.tsv: latitude and longitude of all samples
