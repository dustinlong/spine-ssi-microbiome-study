# Set host machine parameters
host<-as.character(Sys.info()["nodename"])
source(Sys.getenv("LOCAL_PATH_SOURCE_FILE"))
root_dir <- # Set local project directory

# Set export directory
export_dir <- # Set export directory

# Set run name
runname <-"primary"

# Create workbook for holding figure data in tabular format
library(openxlsx)
tabular_data <- createWorkbook()

# 16S sequence data analysis---------
setwd(paste0(root_dir,"/R"))

## Generate metadata sheet from REDCap repository
file.edit("./16S/qiime_metadata_file_generation.R")

## After generating metadata files, run QIIME2 bash scripts
setwd(paste0(root_dir,"/QIIME2"))

### Build feature classifier
file.edit("silva_1388-99_naive-bayes_V3-V4.sh")
 
### Process and merge reads
file.edit("16S_process.sh")
 
### Check metadata and filter reads
file.edit("16S_merge.sh")
file.edit("filter_processed_reads.sh")

### Run classifier on filtered data and generate barplots
file.edit("taxonomic_classification_barplots.sh")

### Generate diversity/phylogeny metrics
file.edit("diversity_phylogeny.sh")

## Import QIIME2 artifacts into R
setwd(paste0(root_dir,"/R"))
file.edit("./16S/import_process_qiime_files.R")

## Curate and organize taxa by clinically relevant groupings
file.edit("./16S/curate_organize_taxonomic_groups.R")

## Initial analyses of 16S sequence data
file.edit("./16S/analysis.R")

## 16S Figure generation
file.edit("./16S/plotting.R")


# GenCap-Seq strain comparison -------- 

## Primary strain comparison analysis
file.edit("./GenCap-Seq/SSI_metagenome_pair_analyses.R")

## E. coli subgroup analysis
file.edit("./GenCap-Seq/Ecoli_interpatatient_comparision.R")

## GenCap-Seq enrichment efficacy
file.edit("./GenCap-Seq/enrichment_efficacy.R")


# AMR analysis ------
## Resistome correlation
file.edit("./AMR/resistome_correlation.R")


# Molecular epidemiology analysis -----
## Generate list of included isolates, analyze WGS distance by species and patient of origin
file.edit("./molecular_epidemiology/genomic_surviellance.R")


# Save workbook of data for figures in tabular format---
saveWorkbook(tabular_data, 
             paste0(export_dir,"/tables/tabular_fig_data/Data_File_S1.xlsx"), 
             overwrite = TRUE)
