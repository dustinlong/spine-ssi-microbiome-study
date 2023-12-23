# Install qiime2r
# if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
# devtools::install_github("jbisanz/qiime2R")

library(tidyverse)
library(qiime2R)
library(xlsx)
library(readxl)
library(reshape2)

metadata <- metadata_sheet

# QIIME import function renames sample-id to SampleID, standardize
colnames(metadata)[which(colnames(metadata)=="sample-id")]<-"SampleID"

# Derive detailed anatomic level, highest, lowest, median
for (sample.temp in metadata$SampleID){
  print(sample.temp)

  # Get values of all detailed anatomic level columns
  levels.all.temp<-as.data.frame(t(
    metadata[metadata$SampleID==sample.temp,
             grepl("anatomic_level_detailed___", colnames(metadata))
    ]
  ))

  # Subset only levels involved in operation
  levels.invovled.temp<-rownames(subset(levels.all.temp, levels.all.temp>0))
  
  if (length(levels.invovled.temp)>0){
  highest_level.temp <-levels.invovled.temp[1]
  lowest_level.temp <-levels.invovled.temp[length(levels.invovled.temp)]
  median_level.temp <- levels.invovled.temp[round(length(levels.invovled.temp)/2)]
  
  print(highest_level.temp)
  metadata[metadata$SampleID==sample.temp,"highest_level"]<-highest_level.temp
  metadata[metadata$SampleID==sample.temp,"lowest_level"]<-lowest_level.temp
  metadata[metadata$SampleID==sample.temp,"median_level"]<-median_level.temp
}
  
}

# Define order for setting factor levels
detailed_anatomic_levels <-c(
  "anatomic_level_detailed___occiput",
  "anatomic_level_detailed___c1",
  "anatomic_level_detailed___c2",
  "anatomic_level_detailed___c3",
  "anatomic_level_detailed___c4",
  "anatomic_level_detailed___c5",
  "anatomic_level_detailed___c6",
  "anatomic_level_detailed___c7",
  "anatomic_level_detailed___t1",
  "anatomic_level_detailed___t2",
  "anatomic_level_detailed___t3",
  "anatomic_level_detailed___t4",
  "anatomic_level_detailed___t5",
  "anatomic_level_detailed___t6",
  "anatomic_level_detailed___t7",
  "anatomic_level_detailed___t8",
  "anatomic_level_detailed___t9",
  "anatomic_level_detailed___t10",
  "anatomic_level_detailed___t11",
  "anatomic_level_detailed___t12",
  "anatomic_level_detailed___l1",
  "anatomic_level_detailed___l2",
  "anatomic_level_detailed___l3",
  "anatomic_level_detailed___l4",
  "anatomic_level_detailed___l5",
  "anatomic_level_detailed___sacrum",
  "anatomic_level_detailed___ileum",
  "anatomic_level_detailed___pelvis"
)


# Set factor levels
metadata$highest_level  <- factor(metadata$highest_level, levels = detailed_anatomic_levels)
metadata$median_level  <- factor(metadata$median_level, levels = detailed_anatomic_levels)
metadata$lowest_level  <- factor(metadata$lowest_level, levels = detailed_anatomic_levels)


# Set taxonomy and table files from QIIME artifacts
qiime_taxonomy_file <-paste0(root_dir,"/reports_intermediatefiles/taxonomy/taxonomy.qza")
qiime_table_file <- paste0(root_dir,"/reports_intermediatefiles/table/table-samplefiltered.qza")

# Read table file
SVs<-read_qza(qiime_table_file)$data

# Parse taxa
taxa_parsed<-read_qza(qiime_taxonomy_file)$data %>% parse_taxonomy()

# Summarize taxa at species level
taxasums<-summarize_taxa(SVs, taxa_parsed)$Species


# Create variable for aim1 skin sampling position
metadata[grepl("skin_pos_[0-9]+",metadata$Source),"aim1_position"]<-
  gsub("skin_pos_([0-9]+)", "\\1", 
       metadata[grepl("skin_pos_[0-9]+",metadata$Source),"Source"]
       , perl=TRUE)

metadata$aim1_position <- as.numeric(metadata$aim1_position)

# import PCA analysis
uwunifrac<-read_qza("./phylogeny/insertion-based/unweighted_unifrac_pcoa_results.qza")


