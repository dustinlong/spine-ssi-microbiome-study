
# If there is a need to update list of clinical taxonomic groups based on new taxa in sequencing table
update_clinical_taxonomy <- F

# Import current clinical taxonomic group mappings
species_list_classified_current <- read_excel(paste0(root_dir,"/taxonomy_curation/species_list_manualclassification_current.xlsx"))

if (update_clinical_taxonomy == T) {

# Count number of samples with at least 10 reads for each OTU
taxasums$count  <- rowSums(taxasums>10, na.rm=TRUE)

# Export list of species to categorize
taxonomy_review <- as.data.frame(rownames(taxasums))

taxonomy_review$count <- taxasums$count

colnames(taxonomy_review) <- c("taxonomy","total")
taxonomy_review <- subset(taxonomy_review, taxonomy_review$total>0)

# Remove column for previous total counts
species_list_classified_current <- species_list_classified_current[,which(!colnames(species_list_classified_current) %in% "total")]

# Merge imported curated list with latest taxonomy counts
taxonomy_review.merged<- merge(species_list_classified_current, taxonomy_review, by="taxonomy", all.x = TRUE, all.y = TRUE)

# sort by count
taxonomy_review.merged <- taxonomy_review.merged[order(-taxonomy_review.merged$total),]

# Rewrite classification table for updated review
write.xlsx(taxonomy_review.merged, file = "/taxonomy_curation/species_list_manualclassification_rawexport.xlsx",   row.names = F)
}



# Merge taxasums data with clinical taxonomic groups
taxamerge<-merge(species_list_classified_current, taxasums, by.x="taxonomy", by.y="row.names", all.x=F, all.y=F)




# Subsetting ----
# Taxonomy (row) subsetting
## Only taxa that have clinical classifaction as SSI-relevant organisms (potential pathogens)
taxamerge.select_taxa <-
  subset(taxamerge, !is.na(taxamerge$group_2))

# Summarize within taxonomic groups
## Aggregate sum of taxa counts within clinical groups for each sample that exists in the QIIME sequencing data
clinicaltaxa_summary <-as.data.frame(
  taxamerge.select_taxa %>%
  group_by(group_2) %>%
  summarise(across(c(
    subset(colnames(taxamerge.select_taxa),
           colnames(taxamerge.select_taxa) %in%
             as.character(metadata[,"SampleID"]))  ), ~ sum(.x, na.rm = TRUE)))
)


# Transpose
clinicaltaxa_summary.transposed <- as.data.frame(t(clinicaltaxa_summary), stringsAsFactors=F)

colnames(clinicaltaxa_summary.transposed) <- clinicaltaxa_summary$group_2

clinicaltaxa_summary.transposed <- clinicaltaxa_summary.transposed[-c(1),]

clinicaltaxa_summary.transposed$sample <-rownames(clinicaltaxa_summary.transposed)


# Calculate total fraction of reads for enteric organisms
clinicaltaxa_summary.transposed$enteric_total <- rowSums(
  sapply(
    clinicaltaxa_summary.transposed[ , c("Enterococcus sp.",
                                                        "Acinetobacter sp.",
                                                        "Stenotrophomonas sp.",
                                                        "Moraxella sp.",
                                                        "Escherichia-Shigella",
                                                        "Enterobacter sp.",
                                                        "Other Enterobacterales",
                                                        "Pseudomonas sp.",
                                                        "Bacteroides sp.",
                                                        "Fusobacterium sp.",
                                                        "Finegoldia sp.",
                                                        "Peptostreptococci",
                                                        "Lactobacilli",
                                                        "Prevotella sp."
                                                        )]
  ,  as.numeric)
)/rowSums(
  sapply(
    clinicaltaxa_summary.transposed[ , which(colnames(clinicaltaxa_summary.transposed) %in% c(
      "Staphylococcus sp.",
      "Streptococcus sp.",
      "Micrococci",
      
      "Cutibacteria",
      "Corynebacteria sp.",
      "Actinomyces",
      "Granulicatella",
      "Nocaria other",
      
      "Enterococcus sp.",
      
      "Acinetobacter sp.",
      "Stenotrophomonas sp.",
      "Moraxella sp.",
      "Escherichia-Shigella",
      "Enterobacter sp.",
      "Other Enterobacterales",
      "Pseudomonas sp.",
      
      "Bacteroides sp.",
      "Fusobacterium sp.",
      "Finegoldia sp.",
      "Peptostreptococci",
      "Lactobacilli",
      "Prevotella sp."
    ))]
    ,  as.numeric)
)


# Melt, keeping sample and enteric total rows, convert organisms columns to rows
clinicaltaxa_summary.transposed.melted <- melt(clinicaltaxa_summary.transposed, id=c("sample","enteric_total"))

clinicaltaxa_summary.transposed.melted$value <- as.numeric(clinicaltaxa_summary.transposed.melted$value)


