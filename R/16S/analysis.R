# Merge processed taxonomic table w metadata
clinicaltaxa_summary.transposed.melted.merged <-merge(clinicaltaxa_summary.transposed.melted, metadata, by.x = "sample", by.y = "SampleID")

# Set factor levels for cervical/thoracic/lumbar opertive level variable
clinicaltaxa_summary.transposed.melted.merged$`Operative Level` <- factor(clinicaltaxa_summary.transposed.melted.merged$`Operative Level`, levels = c("cervical", "thoracic", "lumbar"))

# Store lists of unique samples and taxa
sample.list <- unique(clinicaltaxa_summary.transposed.melted.merged$sample)
taxon.list <- unique(as.character(clinicaltaxa_summary.transposed.melted.merged$variable))

# Patient-level loop
sample.level.summary <- data.frame(sample=sample.list,
                                   stringsAsFactors=FALSE)
for (sample in sample.list){
  print(sample)
  
  # Record number of total reads for clinical taxa for sample
  sample.level.summary[sample.level.summary$sample==sample,"total_read_count.clinical"]<- sum(clinicaltaxa_summary.transposed.melted.merged[clinicaltaxa_summary.transposed.melted.merged$sample==sample,"value"])
  
  # Record total number of reads for all taxa for sample
  sample.level.summary[sample.level.summary$sample==sample,"total_read_count.all"]<- sum(taxasums[,sample])
  
  # Operative level category
  sample.level.summary[sample.level.summary$sample==sample,"Operative Level"]<- as.character(unique(clinicaltaxa_summary.transposed.melted.merged[clinicaltaxa_summary.transposed.melted.merged$sample==sample,"Operative Level"]))
  
  # Median Detailed Operative level
  sample.level.summary[sample.level.summary$sample==sample,"Median Operative Level"]<- as.numeric(unique(clinicaltaxa_summary.transposed.melted.merged[clinicaltaxa_summary.transposed.melted.merged$sample==sample,"median_level"]))
  
  # Highest Detailed Operative level
  sample.level.summary[sample.level.summary$sample==sample,"Highest Operative Level"]<- as.numeric(unique(clinicaltaxa_summary.transposed.melted.merged[clinicaltaxa_summary.transposed.melted.merged$sample==sample,"highest_level"]))
  
  # Lowest Detailed Operative level
  sample.level.summary[sample.level.summary$sample==sample,"Lowest Operative Level"]<- as.numeric(unique(clinicaltaxa_summary.transposed.melted.merged[clinicaltaxa_summary.transposed.melted.merged$sample==sample,"lowest_level"]))
  
  # Sex
  sample.level.summary[sample.level.summary$sample==sample,"Sex"]<- as.character(unique(clinicaltaxa_summary.transposed.melted.merged[clinicaltaxa_summary.transposed.melted.merged$sample==sample,"Sex"]))
  
  # Source
  sample.level.summary[sample.level.summary$sample==sample,"source"]<- as.character(unique(clinicaltaxa_summary.transposed.melted.merged[clinicaltaxa_summary.transposed.melted.merged$sample==sample,"Source"]))
  
  # Sample-taxon-level loop
  for (taxon in taxon.list){
    
    # Record relative abundance of that taxon for that sample
    sample.level.summary[sample.level.summary$sample==sample,taxon]<- clinicaltaxa_summary.transposed.melted.merged[clinicaltaxa_summary.transposed.melted.merged$sample==sample & clinicaltaxa_summary.transposed.melted.merged$variable==taxon,"value"]
    
  } # End sample-taxon-level loop
  
} # End patient-level loop


# Pairwise comparison of taxa
## Create result holders
  taxon.level.summary.continuous <- data.frame()
  
  taxon.level.summary.categorical <- data.frame(taxon=taxon.list,
                                                stringsAsFactors=FALSE)
  

# Subset to skin samples
sample.level.summary.subset<- subset(sample.level.summary, sample.level.summary$source== "Surgical Site")


# Continuous comparison by operative level (median_level) using Spearman rank sum correlation coefficient
for (taxon.temp in taxon.list){
  holder.temp<-data.frame("median_level"=as.numeric(sample.level.summary.subset[,"Median Operative Level"]),
                          "relative_abundance"=as.numeric((sample.level.summary.subset[,taxon.temp]/sample.level.summary.subset[,"total_read_count.all"]))
  )
  
  
  holder.temp <- subset(holder.temp, !is.na(holder.temp$median_level))
  
  print(taxon.temp)
  
  correlation.temp<-suppressWarnings(cor.test(holder.temp$median_level, holder.temp$relative_abundance, method="spearman"))
  
  taxon.level.summary.continuous<-rbind(taxon.level.summary.continuous,
                                        data.frame(
                                          taxon=taxon.temp,
                                          correlation=-correlation.temp$estimate,
                                          pvalue=correlation.temp$p.value,
                                          mean_relativeabundance=mean(((sample.level.summary.subset[,taxon.temp])/sample.level.summary.subset[,"total_read_count.all"])),
                                          color_1=species_list_classified_current[species_list_classified_current$group_2==taxon.temp & !is.na(species_list_classified_current$group_2),"color_1"][[1]][[1]],
                                          color_2=species_list_classified_current[species_list_classified_current$group_2==taxon.temp & !is.na(species_list_classified_current$group_2),"color_2"][[1]][[1]],
                                          species_group=species_list_classified_current[species_list_classified_current$group_2==taxon.temp & !is.na(species_list_classified_current$group_2),"micro_class"][[1]][[1]]
                                        ))
  
  
  
}


# Compare microbiome with SSI
# Load SSI correlation data from NHSN spinal fusion cases
load("/HMC IPC Spine SSI/HMC_spine_SSI_culture_gradient.Rda") 

# Average anatomic correlation values for SSI organism, grouped within their corresponding 16S taxonomic assignent groups
taxon.level.summary.continuous[taxon.level.summary.continuous$taxon=="Staphylococcus sp.","ssi_estimate"]<-
  mean(taxon.level.summary.continuous.ssi[taxon.level.summary.continuous.ssi$taxon %in% c("S. aureus","CoNS"), "correlation"])

taxon.level.summary.continuous[taxon.level.summary.continuous$taxon=="Enterococcus sp.","ssi_estimate"]<-
  mean(taxon.level.summary.continuous.ssi[taxon.level.summary.continuous.ssi$taxon %in% c("Enterococcus"), "correlation"])

taxon.level.summary.continuous[taxon.level.summary.continuous$taxon=="Cutibacteria","ssi_estimate"]<-
  mean(taxon.level.summary.continuous.ssi[taxon.level.summary.continuous.ssi$taxon %in% c("C. acnes"), "correlation"])

taxon.level.summary.continuous[taxon.level.summary.continuous$taxon=="Enterobacter sp.","ssi_estimate"]<-
  mean(taxon.level.summary.continuous.ssi[taxon.level.summary.continuous.ssi$taxon %in% c("Enterobacter"), "correlation"])

taxon.level.summary.continuous[taxon.level.summary.continuous$taxon=="Corynebacteria sp.","ssi_estimate"]<-
  mean(taxon.level.summary.continuous.ssi[taxon.level.summary.continuous.ssi$taxon %in% c("Other GPR"), "correlation"])

taxon.level.summary.continuous[taxon.level.summary.continuous$taxon=="Pseudomonas sp.","ssi_estimate"]<-
  mean(taxon.level.summary.continuous.ssi[taxon.level.summary.continuous.ssi$taxon %in% c("Pseudomonas"), "correlation"])

taxon.level.summary.continuous[taxon.level.summary.continuous$taxon=="Escherichia-Shigella","ssi_estimate"]<-
  mean(taxon.level.summary.continuous.ssi[taxon.level.summary.continuous.ssi$taxon %in% c("E. coli"), "correlation"])

taxon.level.summary.continuous[taxon.level.summary.continuous$taxon=="Streptococcus sp.","ssi_estimate"]<-
  mean(taxon.level.summary.continuous.ssi[taxon.level.summary.continuous.ssi$taxon %in% c("Streptococcus"), "correlation"])

taxon.level.summary.continuous[taxon.level.summary.continuous$taxon=="Other Enterobacterales","ssi_estimate"]<-
  mean(taxon.level.summary.continuous.ssi[taxon.level.summary.continuous.ssi$taxon %in% c("Proteus","Klebsiella","Morganella","Other GNR"), "correlation"])

taxon.level.summary.continuous<-rbind(taxon.level.summary.continuous, data.frame(taxon="Anaerobes",
                                                                                 correlation=mean(taxon.level.summary.continuous[taxon.level.summary.continuous$taxon %in% c("Lactobacilli","Bacteroides sp.","Prevotella sp.","Finegoldia sp.","Peptostreptococci"),"correlation"]),
                                                                                 pvalue=NA,
                                                                                 mean_relativeabundance=NA,
                                                                                 color_1="#800080",
                                                                                 color_2="#800080",
                                                                                 species_group="anaerobes",
                                                                                 ssi_estimate=taxon.level.summary.continuous.ssi[taxon.level.summary.continuous.ssi$taxon %in% c("Anaerobes"), "correlation"])
)

plot(taxon.level.summary.continuous$correlation, taxon.level.summary.continuous$ssi_estimate)
taxon.level.summary.continuous.subset <- subset(taxon.level.summary.continuous, !is.na(taxon.level.summary.continuous$ssi_estimate) & !is.na(taxon.level.summary.continuous$correlation))
taxon.level.summary.continuous.subset[taxon.level.summary.continuous.subset$taxon=="Enterococcus sp.","color_1"]<-"#86d9f0"
taxon.level.summary.continuous.subset<-subset(taxon.level.summary.continuous, !is.na(taxon.level.summary.continuous$correlation) & !is.na(taxon.level.summary.continuous$ssi_estimate))
cor.temp<-cor.test(taxon.level.summary.continuous.subset$correlation, taxon.level.summary.continuous.subset$ssi_estimate)


# Wilcoxon categorical testing
sample.level.summary.subset$wilcox.group <- 0

# Female, male (only surgical site)
sample.level.summary.subset[sample.level.summary.subset$source=="Surgical Site" & sample.level.summary.subset$Sex=="female" & !is.na(sample.level.summary.subset$Sex),"wilcox.group"]<-1
sample.level.summary.subset[sample.level.summary.subset$source=="Surgical Site" & sample.level.summary.subset$Sex=="male" & !is.na(sample.level.summary.subset$Sex),"wilcox.group"]<-2

for (taxon in taxon.list){ # Start taxon loop
  
  wilcox.test.temp<-wilcox.test(
    ((sample.level.summary.subset[sample.level.summary.subset$wilcox.group==1,taxon])/sample.level.summary.subset[sample.level.summary.subset$wilcox.group==1,"total_read_count.all"]),
    ((sample.level.summary.subset[sample.level.summary.subset$wilcox.group==2,taxon])/sample.level.summary.subset[sample.level.summary.subset$wilcox.group==2,"total_read_count.all"])
    ,paired = FALSE,
    conf.int = T)
  
  taxon.level.summary.categorical[taxon.level.summary.categorical$taxon==taxon,"CI_L"]<-wilcox.test.temp[["conf.int"]][1]
  
  
  taxon.level.summary.categorical[taxon.level.summary.categorical$taxon==taxon,"estimate"]<-as.numeric(
    wilcox.test.temp[["estimate"]][["difference in location"]]
  )
  
  taxon.level.summary.categorical[taxon.level.summary.categorical$taxon==taxon,"CI_H"]<-wilcox.test.temp[["conf.int"]][2]
  
  taxon.level.summary.categorical[taxon.level.summary.categorical$taxon==taxon,"pvalue"]<-as.numeric(
    wilcox.test.temp[["p.value"]][1]
  )
  
  
  
  taxon.level.summary.categorical[taxon.level.summary.categorical$taxon==taxon,"relativeabundance.median.group1"] <- 
    median(
      ((sample.level.summary.subset[sample.level.summary.subset$wilcox.group==1,taxon]) /
           sample.level.summary.subset[sample.level.summary.subset$wilcox.group==1,"total_read_count.all"]
        )
      )
  
  taxon.level.summary.categorical[taxon.level.summary.categorical$taxon==taxon,"relativeabundance.median.group2"] <- 
    median(
      ((sample.level.summary.subset[sample.level.summary.subset$wilcox.group==2,taxon]) /
         sample.level.summary.subset[sample.level.summary.subset$wilcox.group==2,"total_read_count.all"]
      )
    )
  
  
  taxon.level.summary.categorical[taxon.level.summary.categorical$taxon==taxon,"mean_relativeabundance"]<-mean(
    ((sample.level.summary.subset[,taxon])/sample.level.summary.subset[,"total_read_count.all"])
  )
  
  taxon.level.summary.categorical[taxon.level.summary.categorical$taxon==taxon,"color_1"] <- species_list_classified_current[species_list_classified_current$group_2==taxon & !is.na(species_list_classified_current$group_2),"color_1"][[1]][[1]]
  taxon.level.summary.categorical[taxon.level.summary.categorical$taxon==taxon,"color_2"] <- species_list_classified_current[species_list_classified_current$group_2==taxon & !is.na(species_list_classified_current$group_2),"color_2"][[1]][[1]]
} # End taxon loop


taxon.level.summary.categorical.export <- data.frame(taxon=taxon.level.summary.categorical$taxon,
                                                     estimate=taxon.level.summary.categorical$estimate,
                                                     CI95_L= taxon.level.summary.categorical$CI_L,
                                                     CI95_H= taxon.level.summary.categorical$CI_H,
                                                     p_value= taxon.level.summary.categorical$pvalue,
                                                     median_relative_abundance.female= taxon.level.summary.categorical$relativeabundance.median.group1,
                                                     median_relative_abundance.male= taxon.level.summary.categorical$relativeabundance.median.group2
                                                     )

# Standardize nomenclature
taxon.level.summary.categorical.export[taxon.level.summary.categorical.export$taxon=="Lactobacilli","taxon"]<- "Lactobacillus sp."
taxon.level.summary.categorical.export[taxon.level.summary.categorical.export$taxon=="Actinomyces","taxon"]<- "Actinomyces sp."
taxon.level.summary.categorical.export[taxon.level.summary.categorical.export$taxon=="Cutibacteria","taxon"]<- "Cutibacterium sp."

taxon.level.summary.categorical.export[taxon.level.summary.categorical.export$taxon=="Nocaria other","taxon"]<- "Nocaria sp."
taxon.level.summary.categorical.export[taxon.level.summary.categorical.export$taxon=="Granulicatella","taxon"]<- "Granulicatella sp."
taxon.level.summary.categorical.export[taxon.level.summary.categorical.export$taxon=="Escherichia-Shigella","taxon"]<- "Escherichia or Shigella sp."
taxon.level.summary.categorical.export[taxon.level.summary.categorical.export$taxon=="Micrococci","taxon"]<- "Micrococcus sp."


taxon.level.summary.categorical.export<-taxon.level.summary.categorical.export[order(taxon.level.summary.categorical.export$p_value),]


write.xlsx(taxon.level.summary.categorical.export, file=paste0(export_dir,"/tables/categorical_sex_differences.xlsx"))
