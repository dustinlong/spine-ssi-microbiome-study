library(ggrepel)
library(openxlsx)

# Set run parameters
## Repeat preop-DOS comparison
run_preop_dos_analysis <-T
# Overwrite figure export files
write_new_figs <-F


# Plot microbiome-SSI anatomic correlation by taxonomic group----
colors <- as.character(taxon.level.summary.continuous.subset$color_1)
names(colors) <- as.character(taxon.level.summary.continuous.subset$color_1)

# Sort clinical classes
taxon.level.summary.continuous.subset$species_group <- factor(taxon.level.summary.continuous.subset$species_group, levels = c( "gram positive cocci",
                                                                                                                               "gram positive  rods",
                                                                                                                               "enteric gram positive cocci",
                                                                                                                               "gram negative rods",
                                                                                                                               "anaerobes"))

SSI_microbiome_correlation<-
  ggplot(data = taxon.level.summary.continuous.subset, aes(x = -ssi_estimate, y = correlation))+ 
  geom_point(size=5,aes(color=species_group), show.legend = TRUE)+
  geom_smooth(method = "glm",alpha=0.1, size=0, level=.95, fullrange=T, linetype="blank")+
  stat_smooth(method="lm",geom="line", alpha=0.1, size=0.5, linetype="dashed", span=10^100) +
  expand_limits(x = c(.8,-0.8), y=c(0.4,-0.4))+
  # Primary analysis
  scale_color_manual(values=c("lightblue","blue","grey","red","purple"),
                     labels=c( "gram positive cocci",
                               "gram positive  rods",
                               "enteric gram positive cocci",
                               "gram negative rods",
                               "anaerobes")
                     )+
  xlab("Surgical Site Infection ~ Spinal Anatomic Level (R^2 Correlation)")+
  ylab("Preoperative Skin Microbiome Overlying Incision ~ Spinal Anatomic Level (R^2 Correlation)")+
  geom_text_repel(aes(label = taxon),
                  box.padding   = 0.5, 
                  point.padding = .02,
                  segment.color = 'grey50',
                  size=3, 
                  segment.alpha = 0.5,
                  min.segment.length = 2,
                  max.overlaps = 100
  ) +
  theme_classic()+
  labs(color = "Clinical Microbiological Group")+
  guides(colour = guide_legend(override.aes = list(size=5)))



SSI_microbiome_correlation


if (write_new_figs==T){
ggsave(
  paste0(export_dir,"/figures/ssimicrobiomecorrelation_", runname, ".pdf"),
  plot=SSI_microbiome_correlation,
  scale = 1,
  width = 30,
  height = 20,
  units = "cm",
  dpi = 300
)
}

# Stage tabular data for export
taxon.level.summary.continuous.subset.export <- data.frame(Taxon=taxon.level.summary.continuous.subset$taxon,
                                                           Microbiome_Anatomy_Correlation=taxon.level.summary.continuous.subset$correlation,
                                                           SSI_Anatomy_Correlation=taxon.level.summary.continuous.subset$ssi_estimate,
                                                           Clinical_group=taxon.level.summary.continuous.subset$species_group
  
)


addWorksheet(tabular_data, sheetName = "Fig 4")
writeDataTable(tabular_data, sheet = "Fig 4", x = taxon.level.summary.continuous.subset.export,
               colNames = T, rowNames = F)


# Patient-level microbiome-anatomy plotting---
# Sort clinical taxa
clinicaltaxa_summary.transposed.melted.merged.sorted <- clinicaltaxa_summary.transposed.melted.merged[order(
  clinicaltaxa_summary.transposed.melted.merged$Source
  , clinicaltaxa_summary.transposed.melted.merged$lowest_level
  , clinicaltaxa_summary.transposed.melted.merged$highest_level
  , clinicaltaxa_summary.transposed.melted.merged$enteric_total
  ),]

# Order taxon factor levels to facilitate fig legend ordering
clinicaltaxa_summary.transposed.melted.merged.sorted$variable <- factor(clinicaltaxa_summary.transposed.melted.merged.sorted$variable, levels=c(
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
    "Pseudomonas sp.",
  "Other Enterobacterales",
  
  "Fusobacterium sp.",
  "Peptostreptococci",
  "Lactobacilli",
  "Finegoldia sp.",
  "Prevotella sp.",
  "Bacteroides sp."
))

# Set colors manually
relative_abundance_palette<-c(
  "#caf0f8", # Staphylococcus sp.
  "#90e0ef", # Streptococcus sp.
  "#48cae4", # Micrococci
  
  "#00b4d8", # Cutibacteria
  "#0096c7", # Corynebacteria sp.
  "#0077b6", # Actinomyces
  "#3c79cf", # Granulicatella
  "#023e8a", # Nocaria other #
  
  "#FFE591", # Enterococcus sp.
  
  "#FFC500", # Acinetobacter sp.
  "#FFAB0A", # Stenotrophomonas sp.
  "#FC7921", # Moraxella sp.
  "#FFA3A4", # Escherichia-Shigella
  "#FF857A", # Enterobacter sp.
  "#FF5963", # Pseudomonas sp.
  "#FF2700", # Other Enterobacterales
  
  "#DC97FF", # Fusobacterium sp.
  "#BD68EE", # Peptostreptococci
  "#AB51E3", # Lactobacilli
  "#8B2FC9", # Finegoldia sp.
  "#733BCE", # Prevotella sp.
  "#5A108F" # Bacteroides sp.
)

# Subset to skin samples
clinicaltaxa_summary.transposed.melted.merged.sorted.subset <- subset(clinicaltaxa_summary.transposed.melted.merged.sorted, clinicaltaxa_summary.transposed.melted.merged.sorted$Source=="Surgical Site")


# Loop through anatomic levels
level_summary_data <- data.frame(taxon=as.character())

for(level.temp in detailed_anatomic_levels){ # Anatomic level loop begin
  
  # Count total number of reads at that level
  total_level_read.temp <-  sum(
    clinicaltaxa_summary.transposed.melted.merged.sorted.subset[
      clinicaltaxa_summary.transposed.melted.merged.sorted.subset[,
                                                                  which(colnames(clinicaltaxa_summary.transposed.melted.merged.sorted.subset)==level.temp)]==1,"value"]
    
    , na.rm = T)
  
  
  # Get median abundance values for all taxa at that level
  vals.temp<-(
    clinicaltaxa_summary.transposed.melted.merged.sorted.subset[clinicaltaxa_summary.transposed.melted.merged.sorted.subset[,
                                                                                                                            which(colnames(clinicaltaxa_summary.transposed.melted.merged.sorted.subset)==level.temp)]==1,] %>%
    group_by(variable) %>%
    summarize(median_abundance=median(value/total_level_read.temp))
  )
  
  # Set the column name to the anatomic level
  colnames(vals.temp)[2]<-level.temp
  
  # Merge into larger anatomic level summary table
  level_summary_data <- merge(level_summary_data, vals.temp, by.x="taxon", by.y = "variable", all = T)
  

  }# Antomic level loop end


# Sort according to operative level group and %enteric
clinicaltaxa_summary.transposed.melted.merged.sorted.subset$sample<-
  factor(clinicaltaxa_summary.transposed.melted.merged.sorted.subset$sample,
         levels = unique(clinicaltaxa_summary.transposed.melted.merged.sorted.subset[order(
                                                                                           clinicaltaxa_summary.transposed.melted.merged.sorted.subset$`Operative Level`,
                                                                                           clinicaltaxa_summary.transposed.melted.merged.sorted.subset$enteric_total
                                                                                           ),"sample"]
         ))



# Plot relative abundance (clinical taxa), patient level, by anatomic group
library(extrafont)
extrafont::loadfonts(quiet = TRUE)


relativeabundance_anatomic_patlevel <-ggplot(clinicaltaxa_summary.transposed.melted.merged.sorted.subset, aes(fill=variable, 
                                                                        x=sample,
                                                                        y=value)) + 
  geom_bar(position="fill", stat="identity") +
  theme_minimal()+
  scale_fill_manual(values=relative_abundance_palette)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_cartesian(ylim=c(0,1), expand = 0, )+
  scale_y_continuous(labels=scales::percent) +
  labs(fill = "Organism")+
  theme(legend.title.align=0.5)+
  #ylab("")+
  theme(axis.title.y=element_blank()
        ,axis.title.x=element_blank()
        ,axis.text.x=element_blank()
        ,axis.ticks.x=element_blank()
        )+
  theme(axis.text.y = element_text(size=12 
                                   )) +
  guides(fill=guide_legend(ncol=1))

relativeabundance_anatomic_patlevel


if (write_new_figs==T){
ggsave(
  paste0(export_dir,"/figures/abundance_by_operativelevelcategory_", runname, ".pdf"),
  plot = relativeabundance_anatomic_patlevel,
  scale = 1,
  width = 40,
  height = 20,
  units = "cm",
  dpi = 300
)
}


# Stage tabular data for export
clinicaltaxa_summary.transposed.melted.merged.sorted.subset.export <- data.frame(
  Sample=clinicaltaxa_summary.transposed.melted.merged.sorted.subset$sample,
  Taxon=clinicaltaxa_summary.transposed.melted.merged.sorted.subset$variable,
  Read_Count=clinicaltaxa_summary.transposed.melted.merged.sorted.subset$value,
  Operative_Level_Group=clinicaltaxa_summary.transposed.melted.merged.sorted.subset$`Operative Level`
  
)


addWorksheet(tabular_data, sheetName = "Fig 2A")
writeDataTable(tabular_data, sheet = "Fig 2A", x = clinicaltaxa_summary.transposed.melted.merged.sorted.subset.export,
               colNames = T, rowNames = F)



# Run a similar loop as above, but as a moving weighted average leveraging data from two adjacent levels to approximate average surgical span of operative levels in cohort
level_summary_data.weighted <- data.frame()
level_summary_data.weighted.iqr <- data.frame()

for(levelno.temp in 1:length(detailed_anatomic_levels)){ # Start loop for moving average by median operative level

  levelno_1above.temp <- levelno.temp-1
  levelno_2above.temp <- levelno.temp-2
  
  levelno_1below.temp <- levelno.temp+1
  levelno_2below.temp <- levelno.temp+2
  
  for (taxon.temp in subset(level_summary_data$taxon, !is.na(level_summary_data$taxon))){

    # Get moving median value with 50% weight reduction per level removed up to 2 levels removed
    level_summary_data.weighted[taxon.temp,detailed_anatomic_levels[levelno.temp]]<- median(c(
      level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno.temp]]
      ,level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno.temp]]
      ,level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno.temp]]
      ,level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno.temp]]
      , if(levelno_1above.temp >=1 & levelno_1above.temp <= length(detailed_anatomic_levels)) {level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno_1above.temp]]} else {NA}
      , if(levelno_1above.temp >=1 & levelno_1above.temp <= length(detailed_anatomic_levels)) {level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno_1above.temp]]} else {NA}
      , if(levelno_1below.temp >=1 & levelno_1below.temp <= length(detailed_anatomic_levels)) {level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno_1below.temp]]} else {NA}
      , if(levelno_1below.temp >=1 & levelno_1below.temp <= length(detailed_anatomic_levels)) {level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno_1below.temp]]} else {NA}
      , if(levelno_2above.temp >=1 & levelno_2above.temp <= length(detailed_anatomic_levels)) {level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno_2above.temp]]} else {NA}
      , if(levelno_2below.temp >=1 & levelno_2below.temp <= length(detailed_anatomic_levels)) {level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno_2below.temp]]} else {NA}
    ), na.rm = T)
    
    
    
    # Calculate IQR in same fashion
    level_summary_data.weighted.iqr[taxon.temp,detailed_anatomic_levels[levelno.temp]]<- IQR(c(
      level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno.temp]]
      ,level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno.temp]]
      ,level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno.temp]]
      ,level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno.temp]]
      , if(levelno_1above.temp >=1 & levelno_1above.temp <= length(detailed_anatomic_levels)) {level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno_1above.temp]]} else {NA}
      , if(levelno_1above.temp >=1 & levelno_1above.temp <= length(detailed_anatomic_levels)) {level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno_1above.temp]]} else {NA}
      , if(levelno_1below.temp >=1 & levelno_1below.temp <= length(detailed_anatomic_levels)) {level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno_1below.temp]]} else {NA}
      , if(levelno_1below.temp >=1 & levelno_1below.temp <= length(detailed_anatomic_levels)) {level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno_1below.temp]]} else {NA}
      , if(levelno_2above.temp >=1 & levelno_2above.temp <= length(detailed_anatomic_levels)) {level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno_2above.temp]]} else {NA}
      , if(levelno_2below.temp >=1 & levelno_2below.temp <= length(detailed_anatomic_levels)) {level_summary_data[level_summary_data$taxon==taxon.temp, detailed_anatomic_levels[levelno_2below.temp]]} else {NA}
    ), na.rm = T)
    
  }
  
} # End for level, moving average

level_summary_data.weighted$taxon <- rownames(level_summary_data.weighted)
level_summary_data.weighted.melted <- melt(level_summary_data.weighted, id=c("taxon"), variable.name="anatomic_level", value.name="median")


level_summary_data.weighted.iqr$taxon <- rownames(level_summary_data.weighted.iqr)
level_summary_data.weighted.iqr.melted <- melt(level_summary_data.weighted.iqr, id=c("taxon"), variable.name="anatomic_level", value.name="IQR")



# Set factor levels
level_summary_data.weighted.melted$taxon <- factor(level_summary_data.weighted.melted$taxon, levels=c(
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
  "Pseudomonas sp.",
  "Other Enterobacterales",
  
  "Fusobacterium sp.",
  "Peptostreptococci",
  "Lactobacilli",
  "Finegoldia sp.",
  "Prevotella sp.",
  "Bacteroides sp."
))



# Plot relative abundance (clinical taxa), by operative level, median weighted
relativeabundance_anatomic_continuous <- 
  ggplot(level_summary_data.weighted.melted, aes(fill=taxon, x=anatomic_level, y=median)) + 
  geom_bar(position="fill", stat="identity") +
  theme_minimal()+
  scale_fill_manual(values=relative_abundance_palette)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_cartesian(ylim=c(0,1), expand = 0, )+
  labs(fill = "Organism")+
  theme(legend.title.align=0.5)+
  ylab("Relative Abundance")+
  xlab("Sample")+
 scale_x_discrete(limits = unique(level_summary_data.weighted.melted$anatomic_level),
                   labels=unique(
                     gsub("anatomic_level_detailed___" , "", as.character(level_summary_data.weighted.melted$anatomic_level))
                                 ))+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none"
  )+
  theme(axis.text.y = element_text(size=12))+
  scale_y_continuous(labels=scales::percent) 

relativeabundance_anatomic_continuous
  
if (write_new_figs==T){
ggsave(
  paste0(export_dir,"/figures/abundance_by_operativelevelcontinuousweighted_", runname, ".pdf"),
  plot = relativeabundance_anatomic_continuous,
  scale = 1,
  width = 30,
  height = 15,
  units = "cm",
  dpi = 300
)
}

# Merge in IQR values, stage for tabular data export
abundance_by_operativelevelcontinuousweighted.valuetable <- left_join(level_summary_data.weighted.melted, level_summary_data.weighted.iqr.melted)


addWorksheet(tabular_data, sheetName = "Fig 2B")
writeDataTable(tabular_data, sheet = "Fig 2B", x = abundance_by_operativelevelcontinuousweighted.valuetable,
               colNames = T, rowNames = F)





# Comparison of preop vs DOS samples from Aim 1 patients ---
if (run_preop_dos_analysis == T){
  
# Subset paired samples for longitudinally sampled patients (spine clinic swab at site corresponding to DOS skin site and matched day of surgery swab)
clinicaltaxa_summary.transposed.melted.merged.sorted.preopdoscomp<- subset(clinicaltaxa_summary.transposed.melted.merged.sorted,
       clinicaltaxa_summary.transposed.melted.merged.sorted$sample %in% c("2-16-dos-skin", "1-14-dos-skin","3-14-dos-skin","5-16-dos-skin")
       | 	(clinicaltaxa_summary.transposed.melted.merged.sorted$sample %in% c("1-10-skin-pos-9","2-8-skin-pos-7","3-9-skin-pos-8","5-9-skin-pos-8"))

       )


# For each patient/taxon combination, calculate the relative abundance for the taxon
for(row.temp in (rownames(clinicaltaxa_summary.transposed.melted.merged.sorted.preopdoscomp))){
  sample.temp <- unique(clinicaltaxa_summary.transposed.melted.merged.sorted.preopdoscomp[row.temp,"sample"])
  clinicaltaxa_summary.transposed.melted.merged.sorted.preopdoscomp[row.temp,"relative_abundance.all"]<-
    clinicaltaxa_summary.transposed.melted.merged.sorted.preopdoscomp[row.temp,"value"]/sample.level.summary[sample.level.summary$sample==sample.temp,"total_read_count.all"]
  
} 


# Classify clinic vs day-of-surgery samples
clinicaltaxa_summary.transposed.melted.merged.sorted.preopdoscomp$timepoint_label<-ifelse(grepl("-dos-", clinicaltaxa_summary.transposed.melted.merged.sorted.preopdoscomp$sample), "Day of Surgery","Clinic")


# Preop clinic vs DOS patient-level plot
preop_dos_plot<-ggplot(
  clinicaltaxa_summary.transposed.melted.merged.sorted.preopdoscomp, 
  aes(
    x      = factor(variable), 
    group  = factor(timepoint_label),
    fill   = factor(timepoint_label),
    y=relative_abundance.all
  )) +
  geom_bar( stat = "identity", position = "dodge") +
  facet_wrap(~record_id,nrow=4) +
  ggtitle('Patient') +
  xlab('Sample') +
  ylab('Count') +
  labs(fill = 'Sampling Encounter')+
  scale_fill_manual(values = c("#99d98c", "#1a759f"), labels=c("Preoperative Clinic", "Day of Surgery"))+
 scale_y_sqrt(labels = scales::percent, breaks=c(0.01,0.1,0.25,0.5,0.75,1))+
  theme_minimal()+
  theme(axis.text.x = element_text(size= 12, angle = -45, vjust = 0.5, hjust=0), axis.text.y = element_text(size=12), axis.title.y = element_blank(), axis.title.x = element_blank()
        , panel.grid.major.y = element_blank())
  
preop_dos_plot

if (write_new_figs==T){
ggsave(
  paste0(export_dir,"/figures/preop_DOS_16S_comparision", runname, ".pdf"),
  plot = preop_dos_plot,
  scale = 1,
  width = 20,
  height = 20,
  units = "cm",
  dpi = 300
)
}

# Prep and export tabular data
clinicaltaxa_summary.transposed.melted.merged.sorted.preopdoscomp.export <- data.frame(
  subject_id=clinicaltaxa_summary.transposed.melted.merged.sorted.preopdoscomp$subject_id,
  taxon=clinicaltaxa_summary.transposed.melted.merged.sorted.preopdoscomp$variable,
  read_count=clinicaltaxa_summary.transposed.melted.merged.sorted.preopdoscomp$value,
  timepoint=clinicaltaxa_summary.transposed.melted.merged.sorted.preopdoscomp$timepoint_label
)

addWorksheet(tabular_data, sheetName = "Fig S2")
writeDataTable(tabular_data, sheet = "Fig S2", 
               x = clinicaltaxa_summary.transposed.melted.merged.sorted.preopdoscomp.export,
               colNames = T, rowNames = F)


}





# For each posterior instrumented fusion case sampled at multiple skin positions in clinic, generate individual plot
aim1.tabulardata<-data.frame()

for (aim1subject.temp in c(1,2,3,5)){

# Set sequential labels for figure
subject.figlabel <- ifelse(aim1subject.temp==5, 4, aim1subject.temp)  
  
# Subset to aim 1 skin samples
aim1.subset<-subset(clinicaltaxa_summary.transposed.melted.merged.sorted, clinicaltaxa_summary.transposed.melted.merged.sorted$aim1==1)

aim1.subset<-subset(aim1.subset, grepl("skin_pos_",aim1.subset$Source))

aim1.subset<-subset(aim1.subset, aim1.subset$subject_id==aim1subject.temp)

# Set factor levels by vertebral anatomic level
aim1_postion_levels<-c("skin_pos_1",
                       "skin_pos_2",
                       "skin_pos_3",
                       "skin_pos_4",
                       "skin_pos_5",
                       "skin_pos_6",
                       "skin_pos_7",
                       "skin_pos_8",
                       "skin_pos_9",
                       "skin_pos_10",
                       "skin_pos_11")

aim1.subset$Source<-factor(aim1.subset$Source, levels=aim1_postion_levels)


# Plot relative abundance (clinical taxa), patient level
aim1_plot.temp<-ggplot(aim1.subset, aes(fill=variable, x=Source, y=value)) +
  geom_bar(position="fill", stat="identity") +
  theme_minimal()+
  scale_fill_manual(values=relative_abundance_palette)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_cartesian(ylim=c(0,1), expand = 0, )+
  scale_y_continuous(labels=scales::percent) +
  labs(fill = "Organism")+
  theme(legend.title.align=0.5)+
  ylab("Relative Abundance")+
  xlab("Sample")+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none"
  )+
  theme(axis.text.y = element_text(size=15)) 

aim1_plot.temp

if (write_new_figs==T){
ggsave(
  paste0(export_dir,"/figures/aim1_subject",aim1subject.temp,"_",runname,".pdf"),
  plot = aim1_plot.temp,
  scale = 1,
  width = 20,
  height = 7,
  units = "cm",
  dpi = 300
)
}

# Append data for all aim 1 subjets in loop
aim1.tabulardata <-rbind(aim1.tabulardata,
                         data.frame(
                           subject_label=subject.figlabel,
                           anatomic_level=aim1.subset$Source,
                           taxon=aim1.subset$variable,
                           read_count=aim1.subset$value
                         ))


} # End for aim 1 subject

# Write tabular data
addWorksheet(tabular_data, sheetName = "Fig 2C")

writeDataTable(tabular_data, sheet = "Fig 2C", 
               x = aim1.tabulardata,
               colNames = T, rowNames = F)




# Anatomic niche plot

# Subset to 16S-derived data by removing null group classifications from SSI-16S correlation steps
taxon.level.summary.continuous <- subset(taxon.level.summary.continuous, !is.na(taxon.level.summary.continuous$pvalue))

# Add labels for common primary anatomic reservoir
taxon.level.summary.continuous[taxon.level.summary.continuous$color_2=="#86d9f0","anatomic_reservoir"]<-"Skin and nasal"
taxon.level.summary.continuous[taxon.level.summary.continuous$color_2=="#3cb5d6","anatomic_reservoir"]<-"Aerodigestive"
taxon.level.summary.continuous[taxon.level.summary.continuous$color_2=="#7400b8","anatomic_reservoir"]<-"Mixed commensals"
taxon.level.summary.continuous[taxon.level.summary.continuous$color_2=="#e63946","anatomic_reservoir"]<-"Enteric"


# Sort by clinical classes for fig legend
taxon.level.summary.continuous$anatomic_reservoir <- factor(taxon.level.summary.continuous$anatomic_reservoir, levels = c( "Skin and nasal",
                                                                                                                  "Aerodigestive",
                                                                                                                  "Mixed commensals",
                                                                                                                 "Enteric"
                                                                                                                 ))

# Generate anatomic niche plot
anatomic_niche_plot <-
ggplot(data = taxon.level.summary.continuous, aes(x = -correlation, y = -sign(correlation)*log10(pvalue)))+ 
  geom_vline(xintercept = 0, color="grey10", size=.1, linetype="dashed")+
  geom_point(aes(size=mean_relativeabundance,color=anatomic_reservoir))+
  scale_color_manual(values=c("#86d9f0","#41a5e8","#a2acb3","#e63946"),
                      labels=c("Skin and nasal flora","Aerodigestive flora","Mixed commensal flora","Enteric flora"))+
  expand_limits(x = -.5, y=7.5)+
  scale_y_continuous(breaks=c(-7.5,-5,-2.5,0,2.5,5,7.5), 
                     labels=c("-7.5","-5","-2.5","0","-2.5","-5","-7.5"))+
  xlab("Correlation of Relative Taxon Abundance at Site of Planned Surgical Skin Incision with Median Surgical Operative Level")+
  ylab("-log10(p value)")+
  geom_text_repel(aes(label = taxon),
                   box.padding   = 0.5, 
                   point.padding = .02,
                   segment.color = 'grey50',
                   size=3, 
                   segment.alpha = 0.2,
                    min.segment.length = 0.3,
                    max.overlaps = 100
                   ) +
  theme_classic()+
  scale_size(guide_legend(title="Average Relative Abundance"),breaks = c(.001,.01,.1), labels=c("0.1%","1%","10%"))+
  labs(color = "Anatomic Niche")+
  guides(colour = guide_legend(override.aes = list(size=5)))+ theme(
    panel.background = element_rect(fill = "transparent",
                                    colour = NA_character_), # necessary to avoid drawing panel outline
    panel.grid.major = element_blank(), # remove major grid
    panel.grid.minor = element_blank(), # remove minor grid
    plot.background = element_rect(fill = "transparent",
                                   colour = NA_character_), # necessary to avoid drawing plot outline
    legend.background = element_rect(fill = "transparent",
                                     colour = NA_character_),
    legend.box.background = element_rect(fill = "transparent",
                                         colour = NA_character_),
    legend.key = element_rect(fill = "transparent",
                              colour = NA_character_)
  )

anatomic_niche_plot

if (write_new_figs==T){
ggsave(
  paste0(export_dir,"/figures/anatomic_niche_plot_",runname,".pdf"),
  plot = anatomic_niche_plot,
  scale = 1,
  width = 30,
  height = 20,
  units = "cm",
  dpi = 300,
  bg = "transparent"
)
}



# Prep and export tabular data
taxon.level.summary.continuous.export <- data.frame(
  taxon=taxon.level.summary.continuous$taxon,
  correlation_value=taxon.level.summary.continuous$correlation,
  P_value=taxon.level.summary.continuous$pvalue,
  relative_abundance=taxon.level.summary.continuous$mean_relativeabundance,
  anatomic_reservoir=taxon.level.summary.continuous$anatomic_reservoir
)

addWorksheet(tabular_data, sheetName = "Fig 3")
writeDataTable(tabular_data, sheet = "Fig 3", 
               x = taxon.level.summary.continuous.export,
               colNames = T, rowNames = F)



