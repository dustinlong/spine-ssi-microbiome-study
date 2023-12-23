library(readxl)
library(readr)
library(jsonlite)
library(RCurl)

# Create list of SSI isolates from repository 1
isolate_list_biorepository1 <- c("A0001"
                      ,"A0002"
                      ,"A0003"
                      ,"A0005"
                      ,"A0006"
                      ,"A0004"
                      ,"A0014"
                      ,"A0013"
                      ,"A0015"
                      ,"A0017"
                      ,"A0018"
                      ,"A0011"
                      ,"A0019"
                      ,"A0023"
                      ,"A0022"
                      ,"A0021"
                      ,"A0024"
                      ,"A0029"
                      ,"A0030"
                      ,"A0016"
                      ,"A0012"
                      ,"A0009"
                      ,"A0010"
                      ,"A0007"
                      ,"A0027"
                      ,"A0025"
                      ,"A0026"
                      ,"A0028"
                      ,"A0020")

# Import clinical metadata from REDCap
redcap_url <- 

Molecular_Epi_of_HAI.json <- postForm(
  uri=redcap_url,
  token=token,
  content='report',
  format='json',
  type='flat',
  csvDelimiter='',
  report_id='217580',
  rawOrLabel='raw',
  rawOrLabelHeaders='raw',
  exportCheckboxLabel='false',
  exportSurveyFields='false',
  exportDataAccessGroups='false',
  returnFormat='json'
)

# Clear token from cache
rm(token)

# Convert JSON to dataframe
clinical_isolates.phase1 <- fromJSON(Molecular_Epi_of_HAI.json)
clinical_isolates.phase1$record_id <- paste0("molepi_",clinical_isolates.phase1$record_id)


# Import clinical metadata from REDCap database 2
sequencing_Microbiome_and_HAI.json <- postForm(
  uri=redcap_url,
  token=token,
  content='report',
  format='json',
  type='flat',
  csvDelimiter='',
  report_id='42920',
  rawOrLabel='raw',
  rawOrLabelHeaders='raw',
  exportCheckboxLabel='false',
  exportSurveyFields='false',
  exportDataAccessGroups='false',
  returnFormat='json'
)

# Clear token from cache
rm(token)

# Convert JSON to dataframe
clinical_data.gcs <- fromJSON(sequencing_Microbiome_and_HAI.json)

clinical_data.gcs$record_id <- paste0("GCS_",clinical_data.gcs$record_id)


# Combine clinical SSI isolate values from multiple fields in REDCap form
ssi_isolate_stock_ids<-c(clinical_isolates.phase1$specimen_location,
clinical_isolates.phase1$specimen_location_2,
clinical_isolates.phase1$specimen_location_2_3,
clinical_isolates.phase1$specimen_location_2_3_4,
clinical_isolates.phase1$specimen_location_2_3_4_5)

# Parse isolate ID
ssi_isolate_stock_ids<-gsub(".*(A0[0-9]+).*","\\1",ssi_isolate_stock_ids)
ssi_isolate_stock_ids<-subset(ssi_isolate_stock_ids, ssi_isolate_stock_ids!="")


# Manually add stock IDs in cases with multiple isolates from a single culture not captured using approach above
ssi_isolate_stock_ids<-c(ssi_isolate_stock_ids,"A0043")

# Import manifest of cases from SSI microbiome sub-aim
GCS_Pair_Table_v4 <- read_excel("./spine_microbiome/GCS Pair Table v4.xlsx")

ssi_isolate_stock_ids<-c(ssi_isolate_stock_ids,subset(GCS_Pair_Table_v4$`SSI Isolate ID`, !is.na(GCS_Pair_Table_v4$`SSI Isolate ID`)))


extraneous_isolates<-c("A0051",
                   "A0093",
                   "A0077",
                   "A0092",
                   "A0086",
                   "A0125",
                   "A0055",
                   "A0056",
                   "A0060",
                   "A0064",
                   "A0072",
                   "A0088",
                   "A0091",
                   "A0113")

# Remove extraneous culture isolates (not SSI pathogens)
ssi_isolate_stock_ids <- subset(ssi_isolate_stock_ids, !(ssi_isolate_stock_ids %in% extraneous_isolates))

# For each included SSI isolate, look up the corresponding patient
ssi_isolate_summary<-data.frame()
for (isolate.temp in ssi_isolate_stock_ids){
 # Check each culture field in REDCap form for a match
  subject.temp<-c(clinical_data.gcs[grepl(isolate.temp,clinical_data.gcs$ssi_1_location), "record_id"],
          clinical_data.gcs[grepl(isolate.temp,clinical_data.gcs$ssi_2_location), "record_id"],
          clinical_data.gcs[grepl(isolate.temp,clinical_data.gcs$ssi_3_location), "record_id"],
          clinical_data.gcs[grepl(isolate.temp,clinical_data.gcs$ssi_4_location), "record_id"],
          clinical_data.gcs[grepl(isolate.temp,clinical_data.gcs$ssi_5_location), "record_id"],
        
          clinical_isolates.phase1[grepl(isolate.temp,clinical_isolates.phase1$specimen_location), "record_id"],
          clinical_isolates.phase1[grepl(isolate.temp,clinical_isolates.phase1$specimen_location_2), "record_id"],
          clinical_isolates.phase1[grepl(isolate.temp,clinical_isolates.phase1$specimen_location_2_3), "record_id"],
          clinical_isolates.phase1[grepl(isolate.temp,clinical_isolates.phase1$specimen_location_2_3_4), "record_id"],
          clinical_isolates.phase1[grepl(isolate.temp,clinical_isolates.phase1$specimen_location_2_3_4_5), "record_id"]
          )
  
  
  ssi_isolate_summary<-rbind(ssi_isolate_summary, 
                             data.frame(subject=subject.temp,
                                        isolate=isolate.temp))
}


library(readxl)
# Import assembly info
SSI_WGS_stats <- read_excel("./sequencing_results/WGS_summary_data/200115_surg_inf_assembly_stats.xlsx")

ssi_isolate_summary.annotated <-   merge(ssi_isolate_summary, SSI_WGS_stats, by.x="isolate", by.y ="Sample_ID...2" , all.x = T, all.y = F)
ssi_isolate_summary.annotated <-subset(ssi_isolate_summary.annotated, !is.na(ssi_isolate_summary.annotated$estimated_coverage))

# Generate summary dataframes for loop
wgs_distance_summary <- data.frame()
included_isolates<-data.frame()
within_vals_summary<-data.frame()

for(organism_temp in c("S. aureus",
                       "S. epidermidis",
                       "E. coli",
                       "P. mirabilis",
                       "P. aeruginosa",
                       "C. striatum",
                       "C. koseri",
                       "E. faecalis",
                       "S. capitis",
                       "K. pneumoniae",
                       "S. lugdunensis"
                       )){
  
if(organism_temp=="E. coli"){

# For each species, define WGS distance grid files, manually store distances from intra-patient pairs and exclude these isolates from inter-patient pair grid
  
# "Escherichia coli"
distance_grid.temp<-"./sequencing_results/WGS_epi_output/Escherichia_coli/Epidemiology/Epidemiology_distance_grid.txt"
exclusion_list.temp <-c("A0009"	,"A0012" ,"A0044" ,"A0126",	"A0127"	,"A0128")
within_vals<-c(51, 31,40, 81, 82, 41,
               53,
               41, 11, 11 ,4 )

}

# S aureus
  if(organism_temp=="S. aureus"){
    distance_grid.temp<-"./sequencing_results/WGS_epi_output/Staphylococcus_aureus/Epidemiology/Epidemiology_distance_grid.txt"
exclusion_list.temp <-c("A0079","A0084","A0067")
within_vals<-c(99, 19, 15 )
  }
  
# Proteus mirabilis
  if(organism_temp=="P. mirabilis"){
    
    distance_grid.temp<-"./sequencing_results/WGS_epi_output/Proteus_mirabilis/Epidemiology/Epidemiology_distance_grid.txt"
exclusion_list.temp <-c("A0007","A0075","A0059","A0087")   

within_vals<-c(13 ,19 ,12 ,33)
  }
  
# Staphylococcus epidermidis
  if(organism_temp=="S. epidermidis"){
    
    distance_grid.temp<-"./sequencing_results/WGS_epi_output/Staphylococcus_epidermidis/Epidemiology/Epidemiology_distance_grid.txt"
exclusion_list.temp <-NA
within_vals<-NA
  }
  
  
# Pseudomonas aeruginosa
  if(organism_temp=="P. aeruginosa"){
    distance_grid.temp<-"./sequencing_results/WGS_epi_output/Pseudomonas_aeruginosa/Epidemiology/Epidemiology_distance_grid.txt"
exclusion_list.temp <-c("A0125")
within_vals<-c(30)
}


  # S lugdunensis
  if(organism_temp=="S. lugdunensis"){
    distance_grid.temp<-"./sequencing_results/WGS_epi_output/Staphylococcus_lugdunensis/Epidemiology/Epidemiology_distance_grid.txt"
    exclusion_list.temp <-NA
    within_vals<-NA
  }

  
  
  # S capitis
  if(organism_temp=="S. capitis"){
    distance_grid.temp<-"./sequencing_results/WGS_epi_output/Staphylococcus_capitis/Epidemiology/Epidemiology_distance_grid.txt"
    exclusion_list.temp <-NA
    within_vals<-NA
  }
  
  
  
  # K  pneumoniae
  if(organism_temp=="K. pneumoniae"){
    distance_grid.temp<-"./sequencing_results/WGS_epi_output/Klebsiella_pneumoniae/Epidemiology/Epidemiology_distance_grid.txt"
    exclusion_list.temp <- NA
    within_vals<-c(13 ,19 ,12 ,33)
  }

  
  
  
  # E  faecalis
  if(organism_temp=="E. faecalis"){
    distance_grid.temp<-"./sequencing_results/WGS_epi_output/Enterococcus_faecalis/Epidemiology/Epidemiology_distance_grid.txt"
    exclusion_list.temp <-NA
    within_vals<-NA
  }
  
  
  # C  striatum
  if(organism_temp=="C. striatum"){
    distance_grid.temp<-"./sequencing_results/WGS_epi_output/Corynebacterium_striatum/Epidemiology/Epidemiology_distance_grid.txt"
    exclusion_list.temp <-NA
    within_vals<-NA
  }
  
  
  
  
  # C  koseri
  if(organism_temp=="C. koseri"){
    distance_grid.temp<-"./sequencing_results/WGS_epi_output/Citrobacter_koseri/Epidemiology/Epidemiology_distance_grid.txt"
    exclusion_list.temp <-NA
    within_vals<-NA
  }
  
  
  
  
# Import distance grid
Epidemiology_distance_grid.temp <- read_delim(distance_grid.temp, 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)


# Clean, format
Epidemiology_distance_grid.temp<-as.data.frame(Epidemiology_distance_grid.temp)
rownames(Epidemiology_distance_grid.temp)<-Epidemiology_distance_grid.temp$...1

Epidemiology_distance_grid.temp<-Epidemiology_distance_grid.temp[,
                                which(colnames(Epidemiology_distance_grid.temp) %in% rownames(Epidemiology_distance_grid.temp))
]


# Remove intra-patient samples in exclusion list
Epidemiology_distance_grid.temp.subset<-Epidemiology_distance_grid.temp[which(
                                                    rownames(Epidemiology_distance_grid.temp) %in% ssi_isolate_summary.annotated$isolate &
                                       !(rownames(Epidemiology_distance_grid.temp) %in% exclusion_list.temp)),
                                      
                                which(
                                  colnames(Epidemiology_distance_grid.temp) %in% ssi_isolate_summary.annotated$isolate &
                                       !(colnames(Epidemiology_distance_grid.temp) %in% exclusion_list.temp))
                                
                                ]


# Convert to matrix
matrix.temp<-as.matrix(Epidemiology_distance_grid.temp.subset)

# Append to summary table
included_isolates<-rbind(included_isolates,
                     data.frame(
                       isolate=colnames(Epidemiology_distance_grid.temp.subset),
                       species=organism_temp
                       )
)

# Subset to non-redundant upper triange  of distance grid, excluding diagonal (self-referent)
distance_values<-data.frame(
  distance=matrix.temp[upper.tri(matrix.temp, diag = F)],
  species=organism_temp)

# Store values
wgs_distance_summary<-rbind(wgs_distance_summary,
                            distance_values)


within_vals_summary<-rbind(within_vals_summary,
                   data.frame(
                     distance=within_vals,
                     species=organism_temp
                   ))

within_vals_summary<-subset(within_vals_summary, !is.na(within_vals_summary$distance))


} # End for organism

# Standardize species nomenclature
included_isolates[included_isolates$species=="S. aureus","species"]<-"Staphylococcus aureus"
included_isolates[included_isolates$species=="S. epidermidis","species"]<-"Staphylococcus epidermidis"
included_isolates[included_isolates$species=="E. coli","species"]<-"Escherichia coli"
included_isolates[included_isolates$species=="P. mirabilis","species"]<-"Proteus mirabilis"
included_isolates[included_isolates$species=="P. aeruginosa","species"]<-"Pseudomonas aeruginosa"
included_isolates[included_isolates$species=="C. striatum","species"]<-"Corynebacterium striatum"
included_isolates[included_isolates$species=="C. koseri","species"]<-"Citrobacter koseri"
included_isolates[included_isolates$species=="E. faecalis","species"]<-"Enterococcus faecalis"
included_isolates[included_isolates$species=="S. capitis","species"]<-"Staphylococcus capitis"
included_isolates[included_isolates$species=="K. pneumoniae","species"]<-"Klebsiella pneumoniae"
included_isolates[included_isolates$species=="S. lugdunensis","species"]<-"Staphylococcus lugdunensis"

# Label pairs as inter-patient or intra-patient, italicize species
wgs_distance_summary$species<-paste0("inter-patient: *",wgs_distance_summary$species,"*")
within_vals_summary$species<-paste0("intra-patient: *",within_vals_summary$species,"*")

# Rbind for combined plotting
distance_summary<-rbind(within_vals_summary,wgs_distance_summary)

library(ggplot2)
library(ggtext)

# Plot WGS distance of intra- vs inter-patient same species pairs
mol_epi<-ggplot(distance_summary, aes(x=distance, fill=species))+
  geom_histogram(bins = 50)+
  scale_y_sqrt(name="Count of Same-Species Pairs Separated by Indicated Genetic Distance", breaks=c(0,1,5,10,25,50,100,200,400,1000))+
  scale_x_log10(name="Same-Species SSI Isolate Genetic Distance (Whole-Genome Variants)",
    limits=c(1,10^5), breaks=c(0,1,2,5,10,25,50,100,250,500,1000,2500,10000,50000),
                labels=c("0","1","2","5","10","25","50","100","250","500","1,000","2,500","10,000","50,000"))+
  geom_vline(aes(xintercept = 36), linetype = "dashed", alpha=0.5)+ # https://pubmed.ncbi.nlm.nih.gov/33296290/
  geom_vline(aes(xintercept = 40), linetype = "dashed", alpha=0.5)+ # https://pubmed.ncbi.nlm.nih.gov/26230489/
  geom_vline(aes(xintercept = 25), linetype = "dashed", alpha=0.5)+ # https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(20)30149-X/fulltext
  scale_fill_manual(name="Species", values=c(
    "#d9ed92", "#b5e48c", "#90be6d", "#76c893", "#52b69a", "#34a0a4", "#2196cc", "#1a759f", "#1e6091", "#184e77","#183877"
    ,"#ff85a1", "#ffb700",  "#f7b267",  "#ff7b00", "#f25c54"
  )
  , aes(labels=species)
  )+
  theme_minimal()+
  theme(legend.text = element_markdown(),panel.grid.minor = element_blank())


mol_epi

ggsave(paste0(export_dir,"/figures/genomic_surveillance_plots", runname, ".pdf"),
       plot=mol_epi,
       width=10,
       height=6,
       units="in")

# Stage data for export in tabular format
molecular_epi.export <- data.frame(
  Patient_species_group=gsub("\\*","",distance_summary$species),
  WGS_distance=distance_summary$distance
)

addWorksheet(tabular_data, sheetName = "Fig 7")
writeDataTable(tabular_data, sheet = "Fig 7", 
               x = molecular_epi.export,
               colNames = T, rowNames = F)

