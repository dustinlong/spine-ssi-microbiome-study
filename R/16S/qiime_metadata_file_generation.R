library(readxl)
library(readr)
library(jsonlite)
library(RCurl)

# Download sequencing manifest as CSV, import
manifest <- read_csv(paste0(root_dir, "/manifests/Manifest.csv"))

# Subset to just project/study of interest
manifest <- subset(manifest, manifest$Study=="Long_Spine_SSI")

# Parse subject number for sample prefix in aim 1 samples
manifest[manifest$aim=="aim_1" & !is.na(manifest$aim) & manifest$Study=="Long_Spine_SSI", "subject_id"] <- gsub("_.*","", 
                                                      as.data.frame(manifest[manifest$aim=="aim_1" & !is.na(manifest$aim) & manifest$Study=="Long_Spine_SSI","Sample Prefix"])[[1]]
)

# Import clinical metadata from REDCap
redcap_url <-

sequencing_Microbiome_and_HAI.json <- postForm(
  uri=redcap_url,
  token=token,
  content='report',
  format='json',
  type='flat',
  csvDelimiter='',
  report_id='94760',
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
clinical_data <- fromJSON(sequencing_Microbiome_and_HAI.json)


# Manually standardize formatting of DOS labels in clinical_data which were inconsistent for early patients
clinical_data[clinical_data$dos_tube_id=="1",'dos_tube_id']<-"0001"
clinical_data[clinical_data$dos_tube_id=="2",'dos_tube_id']<-"0002"
clinical_data[clinical_data$dos_tube_id=="3",'dos_tube_id']<-"0003"
clinical_data[clinical_data$dos_tube_id=="4",'dos_tube_id']<-"0004"
clinical_data[clinical_data$dos_tube_id=="5",'dos_tube_id']<-"0005"
clinical_data[clinical_data$dos_tube_id=="6",'dos_tube_id']<-"0006"
clinical_data[clinical_data$dos_tube_id=="7",'dos_tube_id']<-"0007"
clinical_data[clinical_data$dos_tube_id=="8",'dos_tube_id']<-"0008"
clinical_data[clinical_data$dos_tube_id=="9",'dos_tube_id']<-"0009"
clinical_data[clinical_data$dos_tube_id=="10",'dos_tube_id']<-"0010"
clinical_data[clinical_data$dos_tube_id=="11",'dos_tube_id']<-"0011"


# For aim 1, link by subject no
manifest[manifest$aim=="aim_1" & !is.na(manifest$aim) & manifest$Study=="Long_Spine_SSI", "merge_ID"] <- manifest[manifest$aim=="aim_1" 
                                                                                                          & !is.na(manifest$aim) 
                                                                                                          & manifest$Study=="Long_Spine_SSI",
                                                                                                            "subject_id"]

# For aim 2, link by Sample Prefix (DOS tube ID)
manifest[manifest$aim=="aim_2" & !is.na(manifest$aim) & manifest$Study=="Long_Spine_SSI", "merge_ID"] <- manifest[manifest$aim=="aim_2" 
                                                                                                        & !is.na(manifest$aim) 
                                                                                                        & manifest$Study=="Long_Spine_SSI",
                                                                                                        "Sample Prefix"]


# Clean any merge_IDs in manifest related to manually recorded label inconsistencies
manifest[manifest$merge_ID=="1109B" & !(is.na(manifest$merge_ID)),'merge_ID'] <-"1109"

# Create parallel merge ID for clinical metadata
clinical_data[clinical_data$aim1==1 , "merge_ID"] <- clinical_data[clinical_data$aim1==1 , "record_id"]
clinical_data[clinical_data$aim1==0 , "merge_ID"] <- clinical_data[clinical_data$aim1==0 , "dos_tube_id"]

# Remove rows of clinical metadata that don't have relevant sample IDs
clinical_data <- subset(clinical_data, !is.na(clinical_data$dos_tube_id))

# Check for duplicates in clinical data before merging
# View(clinical_data[which(duplicated(clinical_data$merge_ID)),])
      
# Merge sample data and clinical metadata
metadata_sheet <- merge(manifest, clinical_data, by.x="merge_ID", by.y ="merge_ID", all.x=TRUE, all.y=FALSE )

# Process/summarize metadata info
metadata_sheet[metadata_sheet$sex==1 & !is.na(metadata_sheet$sex),c("sex_expanded")]<-"female"
metadata_sheet[metadata_sheet$sex==2 & !is.na(metadata_sheet$sex),c("sex_expanded")]<-"male"

metadata_sheet[metadata_sheet$level_category___cervical==1 & !is.na(metadata_sheet$level_category___cervical),c("Operative Level")]<-"cervical"
metadata_sheet[metadata_sheet$level_category___lumbosacral==1 & !is.na(metadata_sheet$level_category___lumbosacral),c("Operative Level")]<-"lumbar"
metadata_sheet[metadata_sheet$level_category___thoracic==1 & !is.na(metadata_sheet$level_category___thoracic),c("Operative Level")]<-"thoracic"

# Standardize source labels
metadata_sheet[metadata_sheet$Source=="Rectal_libredo" & !is.na(metadata_sheet$Source),c("Source")]<-"Rectal"
metadata_sheet[metadata_sheet$Source=="Nasal_libredo" & !is.na(metadata_sheet$Source),c("Source")]<-"Nasal"
metadata_sheet[metadata_sheet$Source=="control" & !is.na(metadata_sheet$Source),c("Source")]<-"Control"
metadata_sheet[metadata_sheet$Source=="skin" & !is.na(metadata_sheet$Source),c("Source")]<-"Surgical Site"
metadata_sheet[metadata_sheet$Source=="dos_nasal" & !is.na(metadata_sheet$Source),c("Source")]<-"Nasal"
metadata_sheet[metadata_sheet$Source=="dos_rectal" & !is.na(metadata_sheet$Source),c("Source")]<-"Rectal"
metadata_sheet[metadata_sheet$Source=="dos_skin" & !is.na(metadata_sheet$Source),c("Source")]<-"Surgical Site"
metadata_sheet[metadata_sheet$Source=="dos_skin_2" & !is.na(metadata_sheet$Source),c("Source")]<-"Surgical Site" # Posterior
metadata_sheet[metadata_sheet$Source=="dos_skin_1" & !is.na(metadata_sheet$Source),c("Source")]<-"Surgical Site Anterior"

# Standardize sample names
metadata_sheet$`Sample Name`<- gsub("libredo","redo", metadata_sheet$`Sample Name`)

# Derive combinatorial features
metadata_sheet[metadata_sheet$Source=="Surgical Site" & !is.na(metadata_sheet$Source)
               & metadata_sheet$`Operative Level`=="cervical" & !is.na(metadata_sheet$`Operative Level`)
               ,c("Source - Level")]<-"Surgical Site Cervical"

metadata_sheet[metadata_sheet$Source=="Surgical Site" & !is.na(metadata_sheet$Source)
               & metadata_sheet$`Operative Level`=="thoracic" & !is.na(metadata_sheet$`Operative Level`)
               ,c("Source - Level")]<-"Surgical Site Thoracic"

metadata_sheet[metadata_sheet$Source=="Surgical Site" & !is.na(metadata_sheet$Source)
               & metadata_sheet$`Operative Level`=="lumbar" & !is.na(metadata_sheet$`Operative Level`)
               ,c("Source - Level")]<-"Surgical Site Lumbar"


metadata_sheet[metadata_sheet$Source=="Nasal" & !is.na(metadata_sheet$Source)
               & metadata_sheet$`Operative Level`=="cervical" & !is.na(metadata_sheet$`Operative Level`)
               ,c("Source - Level")]<-"Nasal Cervical"

metadata_sheet[metadata_sheet$Source=="Nasal" & !is.na(metadata_sheet$Source)
               & metadata_sheet$`Operative Level`=="thoracic" & !is.na(metadata_sheet$`Operative Level`)
               ,c("Source - Level")]<-"Nasal Thoracic"

metadata_sheet[metadata_sheet$Source=="Nasal" & !is.na(metadata_sheet$Source)
               & metadata_sheet$`Operative Level`=="lumbar" & !is.na(metadata_sheet$`Operative Level`)
               ,c("Source - Level")]<-"Nasal Lumbar"


metadata_sheet[metadata_sheet$Source=="Rectal" & !is.na(metadata_sheet$Source)
               & metadata_sheet$`Operative Level`=="cervical" & !is.na(metadata_sheet$`Operative Level`)
               ,c("Source - Level")]<-"Rectal Cervical"

metadata_sheet[metadata_sheet$Source=="Rectal" & !is.na(metadata_sheet$Source)
               & metadata_sheet$`Operative Level`=="thoracic" & !is.na(metadata_sheet$`Operative Level`)
               ,c("Source - Level")]<-"Rectal Thoracic"

metadata_sheet[metadata_sheet$Source=="Rectal" & !is.na(metadata_sheet$Source)
               & metadata_sheet$`Operative Level`=="lumbar" & !is.na(metadata_sheet$`Operative Level`)
               ,c("Source - Level")]<-"Rectal Lumbar"


# Rename columns
colnames(metadata_sheet)[which(colnames(metadata_sheet) %in% c("Sample Name"))]<-"sample-id"
colnames(metadata_sheet)[which(colnames(metadata_sheet) %in% c("sex_expanded"))]<-"Sex"

# Standard reformatting for qiime--these restrictions apply only to sample-id?
metadata_sheet$`sample-id` <- gsub("_","-", metadata_sheet$`sample-id`)
metadata_sheet$`sample-id` <- gsub(" ","-", metadata_sheet$`sample-id`)
metadata_sheet$`sample-id` <- gsub("#","", metadata_sheet$`sample-id`)

save(metadata_sheet, file=paste0(root_dir,"/16S_metadata.Rda"))

# Pick fields to export for qiime
metadata_sheet.export <- metadata_sheet[,c("sample-id","Study","Source", "Sex", "Operative Level" , "Source - Level" , "aim1", "record_id","dos_tube_id")]

# Create qiime header
metadata_sheet.export.classified <-
rbind(   c("#q2:types", # Sample ID
          "categorical", #Study
          "categorical", # Source
          "categorical", # Sex expanded
          "categorical", # Operative level
          "categorical", # Source - Level
          "categorical", # aim 1
         "categorical", #  merge_ID
         "categorical", # record_id
         "categorical") # dos_tube_id
         , metadata_sheet.export)

