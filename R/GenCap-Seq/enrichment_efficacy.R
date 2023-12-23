library(scales)
library(readr)

write_plots<-F

# Set Kraken report parameters
kraken_report_columns <- c("read_percent","read_count_recursive","read_count_level","taxonomic_level","taxon_id","taxon_name")
kraken_report_dir <- paste0(shared_dir, "/genome_capture_seq/kraken_reports/")
 
library(readxl)
GCS_Sample_Table <- read_excel(paste0(shared_dir, "spine_microbiome/GCS Sample Table.xlsx"))

# Create summary data holder
enrichment.report<-data.frame()

for(row.temp in rownames(GCS_Sample_Table)){
  subject<-as.character(GCS_Sample_Table[row.temp,"Patient ID"])
  sample<-as.character(GCS_Sample_Table[row.temp,"Library prefix"])
  source<-as.character(GCS_Sample_Table[row.temp,"Sample Source"])
  taxon<-as.character(GCS_Sample_Table[row.temp,"Targeted Taxon"])
  batch<-as.character(GCS_Sample_Table[row.temp,"Batch Date"])

  # Ignore subject without unhyb comparison available
  if( !(subject %in% c("65"))){

# Identify pre-enrichment report
preenrichment_kraken_report<-
  list.files(path=paste0(kraken_report_dir,"subject_",subject,
                         "/sample_",sample,
                         "/",taxon,
                         "/",source,
                         "/preenrichment"), full.names=T, recursive = T)

# Identify pre-enrichment report
postenrichment_kraken_report <- list.files(path=paste0(kraken_report_dir,"subject_",subject,
                       "/sample_",sample,
                       "/",taxon,
                       "/",source,
                       "/postenrichment"), full.names=T, recursive = T)


# Import Kraken reports
enriched <- read_delim(
  postenrichment_kraken_report
  ,
                                                           delim = "\t", escape_double = FALSE,
                                                           col_names = kraken_report_columns, trim_ws = TRUE)

unenriched <- read_delim(
  preenrichment_kraken_report
  ,
                                                           delim = "\t", escape_double = FALSE,
                                                           col_names = kraken_report_columns, trim_ws = TRUE)


total_reads.enriched <- sum(enriched$read_count_level, na.rm = T)
total_reads.unenriched <- sum(unenriched$read_count_level, na.rm = T)

# Percent of overall reads in unenriched / enriched
percent(total_reads.unenriched/total_reads.enriched)

# Merge pre-post-enrichment tables
enrichment_comparison <- merge(unenriched, enriched, by=c("taxon_id","taxon_name"), all.x = T,  all.y = T, suffixes=c(".unenriched",".enriched"))

# Set NAs after merge to zero (no reads)
enrichment_comparison[is.na(enrichment_comparison$read_count_recursive.unenriched),"read_count_recursive.unenriched"]<-0
enrichment_comparison[is.na(enrichment_comparison$read_count_level.unenriched),"read_count_level.unenriched"]<-0

enrichment_comparison[is.na(enrichment_comparison$read_count_level.enriched),"read_count_level.enriched"]<-0
enrichment_comparison[is.na(enrichment_comparison$read_count_level.enriched),"read_count_level.enriched"]<-0

# Add pseudocount
pesudocount <-0.1
enrichment_comparison$proportion_exact_level.unenriched <- (enrichment_comparison$read_count_level.unenriched+pesudocount)/total_reads.unenriched
enrichment_comparison$proportion_exact_recursive.unenriched <- (enrichment_comparison$read_count_recursive.unenriched+pesudocount)/total_reads.unenriched

enrichment_comparison$proportion_exact_level.enriched <- (enrichment_comparison$read_count_level.enriched+pesudocount)/total_reads.enriched
enrichment_comparison$proportion_exact_recursive.enriched <- (enrichment_comparison$read_count_recursive.enriched+pesudocount)/total_reads.enriched

# Summary table
enrichment.summary<-data.frame(
  taxon=enrichment_comparison$taxon_name, 
  fold_enrichment=enrichment_comparison$proportion_exact_recursive.enriched/enrichment_comparison$proportion_exact_recursive.unenriched,
  unenriched_read_count=enrichment_comparison$read_count_recursive.unenriched,
  enriched_read_count=enrichment_comparison$read_count_recursive.enriched)

enrichment.summary <- enrichment.summary[order(-enrichment.summary$fold_enrichment),]

enrichment.report<-rbind(enrichment.report,
      data.frame(
        subject=subject,
        sample=sample,
        source=source,
        taxon=taxon,
        batch=batch,
        read_count.enriched.total.million=total_reads.enriched/1000000,
        read_count.enriched.targettaxon=enrichment.summary[enrichment.summary$taxon==taxon,"enriched_read_count"],
        relative_abundance.unenriched.targettaxon=enrichment_comparison[enrichment_comparison$taxon_name==taxon, "proportion_exact_recursive.unenriched"],
        fold_enrichment.targettaxon=round(enrichment.summary[enrichment.summary$taxon==taxon,"fold_enrichment"],2)
        )
      )
}
} # End for sample


# Subset samples without protocol failure
enrichment.report.included<-subset(enrichment.report, enrichment.report$fold_enrichment.targettaxon > 0.25)

# Label according to protocol type
enrichment.report.included[enrichment.report.included$batch %in% c("Aug 22","April 22"),"GCS_protocol"]<-"modified"
enrichment.report.included[!enrichment.report.included$batch %in% c("Aug 22","April 22"),"GCS_protocol"]<-"standard"

# Plot enrichment efficacy by baseline abundance
enrichment.report.plot<-ggplot(enrichment.report.included, aes(x=1/relative_abundance.unenriched.targettaxon, y=fold_enrichment.targettaxon))+
  geom_point()+
  theme_minimal()+
  geom_smooth(method="loess")+
  geom_point(aes(
                color=source, shape=GCS_protocol
                ), size=4)+
  scale_x_log10(name="Relative Abundance",
                breaks = c(10,100,1000,10000,100000,1000000,10000000),
                labels=c("1/10","1/100","1/1,000","1/10,000","1/1000,000","1/1,000,000","1/10,000,000")
                )+
  scale_y_log10(name="Target-Taxon Enrichment",
                breaks = c(0.5,1,2,5,10,25,50,100,250,500,1000,2500,10000,25000),
                labels=c("0.5x","1x","2x","5x","10x","25x","50x","100x","250x","500x","1,000x","2,500x","10,000x","25,000x"))+
  labs(shape="GenCap-Seq Protocol", color="Sample Source"
       )+
  scale_color_manual(values=c("#00b4d8", "#5a189a", "#f4d35e"))

enrichment.report.plot

if(write_plots==T){
ggsave(
       paste0(export_dir,"/figures/GenCap-Seq Enrichment_",runname,".pdf"),
       plot=enrichment.report.plot,
       height = 6,
       width=12,
       units = "in")
}

# Stage tabular data export
enrichment.report.export <- data.frame(
                                       Sample_ID=enrichment.report.included$sample,
                                       Source_type=enrichment.report.included$source,
                                       Protocol=enrichment.report.included$GCS_protocol,
                                       Baseline_abundance=enrichment.report.included$relative_abundance.unenriched.targettaxon,
                                       Fold_enrichment=enrichment.report.included$fold_enrichment.targettaxon)

addWorksheet(tabular_data, sheetName = "Fig S3")

writeDataTable(tabular_data, sheet = "Fig S3", 
               x = enrichment.report.export,
               colNames = T, rowNames = F)






