library(readr)
library(vcfR)
library(tidyverse)
library(readxl)


max_snapshot_loci <- 50000
pileupsites <- 50000
parseandrun_IGVcommand <- F
export_reports <- F
export_plots <- F

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Rsamtools")
# 
# library(Rsamtools)

GCS_Pair_Table <- read_excel(paste0(onedrive_dir,"spine_microbiome/GCS Pair Table v2.xlsx"))

GCS_Pair_Table<-subset(GCS_Pair_Table, !is.na(GCS_Pair_Table$`Case Number`))


allele_support.summary <-data.frame()

# For each microbiome-SSI isolate pair in manifest
for(row.temp in as.numeric(rownames(GCS_Pair_Table))){
  # Wait and clean up environment
  Sys.sleep(10)
  rm(list=setdiff(ls(), c("pileupsites","allele_support.summary","GCS_Pair_Table","row.temp", "share_dir","github_dir","host", "onedrive_dir", "max_snapshot_loci","parseandrun_IGVcommand")))

  # Store case-level parameters
  subject<-as.character(GCS_Pair_Table[as.character(row.temp),"Patient ID"])
  print(subject)
  
  sample<-as.character(GCS_Pair_Table[as.character(row.temp),"Library ID"])
  print(sample)
  
  reference_genome<-as.character(GCS_Pair_Table [as.character(row.temp),"SSI Isolate Closest Reference Genome"])
  print(reference_genome)
  
  manual_call<-as.character(GCS_Pair_Table [as.character(row.temp),"SSI Strain Identified in Preoperative Microbiome Specimens?"])
  print(manual_call)
  
  ssi_isolate<-as.character(GCS_Pair_Table [as.character(row.temp),"SSI Isolate ID"])
  print(ssi_isolate)
  
# Get list of organism directories for this pair
organisms<-
  list.files(path=paste0(share_dir,"sequencing_results/genome_capture_seq/variants_alignments/subject_",subject,
                         "/sample_",sample
                         ), full.names=T)

# Identify organism dir for the given reference in this loop
organism_dir<-
  organisms[grepl(reference_genome,organisms)]

# List all the files for this sample-probe set
list<-list.files(path=organism_dir, full.names=T, recursive = T)

# Get SSI VCF from this list
ssi_vcf <- list[grepl("ssi_alignments.*vcf$",list, perl = T)]

# Get SSI BAM from this list
ssi_bam <- list[grepl(paste0(".*","ssi_alignments",".*bam$"),list, perl = T)]

# Identify reference genome sequence file
reference_file <- list[grepl(".*fa$",list, perl = T)]

# Get sample sources for the patient
sample_names <-
  list.dirs(path=paste0(organism_dir, "/same_patient/"), full.names=F, recursive = F)

gcs_samples <- data.frame(sample_name=sample_names)

# For each preoperative swab from this pair
for (sample_name.temp in gcs_samples$sample_name){
  
  # Identify BAM file
  gcs_samples[gcs_samples$sample_name==sample_name.temp,"bam_file"] <-  list[grepl(paste0(".*/",sample_name.temp,"/.*bam$"),list, perl = T)]
  
}


# Import and parse SSI isolate WGS VCF
ssi_import <- read.vcfR(file= ssi_vcf, verbose = F)

ssi_genotypes <- as.data.frame(extract.gt(ssi_import, element = "GT"))

ssi_alleles <- as.data.frame(ssi_import@fix)

ssi_variants <- data.frame("chromosome"=ssi_alleles$CHROM, 
                           "position"=as.numeric(ssi_alleles$POS), 
                           "reference_allele"=ssi_alleles$REF,
                           "alternate_allele"=ssi_alleles$ALT,
                           "genotype.ssi"=ssi_genotypes[[1]])



# Import GenCap-Seq depth of coverage info
for (gcs_sample.row in rownames(gcs_samples)){
  sample_name.temp <- gcs_samples[gcs_sample.row,"sample_name"]
  
  print(sample_name.temp)
  
  # Check if coverage file exists, if not generate
  if(!(file.exists(paste0(gcs_samples[gcs_sample.row,"bam_file"],".coverage")))){
    print("Generating coverage file")
    system(paste0("samtools depth -a ",
                  paste0(gcs_samples[gcs_sample.row,"bam_file"])
                  ," > "
                  ,paste0(gcs_samples[gcs_sample.row,"bam_file"],".coverage"))
    )
  }
  
  # Import coverage file
  coverage_file.temp <-    as.data.frame(read_delim(paste0(gcs_samples[gcs_sample.row,"bam_file"],".coverage"), 
                                      delim = "\t", escape_double = FALSE, 
                                      col_names = FALSE, trim_ws = TRUE, col_types = cols(),
                                      progress = F))

  
  
  
  if (nrow(coverage_file.temp)>0){
  # Clean coverage table
  coverage_file.temp<-coverage_file.temp[-1]
  colnames(coverage_file.temp)<- c("position",paste0(sample_name.temp,"_coverage"))
  
  # Name and store temp coverage file
  assign(paste0(sample_name.temp,".coverage"), coverage_file.temp)
    
}

}


# Merge coverage data into one df
coverage.combined <- Reduce(function(x, y) merge(x, y, by="position", all.x=T, all.y=F), mget(
  c("ssi_variants",ls(pattern="\\.coverage$"))
  ))




# Review SNP sites w coverage

# All sites with variant in SSI
all_overlaping_sites <- 
  subset(coverage.combined, 
           coverage.combined$genotype.ssi==1
)




samples<-colnames(all_overlaping_sites[6:ncol(all_overlaping_sites)])


# Count number of samples with coverage
all_overlaping_sites$samples_covered <- rowSums(all_overlaping_sites[,which(
  colnames(all_overlaping_sites) %in% paste0(gcs_samples$sample_name,"_coverage")
  & !(grepl("postop",colnames(all_overlaping_sites)))
)
]>0)


# Count total reads
all_overlaping_sites$total_coverage <- rowSums(
  all_overlaping_sites[,which(
  colnames(all_overlaping_sites) %in% paste0(gcs_samples$sample_name,"_coverage")
  & !(grepl("postop",colnames(all_overlaping_sites)))
)]
)

# Store strain-informative loci
all_overlaping_sites<-all_overlaping_sites[order(-all_overlaping_sites$samples_covered, -all_overlaping_sites$total_coverage),]



# IGV alignment batch script generation
 if (parseandrun_IGVcommand==T){


raw_text_template<-read_file( paste0(github_dir,"SSI_GenomeCaptureSeq/IGV_batch_scripting/IGV_batch_template.txt"))

raw_text <- raw_text_template


raw_text<-gsub("<REFERENCE_GENOME>",
               reference_file,
               raw_text)


raw_text<-gsub("<SSI_BAM>",
               ssi_bam,
               raw_text)

raw_text<-gsub("<SNAPSHOT_PATH>",
               paste0(share_dir,"/genome_capture_seq/variant_alignment_snapshots/","subject",subject,"_sample",sample,"_",reference_genome),
               raw_text)


# Set track order
gcs_samples$order <-100

gcs_samples[gcs_samples$sample_name=="skin","order"]<-1
gcs_samples[gcs_samples$sample_name=="nares","order"]<-2
gcs_samples[gcs_samples$sample_name=="nasal","order"]<-2
gcs_samples[gcs_samples$sample_name=="rectal","order"]<-3


  
for(sample_source in gcs_samples[order(gcs_samples$order),"sample_name"]){
  print(sample_source)
  
  raw_text<- paste(raw_text, 
                      paste0("load ", gcs_samples[gcs_samples$sample_name==sample_source,"bam_file"], " name=",sample_source)
                      , sep = "\n")
  
}

raw_text<- paste(raw_text, 
                 "squish"
                 , sep = "\n")

any_overlaping_sites<-subset(all_overlaping_sites, all_overlaping_sites$samples_covered >=1)

for(covered_variant_pos in any_overlaping_sites$position[1:min(nrow(any_overlaping_sites),max_snapshot_loci)]){
  print(covered_variant_pos)

  raw_text<- paste(raw_text, 
                   paste0("goto ", reference_genome,":", covered_variant_pos)
                   , sep = "\n")
  
  raw_text<- paste(raw_text, 
                   "snapshot",
                   "squish"
                   , sep = "\n")
  
  
  
}


raw_text<- paste(raw_text, 
                 paste0("saveSession ",share_dir,"UW/research/projects/SSI/HMC_spine_SSI/sequencing_results/genome_capture_seq/variant_alignment_snapshots/","subject",subject,"_sample",sample,"_",reference_genome,"/session")
                 , sep = "\n")



raw_text<- paste(raw_text, 
                 "exit"
                 , sep = "\n")


write(raw_text,                                            # Write new line to file
      file = paste0(github_dir,"SSI_GenomeCaptureSeq/IGV_batch_scripting/IGV_batch_current.txt"),
      append = FALSE)



# Run IGV batch command
system(paste0("/Users/dustinlong/Downloads/IGV_2.14.0/igv.command -b ",
              github_dir, "SSI_GenomeCaptureSeq/IGV_batch_scripting/IGV_batch_current.txt"))


}


# Navigate BAM for GenCap-Seq metagenomic reads
 library(Rsamtools)

 PileupParam(max_depth=250, min_base_quality=13, min_mapq=0,
            min_nucleotide_depth=1, min_minor_allele_depth=0,
            distinguish_strands=FALSE, distinguish_nucleotides=TRUE,
            ignore_query_Ns=TRUE, 
            include_deletions=FALSE, include_insertions=FALSE,
            left_bins=NULL, query_bins=NULL, cycle_bins=NULL)

 for(row.temp in as.numeric(rownames(gcs_samples))){

  source.temp<-gcs_samples[as.character(row.temp),"sample_name"]

  pileup.temp<-pileup(gcs_samples[as.character(row.temp),"bam_file"],  pileupParam=PileupParam())

  proportion.temp<-data.frame(position=as.character(), fraction=as.numeric())

  if (paste0(source.temp,"_coverage") %in% colnames(all_overlaping_sites)){

    nonzero.temp<-all_overlaping_sites[,paste0(source.temp,"_coverage")]>0
    
    # For each strain-informative site (SSI variant from closest reference covered by taxon-specific metagenomic read)
    for(covered_variant_pos in all_overlaping_sites[nonzero.temp,"position"][1:min(nrow(
      all_overlaping_sites[nonzero.temp,]),pileupsites)]){

      # Compare alleles in metagenomic reads at this position to SSI VCF/reference
      ref.temp<-all_overlaping_sites$reference_allele[all_overlaping_sites$position==covered_variant_pos]
      alt.temp<-all_overlaping_sites$alternate_allele[all_overlaping_sites$position==covered_variant_pos]
      alleles.temp<-pileup.temp[pileup.temp$pos==covered_variant_pos,]
      refcount.temp<-sum(subset(alleles.temp$count, alleles.temp$nucleotide==ref.temp))
      altcount.temp<-sum(subset(alleles.temp$count, alleles.temp$nucleotide==alt.temp))
      totalcount.temp <- sum(alleles.temp$count)
      fraction.temp<-altcount.temp/totalcount.temp
      binary.temp <- ifelse(fraction.temp >0, 1,0)
      proportion.temp<-rbind(proportion.temp, data.frame(position=covered_variant_pos, fraction=binary.temp))
}
    
    assign(paste0("allele_support.",subject,".",reference_genome,".", source.temp), proportion.temp)
    
    # Append locus information to summary dataframe for sample
    allele_support.summary<-rbind(allele_support.summary, data.frame(
      subject=subject,
      sample=sample,
      source=source.temp,
      ssi_isolate=ssi_isolate,
      reference_genome=reference_genome,
      ssi_variants=nrow(subset(ssi_variants, ssi_variants$genotype.ssi==1)),
      covered_sites=length(subset(nonzero.temp, nonzero.temp==TRUE)),
      proportion=mean(proportion.temp$fraction, na.rm=T),
      manual_call=manual_call)
    )
    
    print(source.temp)
    print(mean(proportion.temp$fraction))
}
}

}

# Import an annotate microbiome-SSI pair manifest with concordance values
GCS_Pair_Table_v5 <- read_excel("/Users/drlong/Library/CloudStorage/OneDrive-UW/projects/spine_microbiome/GCS Pair Table v5.xlsx", skip = 1)

sites.function <- function(x) {
  if(identical(x, integer(0))) 0 else x
}

match.function <- function(x) {
  if(identical(x, numeric(0))) NA else x
}

allelic_support_transposed<-data.frame()
#ssi_isolate<-"A0008"
for (ssi_isolate in GCS_Pair_Table_v5$`SSI Study Isolate and SRA ID`){
    allelic_support_transposed<-rbind(allelic_support_transposed,
                                      data.frame(
                                        isolate=ssi_isolate,
                                        rectal_sites=sites.function(allele_support.summary[allele_support.summary$ssi_isolate==ssi_isolate & allele_support.summary$source=="rectal","covered_sites"]),
                                        rectal_match=match.function(allele_support.summary[allele_support.summary$ssi_isolate==ssi_isolate & allele_support.summary$source=="rectal","proportion"]),
                                        nasal_sites=sites.function(allele_support.summary[allele_support.summary$ssi_isolate==ssi_isolate & (allele_support.summary$source=="nasal"|allele_support.summary$source=="nares"),"covered_sites"]),
                                        nasal_match=match.function(allele_support.summary[allele_support.summary$ssi_isolate==ssi_isolate & (allele_support.summary$source=="nasal"|allele_support.summary$source=="nares"),"proportion"]),
                                        skin_sites=sites.function(allele_support.summary[allele_support.summary$ssi_isolate==ssi_isolate & (allele_support.summary$source=="skin" | allele_support.summary$source=="skin_cervicothoracic"),"covered_sites"]),
                                        skin_match=match.function(allele_support.summary[allele_support.summary$ssi_isolate==ssi_isolate & (allele_support.summary$source=="skin" | allele_support.summary$source=="skin_cervicothoracic"),"proportion"])
                                      ))
    
  }


if (export_reports==T){
write.csv(allelic_support_transposed, file="../reports_intermediatefiles/allelic_support_summary.csv")
}






# Visualize and analyze output -----

library(ggplot2)

library(mclust)


# Subset to cases with intermedite levels of allelic concordance
allele_support.summary.subset<-subset(allele_support.summary, allele_support.summary$proportion >0 & allele_support.summary$proportion<1)

# Generate mixture model of intermediate values
gmm_allelicconcordance = Mclust(
  log10(1-(allele_support.summary.subset$proportion))
, G=1:2)

# Annotate with mixture model calls
gmm_calls<-data.frame(value=gmm_allelicconcordance[["data"]],
           call=gmm_allelicconcordance[["classification"]],
           uncertainty=gmm_allelicconcordance[["uncertainty"]])
       
allele_support.summary[(allele_support.summary$proportion >0 
                       & allele_support.summary$proportion<1)
                       & !is.nan(allele_support.summary$proportion),"uncertainty"]<-gmm_calls$uncertainty

allele_support.summary[(allele_support.summary$proportion >0 
                        & allele_support.summary$proportion<1)
                       & !is.nan(allele_support.summary$proportion),"gmm_call"]<-gmm_calls$call


allele_support.summary[allele_support.summary$proportion==1
                       & !is.nan(allele_support.summary$proportion)
                       ,"gmm_call"]<-1

allele_support.summary[allele_support.summary$proportion==1
                       & !is.nan(allele_support.summary$proportion)
                       ,"uncertainty"]<-0

# Depth filter
allele_support.summary.subset<-subset(allele_support.summary, allele_support.summary$covered_sites>10)

# Plot mixture model classifications and uncertainty
mixture_model_plot<-ggplot(allele_support.summary.subset) + 
  geom_smooth(aes(y=uncertainty*10, x=proportion), se=F, fullrange=T, method = "loess", span=0.1, color="grey")+
  geom_histogram(aes(x=proportion, fill=as.factor(gmm_call)), bins=100)+
  scale_y_continuous(
    # Features of the first axis
    name = "Microbiome-SSI Pair Count",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*.1, name="Mixture Model Uncertainty")
  )+
  scale_x_continuous(labels = scales::percent, name="Allelic Concordance at Strain-Specific Sites")+
  scale_color_manual(values=c("#1e6091","#d9ed92"), aesthetics = c("colour","fill"), name="Classification",labels=c("Endogenous", "Exogenous"))+
  theme_minimal()

mixture_model_plot

if(export_plots==T){
ggsave("/figures/mixture_model_plot.png",
       plot=mixture_model_plot, width=10, height=5, units="in", dpi = 500, bg="white")
}

# Stage tabular data for export
allele_support.summary.subset.export <-data.frame(
  Subjet_ID=allele_support.summary.subset$subject,
  Sample_ID=allele_support.summary.subset$sample,
  Sample_source_type=allele_support.summary.subset$source,
  Allelic_concordance_proportion=allele_support.summary.subset$proportion,
  Uncertainly=allele_support.summary.subset$uncertainty,
  Classification=allele_support.summary.subset$gmm_call
)


allele_support.summary.subset.export[allele_support.summary.subset.export$Classification==1,"Classification"]<-"Endogenous"
allele_support.summary.subset.export[allele_support.summary.subset.export$Classification==2,"Classification"]<-"Exogenous"


addWorksheet(tabular_data, sheetName = "Fig S4")

writeDataTable(tabular_data, sheet = "Fig S4", 
               x = allele_support.summary.subset.export,
               colNames = T, rowNames = F)


# 
# ggsave("/Users/drlong/Library/CloudStorage/Dropbox/UW/research/projects/SSI/HMC_spine_SSI/manuscripts/microbiome/figures/mixture_model_plot.pdf", 
#        plot=, width=10, height=5, units="in", dpi = 500, bg="white")
# 
# 
# 
# ggplot(data, aes(x=day)) +
#   
#   geom_line( aes(y=temperature), size=2, color="blue") + 
#   geom_line( aes(y=price / coeff), size=2, color="red") +
#   
#   scale_y_continuous(
#     
#     # Features of the first axis
#     name = "Temperature (Celsius °)",
#     
#     # Add a second axis and specify its features
#     sec.axis = sec_axis(~.*coeff, name="Price ($)")
#   ) + 
#   
#   theme_ipsum() +
#   
#   theme(
#     axis.title.y = element_text(color = temperatureColor, size=13),
#     axis.title.y.right = element_text(color = priceColor, size=13)
#   ) +
#   
#   ggtitle("Temperature down, price up")
# 
# 
# 
# 
# 
# 
# 
# 
# ggplot(allele_support.summary, aes(y=uncertainty, x=proportion)) + 
#   theme_minimal()+
#   geom_smooth()
# 
# 
# ggplot(allele_support.summary, aes(x=proportion, fill=as.factor(gmm_call))) + 
#   geom_histogram(bins=100)+
#   theme_minimal()+
#   geom_density(aes(x=uncertainty))
# 
# 
# merge(allele_support.summary$proportion,gmm_calls$value, by.x = , by.y= , all.x=T, all.y = F)
# 
# 
# 
# # nrow(subset(allele_support.summary, allele_support.summary$covered_sites>10 & allele_support.summary$proportion>=0.8))
# # nrow(subset(allele_support.summary, allele_support.summary$covered_sites>10 & allele_support.summary$proportion<0.8))