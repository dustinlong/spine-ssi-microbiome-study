library(readr)
library(vcfR)
library(tidyverse)
library(readxl)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Rsamtools")

 library(Rsamtools)

# Import manifest
GCS_Pair_Table <- read_excel(paste0(onedrive_dir,"spine_microbiome/GCS Pair Table v2.xlsx"))

# Subset to E coli SSI cases
GCS_Pair_Table<-subset(GCS_Pair_Table, GCS_Pair_Table$`SSI Isolate Species`=="Escherichia coli")


# Generate summary dataframe
comparison.summary<-data.frame()

# For each E coli SSI
for(row.temp in as.numeric(rownames(GCS_Pair_Table))){
  # Wait and clean up environment
  Sys.sleep(10)
  rm(list=setdiff(ls(), c("comparison.summary","pileupsites","GCS_Pair_Table","row.temp", "dropbox_dir","github_dir","host", "onedrive_dir")))
  
  # Store case parameters
  subject<-as.character(GCS_Pair_Table[as.character(row.temp),"Patient ID"])
  print(subject)
  
  sample<-as.character(GCS_Pair_Table[as.character(row.temp),"Library ID"])
  print(sample)
  
  reference_genome<-as.character(GCS_Pair_Table [as.character(row.temp),"SSI Isolate Closest Reference Genome"])
  print(reference_genome)


# Identify files for preop-SSI "case" pair
organisms<-
  list.files(path=paste0(dropbox_dir,"UW/research/projects/SSI/HMC_spine_SSI/sequencing_results/genome_capture_seq/variants_alignments/subject_",subject,
                         "/sample_",sample,
                         "/"), full.names=T)

organism_dir<-
  organisms[grepl(reference_genome,organisms)]


list<-list.files(path=organism_dir, full.names=T, recursive = T)

ssi_vcf <- list[grepl(".*ssi_alignments.*vcf$",list, perl = T)]

ssi_bam <- list[grepl(paste0(".*","ssi_alignments",".*bam$"),list, perl = T)]

reference_file <- list[grepl(".*fa$",list, perl = T)]

sample_names <-
  list.dirs(path=paste0(organism_dir, "/same_patient/"), full.names=F, recursive = F)


# Identify corresponding files for other patients also with E coli SSI ("controls"), aligned to closest matching reference genome for case SSI
control_bams<-
  list.files(path=paste0(organism_dir, "/other_patients/"), pattern = "*.shortread_genome.bam$",  full.names=T, recursive = T)

control_bams<- subset(control_bams, !grepl("archive",control_bams))

gcs_samples <- data.frame(sample_name=sample_names, class="case")

gcs_samples <- subset(gcs_samples, gcs_samples$sample_name=="rectal"
                      | gcs_samples$sample_name=="postop_stool")

for (sample_name.temp in gcs_samples$sample_name){

  gcs_samples[gcs_samples$sample_name==sample_name.temp,"bam_file"]<-  list[grepl(paste0(".*/",sample_name.temp,"/.*bam$"),list, perl = T)]
  
}


gcs_controls <- data.frame(bam_file=control_bams, class="control")

gcs_controls$sample_name <- gsub(".shortread_genome.bam","",basename(control_bams))


gcs_case_control<-rbind(gcs_samples,gcs_controls)


# Set track order
gcs_case_control[gcs_case_control$class=="case","order"]<-100
gcs_case_control[gcs_case_control$class=="control","order"]<-200

gcs_case_control[gcs_case_control$sample_name=="skin","order"]<-1
gcs_case_control[gcs_case_control$sample_name=="nares","order"]<-2
gcs_case_control[gcs_case_control$sample_name=="nasal","order"]<-2
gcs_case_control[gcs_case_control$sample_name=="rectal","order"]<-3


# Import SSI VCF
ssi_import <- read.vcfR(file= ssi_vcf, verbose = FALSE)

ssi_genotypes <- as.data.frame(extract.gt(ssi_import, element = "GT"))

ssi_alleles <- as.data.frame(ssi_import@fix)

ssi_variants <- data.frame("chromosome"=ssi_alleles$CHROM, 
                           "position"=as.numeric(ssi_alleles$POS), 
                           "reference_allele"=ssi_alleles$REF,
                           "alternate_allele"=ssi_alleles$ALT,
                           "genotype.ssi"=ssi_genotypes[[1]])



# Import GenCap-Seq depth of coverage info
for (gcs_sample.row in rownames(gcs_case_control)){
  sample_name.temp <- gcs_case_control[gcs_sample.row,"sample_name"]
  
  # Generate coverage file if does not exist
  if(!(file.exists(paste0(gcs_case_control[gcs_sample.row,"bam_file"],".coverage")))){
    print("Generating coverage file")
    system(paste0("samtools depth -a ",
                  paste0(gcs_case_control[gcs_sample.row,"bam_file"])
                  ," > "
                  ,paste0(gcs_case_control[gcs_sample.row,"bam_file"],".coverage"))
    )
  }
  
  
  # Import coverage file
  coverage_file.temp <-    as.data.frame(read_delim(paste0(gcs_case_control[gcs_sample.row,"bam_file"],".coverage"), 
                                      delim = "\t", escape_double = FALSE, 
                                      col_names = FALSE, trim_ws = TRUE, col_types = cols()))
  
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


# Navigate BAM pileup

PileupParam(max_depth=250, min_base_quality=13, min_mapq=0,
            min_nucleotide_depth=1, min_minor_allele_depth=0,
            distinguish_strands=FALSE, distinguish_nucleotides=TRUE,
            ignore_query_Ns=TRUE, include_deletions=TRUE, include_insertions=TRUE,
            left_bins=NULL, query_bins=NULL, cycle_bins=NULL)

for(row.temp in as.numeric(rownames(gcs_case_control))){

  source.temp<-gcs_case_control[as.character(row.temp),"sample_name"]
  
  print(source.temp)

  
  
  pileup.temp<-pileup(gcs_case_control[as.character(row.temp),"bam_file"],  pileupParam=PileupParam())
  
  proportion.temp<-data.frame(position=as.character(), fraction=as.numeric())
  
  if (paste0(source.temp,"_coverage") %in% colnames(all_overlaping_sites)){
  
  nonzero.temp<-all_overlaping_sites[,paste0(source.temp,"_coverage")]>0
  
  # For each strain-informative site, check for variant SSI strain allele in taxon-specific metagenomic reads mapped to site
  for(covered_variant_pos in all_overlaping_sites[nonzero.temp,"position"][1:min(nrow(
    all_overlaping_sites[nonzero.temp,]),5000)]){
    
    ref.temp<-all_overlaping_sites$reference_allele[all_overlaping_sites$position==covered_variant_pos]
    alt.temp<-all_overlaping_sites$alternate_allele[all_overlaping_sites$position==covered_variant_pos]
    alleles.temp<-pileup.temp[pileup.temp$pos==covered_variant_pos,]
    refcount.temp<-sum(subset(alleles.temp$count, alleles.temp$nucleotide==ref.temp))
    altcount.temp<-sum(subset(alleles.temp$count, alleles.temp$nucleotide==alt.temp))
    totalcount.temp <- sum(alleles.temp$count)
    fraction.temp<-altcount.temp/totalcount.temp
    binary.temp <- ifelse(fraction.temp >0, 1,0)
    proportion.temp<-rbind(proportion.temp, data.frame(position=covered_variant_pos, fraction=binary.temp))
    
    assign(paste0("allele_support.",source.temp), proportion.temp)
  
  }
  
  print(mean(proportion.temp$fraction, na.rm=T))
  gcs_case_control[gcs_case_control$sample_name==source.temp,"case_subject"]<-subject
  gcs_case_control[gcs_case_control$sample_name==source.temp,"allelic_support_fraction"]<-mean(proportion.temp$fraction, na.rm=T)
    
  }
}

# Summary table of allelic concordance of reads from control patients
for(comparator.temp in subset(gcs_case_control$sample_name, gcs_case_control$class=="control")){

comparator.alleles.temp<-paste0("allele_support.",
                   comparator.temp)
  
if (exists(comparator.alleles.temp)){

comparison.summary<-rbind(comparison.summary,
                          gcs_case_control)
}
}

# Calculate relative metagenome-SSI strain similarity distances for each SSI node 
for(case_subject.temp in unique(comparison.summary$case_subject)){
  comparison.summary[comparison.summary$case_subject==case_subject.temp,"scaled_distance"] <-
  (1-comparison.summary[comparison.summary$case_subject==case_subject.temp,"allelic_support_fraction"]) /
  sum(1-comparison.summary[comparison.summary$case_subject==case_subject.temp,"allelic_support_fraction"])*100
    
}
}
