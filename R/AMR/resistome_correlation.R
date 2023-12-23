library(readr)
library(vcfR)
library(tidyverse)
library(readxl)

write_figures<-F # Overwrite figures
sccmec_stringent<-F # optional sensitivity analysis for analysis of sccmec IV hypervariable region rather than mecA gene

# Import manifest of preop microbiome - postop SSI pairs
GCS_Sample_Table <- read_excel(paste0(share_dir, "spine_microbiome/GCS Sample Table.xlsx"))

# Import list of implicated AMR genes in cases of prophylaxis resistant infection
isolate_amr_genes <- read_excel(paste0(share_dir,"spine_microbiome/isolate_amr_genes.xlsx"))
isolate_amr_genes <- subset(isolate_amr_genes, isolate_amr_genes$accessory_gene_implicated==1)

GCS_Sample_Table <- subset(GCS_Sample_Table, GCS_Sample_Table$`Patient ID` %in% isolate_amr_genes$patient)

# Import table of total read counts
read_counts <- read_excel("./read_counts/221213_all_samples_readcount_after_trim.xlsx", 
                          col_names = c("readcount","sampleprefix"))


# Generate list of AMR genes based on filenames
amr_genes<-
  list.dirs(path=paste0(share_dir,"/sequencing_results/genome_capture_seq/blast_20_genes"), full.names=F)

amr_genes<-subset(amr_genes, amr_genes!="")

# Subset to genes most directly implicated in SSI isolates and not regulated by de novo mutation
amr_genes <- subset(amr_genes,  !(amr_genes %in% c("ampC",
                                                   "icaC",
                                                   "blaEC",
                                                   "blaEC_b",
                                                   "KpnH",
                                                   "vanR",
                                                   "vanS",
                                                   "vanT",
                                                   "blaACT",
                                                   "SCCmec_IV" # from sensitivity test for SCCmec type IV hyper-variable region
                                                   )))
# Import list of AMR gene lengths
resistome_gene_lengths <- read_excel("./resistome_gene_lengths.xlsx")

amr_correlation.summary <- data.frame()
evalue.holder <- data.frame()

# For each AMR gene
for (amr_gene.temp in amr_genes){
  
  # Store associated file list
  amr_results.temp <-
    list.files(path=
                 paste0(share_dir,"/blast_20_genes/",amr_gene.temp)
               , full.names=T)
  
  # For each preop resistome sample
  for(fastqprefix.temp in unique(GCS_Sample_Table$FASTQ_prefix_AMR)){
    print(fastqprefix.temp)
    
    # Store characteristics
    subject.temp <- as.character(unique(GCS_Sample_Table[GCS_Sample_Table$FASTQ_prefix_AMR==fastqprefix.temp,"Patient ID"]))
    source.temp <- as.character(unique(GCS_Sample_Table[GCS_Sample_Table$FASTQ_prefix_AMR==fastqprefix.temp,"Sample Source"]))
    
    # Identify AMR blast result files by gene name and prefix
    blasthits.R1.temp <- 
      amr_results.temp[grepl(paste0("/",fastqprefix.temp,".*R1"),amr_results.temp)]
    
    blasthits.R2.temp <- 
      amr_results.temp[grepl(paste0("/",fastqprefix.temp,".*R2"),amr_results.temp)]
    
    # Import AMR blast result files
    R1blastresults.temp <- read_delim(blasthits.R1.temp, 
                                      delim = "\t", escape_double = FALSE, 
                                      col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"), 
                                      trim_ws = TRUE,
                                      show_col_types = F)  
    
    R2blastresults.temp <- read_delim(blasthits.R2.temp, 
                                      delim = "\t", escape_double = FALSE, 
                                      col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"), 
                                      trim_ws = TRUE,
                                      show_col_types = F)  
    
    combined.blastresults.temp<- rbind(R1blastresults.temp, R2blastresults.temp)
    
    # Option to use more stringent criteria for SCCmec (limit to hyper-variable gene region, use more restrictive E-value threshold, require longer length)
    if(sccmec_stringent==T){
    if(amr_gene.temp=="SCCmec_IV"){
      combined.blastresults.temp<-subset(combined.blastresults.temp, 
                                         combined.blastresults.temp$sstart >19698
                                & combined.blastresults.temp$send <20272
                                & combined.blastresults.temp$length > 140
                                & combined.blastresults.temp$evalue < 10^-40)
    }
    }
    
    # Read counts
    readcount.temp <- as.numeric(read_counts[read_counts$sampleprefix==fastqprefix.temp,"readcount"])
    print(readcount.temp)
    
    # Append to summary dataframe
    if(nrow(combined.blastresults.temp)>0){
      evalue.holder<-append(evalue.holder, combined.blastresults.temp$evalue)
    }
    
    # Post hoc conversion to RPKM using GMM threshold derived at later step
    passing_hits.temp <-  nrow(subset(combined.blastresults.temp, combined.blastresults.temp$evalue < 10^-26)) /
      ( as.numeric(resistome_gene_lengths[resistome_gene_lengths$gene==amr_gene.temp,"length"])/1000 * readcount.temp/1000000 )

    # Check whether this particular AMR gene implicated in SSI for corresponding patient
    amrgene_in_SSI.temp <- as.logical(isolate_amr_genes[isolate_amr_genes$patient==subject.temp, amr_gene.temp])
    
    # Store in dataframe
    amr_correlation.summary<-rbind(amr_correlation.summary, data.frame(subject=subject.temp, 
                                                                       sample=source.temp,
                                                                       gene=amr_gene.temp,
                                                                       hits=passing_hits.temp,
                                                                       gene_in_SSI=amrgene_in_SSI.temp))
  } # End for preoperative resistome FASTQ
} # End for AMR gene


# Correct 1090 sample labels
amr_correlation.summary[amr_correlation.summary$subject=="100" & amr_correlation.summary$sample=="nasal","sample"]<-"rectal_true"
amr_correlation.summary[amr_correlation.summary$subject=="100" & amr_correlation.summary$sample=="rectal","sample"]<-"nasal_true"

amr_correlation.summary[amr_correlation.summary$subject=="100" & amr_correlation.summary$sample=="rectal_true","sample"]<-"rectal"
amr_correlation.summary[amr_correlation.summary$subject=="100" & amr_correlation.summary$sample=="nasal_true","sample"]<-"nasal"

# Transform E values
evalue.holder2<- unlist(evalue.holder)
evalue.holder3<- -log10(evalue.holder2)

# Mixture modeling for E values to derive high-specificity threshold
library(mclust)

gmm_fit_blast = Mclust(
  evalue.holder3
  , G=1:2)

# Store model assignments
gmm_calls<-data.frame(value=gmm_fit_blast[["data"]],
                      call=gmm_fit_blast[["classification"]],
                      uncertainty=gmm_fit_blast[["uncertainty"]])



evalue_gmm<-data.frame(evaluex=evalue.holder3,
                       evaluegmm=gmm_calls$value,
                       call=gmm_calls$call,
                       uncertainty=gmm_calls$uncertainty)


# Review and set high/low specificity BLAST hits based on gmm calls for E value threshold derived above
evalue_gmm$specificity <- "high"
evalue_gmm[evalue_gmm$evaluex <26, "specificity"] <- "low"

evalue_gmm$specificity <- factor(evalue_gmm$specificity, levels=c("low","high"))

# Generate histogram showing threshold for low/high-specificity BLAST hits to AMR genes
mixture_model_plot<-
  ggplot(evalue_gmm) + 
  geom_histogram(aes(x=evaluex, fill=as.factor(specificity)), bins=100)+
  scale_y_continuous(
    name = "Read Hit Count",
  )+
  scale_x_continuous(name="AMR Gene BLAST E-value")+
  scale_color_manual(values=c("#e63946","#0096c7"), aesthetics = c("colour","fill"), name="Classification",labels=c("Low-Specificity Hits", "High-Specificity Hits"))+
  theme_minimal()


mixture_model_plot


if (write_figures==T){
ggsave(
  paste0(share_dir,"/figures/resistome_blast_gmm.pdf")
  , plot=mixture_model_plot, 
  width=12,
  height = 6,
  units = "in",
  dpi = 600)
}

# Sample 1 million reads for export in tabular format
evalue_gmm.reduced<-evalue_gmm[sample(rownames(evalue_gmm), 1000000),]
evalue_gmm.reduced.export<- data.frame(
  Evalue_log10=evalue_gmm.reduced$evaluex,
  Specificity_classification=evalue_gmm.reduced$specificity
)

addWorksheet(tabular_data, sheetName = "Fig 6A")

writeDataTable(tabular_data, sheet = "Fig 6A", 
               x = evalue_gmm.reduced.export,
               colNames = T, rowNames = F)




# Prepare data for resistome-SSI grid plot
# Derive label
amr_correlation.summary$label<- paste0(amr_correlation.summary$subject, " ", amr_correlation.summary$sample)

# -log10 transformation
amr_correlation.summary$hits_transformed <-log10(amr_correlation.summary$hits+10)-1

# set plotting shape
amr_correlation.summary[amr_correlation.summary$hits == 0, "shape"] <- as.character(1)
amr_correlation.summary[amr_correlation.summary$hits > 0, "shape"] <- as.character(3)

# Generate resistome-SSI grid plot
amr_correlation_matrix<-
  ggplot(amr_correlation.summary, aes(fct_relevel(gene,
                                                  #'SCCmec_IV', # Sensitivity analysis for SCCmec type IV hyper-variable region
                                                  'mecA',
                                                  'lukF-PV','lukS-PV','qacA','qacR','vanC','blaTEM-1','blaACC-1a','blaCKO','blaPDC-5','blaCMY-2', 'blaSHV-187','oqxB32','oqxA'), fct_relevel(sample,c('rectal','skin','nasal')), color=gene_in_SSI, size = hits_transformed )) + 
  geom_point(aes(shape=shape))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_shape_manual(values = c(1,16))+
  facet_grid(fct_relevel(subject,'44','38','91','65','26','141','137','110','100') ~ ., scales = "free", space = "free", switch="y") +
  theme(strip.text.y = element_text(angle = 0))+
  scale_color_manual(values=c("#e63946","#0096c7"))+
  scale_size_continuous(breaks=c(log10(1),log10(2),log10(5),log10(10),log10(50),log10(100),log10(500),log10(1000)), labels=c("1","2","5","10","50","100","500","1,000"))+
  xlab("SSI Antimicrobial Resistance Gene") + 
  ylab("Preoperative Microbiome Sample")+
  theme(legend.position="none")

amr_correlation_matrix

if (write_figures==T){
ggsave(
  paste0(share_dir,"/figures/resistome_correlation_matrix_tag.pdf")
  , plot=amr_correlation_matrix, 
  width=6,
  height = 9,
  units = "in",
  dpi = 600)
}

# Stage data for export in tabular format
amr_correlation.summary.export<-data.frame(
  Subject=amr_correlation.summary$subject,
  Sample_source=amr_correlation.summary$sample,
  AMR_gene=amr_correlation.summary$gene,
  Resistome_BLAST_hit_count_transformed=amr_correlation.summary$hits_transformed,
  Gene_in_SSI=amr_correlation.summary$gene_in_SSI
)

addWorksheet(tabular_data, sheetName = "Fig 6B")

writeDataTable(tabular_data, sheet = "Fig 6B", 
               x = amr_correlation.summary.export,
               colNames = T, rowNames = F)



# Store hit counts for resistance-gene concordant vs discordant pairs
AMRpresent.vals <- data.frame(values=subset(amr_correlation.summary$hits_transformed, amr_correlation.summary$gene_in_SSI==TRUE))
AMRabsent.vals <- data.frame(values=subset(amr_correlation.summary$hits_transformed, amr_correlation.summary$gene_in_SSI==FALSE))

# Resistome AMR gene dose - SSI genotype comparison in non-zero cases
amr_correlation.summary.nonzero<-subset(amr_correlation.summary, amr_correlation.summary$hits_transformed>0)

null_proportion<-nrow(subset(amr_correlation.summary.nonzero, amr_correlation.summary.nonzero$gene_in_SSI==TRUE))/
  nrow(amr_correlation.summary.nonzero)

# Density plot of AMR read count and SSI genotype concordance
resistance_distribution_density<-
  ggplot(amr_correlation.summary.nonzero, aes(x=hits, fill=gene_in_SSI)) + 
  geom_density(alpha=1, position="fill", adjust=1, color="white")+
  geom_hline(yintercept=null_proportion, linetype=2, color="white", size=1)+
  scale_fill_manual(values=c("#e63946","#0096c7"))+
  scale_x_log10(breaks=c(0.01,0.1,1,10,100,1000))+
  scale_y_continuous(labels = scales::percent)+
  theme_minimal()+
  xlab("Antimicrobial Resistance Gene BLAST Hits (RPKM)") + 
  ylab("Proportion of Sample-AMR Gene Pairs")+
  theme(legend.position = "none")

if (write_figures==T){
ggsave(
  paste0(share_dir,"/figures/resistance_distribution_density.pdf")
  , plot=resistance_distribution_density, 
  width=8,
  height = 4,
  units = "in",
  dpi = 600)
}

# Stage for export in tabular format
amr_correlation.summary.nonzero.export<-data.frame(
  Resistome_BLAST_hit_count=amr_correlation.summary.nonzero$hits,
  Gene_present_in_SSI=amr_correlation.summary.nonzero$gene_in_SSI
)

addWorksheet(tabular_data, sheetName = "Fig 6D")

writeDataTable(tabular_data, sheet = "Fig 6D", 
               x = amr_correlation.summary.nonzero.export,
               colNames = T, rowNames = F)


# Option for subsetting (not used)
amr_correlation.summary.subset <- amr_correlation.summary

# Create summary dataframe to simulate combing reads from all sources within individidual patients
combined.amr.reads <- data.frame()

# Simulate combining AMR read information from all samples per subject
for (subject.temp in unique(amr_correlation.summary$subject)){
  for(gene.temp in unique(amr_correlation.summary$gene)){
    
    gene_in_SSI.temp<-unique(amr_correlation.summary[amr_correlation.summary$subject==subject.temp &
                                                       amr_correlation.summary$gene==gene.temp,
                                                     "gene_in_SSI"])
    
    combined_count.temp<-sum(amr_correlation.summary[amr_correlation.summary$subject==subject.temp &
                                                       amr_correlation.summary$gene==gene.temp,
                                                     "hits_transformed"])
    
    combined.amr.reads<-rbind(combined.amr.reads,
                              data.frame(subject=subject.temp,
                                         gene=gene.temp,
                                         gene_in_SSI=gene_in_SSI.temp,
                                         combined_count=combined_count.temp))
    
  }
}


# For SSI associated AMR genes, what was the sensitivity for detection of each swab type?
amr_correlation.summary.subset <- subset(amr_correlation.summary, amr_correlation.summary$gene_in_SSI==TRUE)

nrow(subset(amr_correlation.summary.subset, amr_correlation.summary.subset$sample=="skin" & amr_correlation.summary.subset$hits >0))/
  nrow(subset(amr_correlation.summary.subset, amr_correlation.summary.subset$sample=="skin"))

nrow(subset(amr_correlation.summary.subset, amr_correlation.summary.subset$sample=="rectal" & amr_correlation.summary.subset$hits >0))/
  nrow(subset(amr_correlation.summary.subset, amr_correlation.summary.subset$sample=="rectal"))

nrow(subset(amr_correlation.summary.subset, amr_correlation.summary.subset$sample=="nasal" & amr_correlation.summary.subset$hits >0))/
  nrow(subset(amr_correlation.summary.subset, amr_correlation.summary.subset$sample=="nasal"))


# Combined sensitivity
nrow(subset(combined.amr.reads, combined.amr.reads$gene_in_SSI==TRUE & combined.amr.reads$combined_count >0))/
  nrow(subset(combined.amr.reads, combined.amr.reads$gene_in_SSI==TRUE))

# Observed difference in preoperative read counts for AMR genes that were vs were not implicated in postoperative SSI resistance for same patient by WGS
wilcox_all.temp<-wilcox.test(subset(amr_correlation.summary$hits_transformed, amr_correlation.summary$gene_in_SSI==TRUE),
                             subset(amr_correlation.summary$hits_transformed, amr_correlation.summary$gene_in_SSI==FALSE))


# Compare this values with simulate random shuffling of preop microbiome samples between patients to derive bootstrapped distribution of concordance with SSI genotype
bootstrap_holder<-data.frame()

subject.sample.list<-unique(amr_correlation.summary$subject)


for (i in 1:10^5.5){ #Set number of iterations
  # Make shuffled sample cross-table
  subject.shufflekey.temp<-data.frame(observed=subject.sample.list,
                                      shuffled=sample(subject.sample.list))
  
  # If all subjects are shuffled
  if(length(which(subject.shufflekey.temp$observed==subject.shufflekey.temp$shuffled))==0){
    
    # Merge with random crosstable
    amr_correlation.summary.shuffled<-merge(amr_correlation.summary, subject.shufflekey.temp, by.x="subject", by.y="observed")
    
    for (subject.temp in subject.sample.list){
      # Get the shuffled subject
      shuffle_subj.temp<-unique(amr_correlation.summary.shuffled[amr_correlation.summary.shuffled$subject==subject.temp,"shuffled"])
      
      for(gene.temp in unique(amr_correlation.summary.shuffled$gene)){
        # Store whether shuffle subject resistome has gene...
        shuffle_subj_genepresent.temp<- unique(amr_correlation.summary.shuffled[amr_correlation.summary.shuffled$subject==shuffle_subj.temp
                                                                                & amr_correlation.summary.shuffled$gene==gene.temp
                                                                                ,"gene_in_SSI"]) 
        
        #...and replace that value in the rows for the true subject
        amr_correlation.summary.shuffled[amr_correlation.summary.shuffled$shuffled==subject.temp
                                         & amr_correlation.summary.shuffled$gene==gene.temp
                                         ,"gene_in_SSI_shuffled"]<-shuffle_subj_genepresent.temp
        
      }
      
    }
    
    # Perform corresponding Wilcoxon test for reads count by SSI AMR genotype concordance for shuffled dataset
    wilcox.temp<-wilcox.test(subset(amr_correlation.summary.shuffled$hits_transformed, amr_correlation.summary.shuffled$gene_in_SSI_shuffled==TRUE),
                             subset(amr_correlation.summary.shuffled$hits_transformed, amr_correlation.summary.shuffled$gene_in_SSI_shuffled==FALSE))
    
    # Store W statistic for iteration
    bootstrap_holder[i,"measure"]<-as.numeric(wilcox.temp$statistic)
  } #end if all sample shuffled
  
}


# Store dataframe copy
bootstrap_holder<-subset(bootstrap_holder, !is.na(bootstrap_holder$measure))
bootstrap_holder.backup<-bootstrap_holder

# Sample 100,000 iterations from pool
bootstrap_holder<-data.frame(measure=sample(bootstrap_holder$measure, 10^5))

# Estimate bootstrapped P value, 95% CIs
pvalue_bootstrapped<-as.numeric(2*pnorm(q=(wilcox_all.temp$statistic-mean(bootstrap_holder$measure))/sd(bootstrap_holder$measure), lower.tail=FALSE))
CI95H<-quantile(bootstrap_holder$measure,(97.5/100), na.rm = T)
CI95L<-quantile(bootstrap_holder$measure,(2.5/100), na.rm = T)

ggplot(bootstrap_holder, aes(x=measure)) +
  geom_histogram(fill="#0096c7", color="black", bins=100)+
  scale_x_continuous(limits = c(5000,16000), breaks = c(3000,5000,7000,9000,11000,13000,15000))+
  geom_vline(aes(xintercept=CI95H), color="red", linetype="dashed")+
  geom_vline(aes(xintercept=CI95L), color="red", linetype="dashed")+
  geom_vline(aes(xintercept=wilcox_all.temp$statistic), color="blue",linetype="dashed")+
  labs(title="", #Observed versus Expected Distribution of Resistome-SSI Correlation
       x="Wilcoxon signed-rank test statistic (W)", 
       y = "Bootstrap Iteration Count")+
  theme_minimal()

if (write_figures==T){
ggsave("./figures/resistome_correlation_distribution.pdf", height=4,width = 10)

}

# Export verion without CI95 markings
ggplot(bootstrap_holder, aes(x=measure)) +
  geom_histogram(fill="#e63946", color="black", bins=100)+
  scale_x_continuous(limits = c(5000,16000), breaks = c(3000,5000,7000,9000,11000,13000,15000))+
  labs(title="",
       x="Wilcoxon signed-rank test statistic (W)", 
       y = "Bootstrap Iteration Count")+
  theme_minimal()

if (write_figures==T){
ggsave("./figures/resistome_correlation_distribution_nomarks.pdf", height=4,width = 10)

}

# Stage data for export in tabular format
bootstrap_holder.export <- data.frame(
  iteration=rownames(bootstrap_holder),
  W_statistic=bootstrap_holder$measure
)


addWorksheet(tabular_data, sheetName = "Fig 6C")

writeDataTable(tabular_data, sheet = "Fig 6C", 
               x = bootstrap_holder.export,
               colNames = T, rowNames = F)


