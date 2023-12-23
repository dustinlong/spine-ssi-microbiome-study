conda activate qiime2-2019.10

qiime2_resource_dir=/resources/qiime2

# Run parameters
working_dir=
fastq_dir=
fastq_pattern=*_L001_R*_001.fastq.gz
classifier_file=$qiime2_resource_dir/classifiers/classifier-347to803-var50.qza
n_cores=

# Check pattern
ls $fastq_dir/$fastq_pattern

# Make subdir for just fastqs of interest
mkdir $fastq_dir/select_fastqs

# Copy to subdir 
cp $fastq_dir/$fastq_pattern $fastq_dir/select_fastqs

# Remove undetermined reads
rm $fastq_dir/select_fastqs/Undetermined_*

cd $working_dir

# Import reads
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $fastq_dir/select_fastqs \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path $working_dir/demux-paired-end.qza

# Check summary of imported, demuxed reads (takes about 5 min)
  qiime demux summarize \
  --i-data $working_dir/demux-paired-end.qza \
  --o-visualization $working_dir/demux.qzv

# Run DADA2
 qiime dada2 denoise-paired \
 --i-demultiplexed-seqs $working_dir/demux-paired-end.qza \
  --p-trunc-len-f 265 \
  --p-trunc-len-r 250 \
  --p-n-threads $n_cores \
  --o-representative-sequences $working_dir/rep-seqs-dada2.qza \
  --o-table $working_dir/table-dada2.qza \
  --o-denoising-stats $working_dir/stats-dada2.qza
  
 #Generate DADA2 stats 
  qiime metadata tabulate \
  --m-input-file $working_dir/stats-dada2.qza \
  --o-visualization $working_dir/stats-dada2.qzv
  

 # conda deactivate
  