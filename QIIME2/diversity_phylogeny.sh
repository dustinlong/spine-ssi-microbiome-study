conda activate qiime2-2019.10

qiime2_resource_dir=

# Run parameters
working_dir=

n_cores=50

# SEPP
  qiime fragment-insertion sepp \
  --i-representative-sequences $working_dir/seqs-filtered.qza \
  --i-reference-database $qiime2_resource_dir/sepp-refs-fragment-phylogeny/sepp-refs-silva-128.qza \
  --p-threads $n_cores \
  --o-tree $working_dir/insertion-tree.qza \
  --o-placements $working_dir/insertion-placements.qza
  
# Alpha rarefaction analysis
qiime diversity alpha-rarefaction \
  --i-table $working_dir/table-samplefiltered.qza \
  --p-min-depth 10 \
  --p-steps 30 \
  --p-max-depth 60000 \
  --m-metadata-file $metadata_file \
  --o-visualization $working_dir/alpha_rarefaction_curves.qzv 

## Review and select value that retains majority of samples and captures the majority of OTUs
rarefaction_sampling_depth=7500

# Core diversity metrics using insertion-derived tree
  qiime diversity core-metrics-phylogenetic \
  --i-phylogeny $working_dir/insertion-tree.qza \
  --i-table $working_dir/table-samplefiltered.qza \
  --p-n-jobs $n_cores \
  --p-sampling-depth $rarefaction_sampling_depth \
  --m-metadata-file $metadata_file \
  --output-dir $working_dir/core-metrics-results-insertion