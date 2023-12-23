conda activate qiime2-2019.10

qiime2_resource_dir=

# Run parameters
working_dir=
metadata_file=
classifier_file=$qiime2_resource_dir/classifiers/classifier-347to803-var50.qza
n_cores=50

 # Run the feature classifier on filtered
qiime feature-classifier classify-sklearn \
  --i-classifier $classifier_file \
  --i-reads $working_dir/seqs-filtered.qza \
  --p-n-jobs $n_cores \
  --o-classification $working_dir/taxonomy.qza


# Generate feature classification summary data
qiime metadata tabulate \
  --m-input-file $working_dir/taxonomy.qza \
  --o-visualization $working_dir/taxonomy.qzv

 # Generate barplots for initial review
 qiime taxa barplot \
 --i-table $working_dir/table-samplefiltered.qza \
 --i-taxonomy $working_dir/taxonomy.qza \
 --m-metadata-file $metadata_file \
 --o-visualization $working_dir/taxa-bar-plots.qzv


  
  