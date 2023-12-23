conda activate qiime2-2019.10

qiime2_resource_dir=

# Train classifier
  qiime feature-classifier extract-reads \
  --i-sequences $qiime2_resource_dir/silva-138-99-seqs.qza \
  --p-f-primer GGAGGCAGCAGTRRGGAAT \
  --p-r-primer CTACCRGGGTATCTAATCC \
  --p-min-length 406 \
  --p-max-length 506 \
  --p-n-jobs 90 \
  --o-reads $qiime2_resource_dir/classifiers/ref-seqs-347to803-var50.qza

# Train the classifier using these sequences
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $qiime2_resource_dir/classifiers/ref-seqs-347to803-var50.qza \
  --i-reference-taxonomy $qiime2_resource_dir/silva-138-99-tax.qza \
  --o-classifier $qiime2_resource_dir/classifiers/classifier-347to803-var50.qza

 # conda deactivate
  