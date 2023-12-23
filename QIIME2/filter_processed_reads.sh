conda activate qiime2-2019.10

qiime2_resource_dir=
working_dir=
metadata_file=

mkdir $working_dir

cd $working_dir

# Check metadata formatting
qiime metadata tabulate \
  --m-input-file $metadata_file \
  --o-visualization $working_dir/tabulated-sample-metadata.qzv
  
# Sequence-level filters:
# Filter out sequences based on length
	#--p-global-max 550 \
qiime rescript filter-seqs-length \
	--i-sequences $working_dir/rep-seqs-merged.qza \
	--p-global-min 406 \
	--o-filtered-seqs $working_dir/seqs-filtered.qza \
	--o-discarded-seqs $working_dir/seqs-discarded.qza

# Filter feature table to just sequences that passed sequence-level filters
qiime feature-table filter-features \
  --i-table $working_dir/table-merged.qza \
  --m-metadata-file $working_dir/seqs-filtered.qza \
  --o-filtered-table $working_dir/table-seqfiltered.qza

# Sample-level filters:
# Filter feature table to just samples in experiment/study
qiime feature-table filter-samples \
  --i-table $working_dir/table-seqfiltered.qza \
  --m-metadata-file $metadata_file \
  --p-where "[Study]='Long_Spine_SSI'" \
  --o-filtered-table $working_dir/table-samplefiltered.qza
  

conda deactivate