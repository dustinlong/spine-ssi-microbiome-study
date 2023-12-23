conda activate qiime2-2019.10

qiime2_resource_dir=

working_dir=

mkdir $working_dir

cd $working_dir

# Merge DADA2-processed reads
## Merge tables
qiime feature-table merge \
 --i-tables $working_dir/*
 --o-merged-table $working_dir/table-merged.qza
 
 ## Merge representative sequences
 qiime feature-table merge-seqs \
 --i-data $working_dir/*
 --o-merged-data $working_dir/rep-seqs-merged.qza
 
 