### 16S rRNA metagenomics

#### Installation of the necessary tools
```bash
# Installation of seqkit (for correcting paired-end reads)
conda install -c bioconda seqkit

# Installation of PEAR (for pairing paired-end reads)
conda install -c bioconda pear
```

#### Pairing the paired-end sequence reads
```bash
# Prepare the paired-end sequence reads
seqkit sana FL_R1_001.fastq.gz -o FL_R1_001.cleaned.fastq 2>FL_R1_001.cleaned.fastq.log
gzip FL_R1_001.cleaned.fastq
seqkit sana FL_R2_001.fastq.gz -o FL_R2_001.cleaned.fastq 2>FL_R2_001.cleaned.fastq.log
gzip FL_R2_001.cleaned.fastq
seqkit pair -1 FL_R1_001.cleaned.fastq.gz -2 FL_R2_001.cleaned.fastq.gz 2>FL_R2_001.cleaned.paired.log

# Paired-end sequence reads
pear -f FL_R1_001.cleaned.paired.fastq.gz -r FL_R2_001.cleaned.paired.fastq.gz -o 16S_rRNA_seq_30_1302373217__FL
for i in 16S*.fastq; do gzip $i; done
```

#### Installation of QIIME2
```bash
# Installation of QIIME2
conda install -n base -c conda-forge conda=25.11.1
conda env remove -n qiime2-amplicon-2026.1 -y
conda env create -n qiime2-amplicon-2026.1 --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml

# Testing installation of QIIME2
conda activate qiime2-amplicon-2026.1
qiime info
```

#### Main QIIME2 analysis
```bash
# Create a TSV-formatted manifest file
printf "sample-id\tabsolute-filepath\n" > manifest.tsv
printf "FossilLake\t%s\n" "$(realpath 16S_rRNA_seq_30_1302373217__FL.assembled.fastq.gz)" >> manifest.tsv

# Import the FASTQ file
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path demux-single-end.qza \
  --input-format SingleEndFastqManifestPhred33V2

# Summarize read quality
qiime demux summarize \
  --i-data demux-single-end.qza \
  --o-visualization demux-summary.qzv

# Denoise single reads with DADA2
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux-single-end.qza \
  --p-trunc-len 250 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza \
  --o-base-transition-stats base-transition-stats.qza

# Conduct filtering
qiime feature-table summarize \
  --i-table table.qza \
  --o-feature-frequencies feature-frequencies.qza \
  --o-sample-frequencies sample-frequencies.qza \
  --o-summary table-summary.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

# Obtain a classifer for the common 515F/806R V4 region
wget \
  -O "gg-13-8-99-515-806-nb-classifier.qza" \
  "https://data.qiime2.org/classifiers/sklearn-1.4.2/greengenes/gg-13-8-99-515-806-nb-classifier.qza"

# Assign taxonomy to the sequences
qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

# View the taxonomy table
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# Make a taxa composition plot
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file manifest.tsv \
  --o-visualization taxa-bar-plots.qzv

# Export the taxonomy table to plain text
qiime tools export \
  --input-path taxonomy.qza \
  --output-path exported-taxonomy

qiime tools export \
  --input-path table.qza \
  --output-path exported-table
```


#### Merging features sharing the same taxonomic label into one row
```bash
# Make a readable taxonomy table
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# Combine sequences with taxonomy
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --o-visualization rep-seqs-with-taxonomy.qzv
```

#### Make a genus-level abundance table
```bash
# Make a genus-level abundance table
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table genus-table.qza

# Summarize the genus-level abundance table
qiime feature-table summarize \
  --i-table genus-table.qza \
  --o-feature-frequencies genus-feature-frequencies.qza \
  --o-sample-frequencies genus-sample-frequencies.qza \
  --o-summary genus-table-summary.qzv

# Export the genus table to a plain text file
qiime tools export \
  --input-path genus-table.qza \
  --output-path exported-genus-table
biom convert \
  -i exported-genus-table/feature-table.biom \
  -o genus-table.tsv \
  --to-tsv
```

#### Make a family-level abundance table
```bash
# Make a family-level abundance table
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table family-table.qza

# Summarize the family-level abundance table
qiime feature-table summarize \
  --i-table family-table.qza \
  --o-feature-frequencies family-feature-frequencies.qza \
  --o-sample-frequencies family-sample-frequencies.qza \
  --o-summary family-table-summary.qzv

# Export the family-level table to a plain text file
qiime tools export \
  --input-path family-table.qza \
  --output-path exported-family-table

biom convert \
  -i exported-family-table/feature-table.biom \
  -o family-table.tsv \
  --to-tsv
```