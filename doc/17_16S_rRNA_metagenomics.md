## 16S rRNA metagenomics

### Installation of the necessary tools
```bash
# Installation of seqkit (for correcting paired-end reads)
conda install -c bioconda seqkit

# Installation of PEAR (for pairing paired-end reads)
conda install -c bioconda pear
```

### Pairing the paired-end sequence reads
```bash

SAMPLE=30_1297791016_H1

# Prepare the paired-end sequence reads
seqkit sana ${SAMPLE}_R1_001.fastq.gz -o ${SAMPLE}_R1_001.cleaned.fastq 2>${SAMPLE}_R1_001.cleaned.fastq.log
gzip ${SAMPLE}_R1_001.cleaned.fastq
seqkit sana ${SAMPLE}_R2_001.fastq.gz -o ${SAMPLE}_R2_001.cleaned.fastq 2>${SAMPLE}_R2_001.cleaned.fastq.log
gzip ${SAMPLE}_R2_001.cleaned.fastq
seqkit pair -1 ${SAMPLE}_R1_001.cleaned.fastq.gz -2 ${SAMPLE}_R2_001.cleaned.fastq.gz 2>${SAMPLE}_R2_001.cleaned.paired.log

# Paired-end sequence reads
pear -f ${SAMPLE}_R1_001.cleaned.paired.fastq.gz -r ${SAMPLE}_R2_001.cleaned.paired.fastq.gz -o 16S_rRNA_seq_${SAMPLE}
for i in 16S*_${SAMPLE}*.fastq; do gzip $i; done
```

### Installation of QIIME2
```bash
# Installation of QIIME2
conda install -n base -c conda-forge conda=25.11.1
conda env remove -n qiime2-amplicon-2026.1 -y
conda env create -n qiime2-amplicon-2026.1 --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml

# Testing installation of QIIME2
conda activate qiime2-amplicon-2026.1
qiime info
```

### Main QIIME2 analysis

#### Preparatory steps
```bash
# Define variables
LOCATION="HorsethiefReservoir"   # CHANGE AS NEEDED!
SAMPLE="30_1302373217_HT"        # CHANGE AS NEEDED!

# Make sample-specific folder and operate only in the folder
mkdir -p 16S_rRNA_seq_${SAMPLE}
mv 16S_rRNA_seq_${SAMPLE}.assembled.fastq.gz 16S_rRNA_seq_${SAMPLE}
cd 16S_rRNA_seq_${SAMPLE}

# Create a TSV-formatted manifest file
printf "sample-id\tabsolute-filepath\n" > manifest.tsv
printf "%s\t%s\n" "$LOCATION" "$(realpath 16S_rRNA_seq_${SAMPLE}.assembled.fastq.gz)" >> manifest.tsv
```

#### Importing FASTQ, assessing read quality
```bash
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

qiime tools view demux-summary.qzv
```

#### Denoising reads
```bash
# Denoise single reads with DADA2
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux-single-end.qza \
  --p-trunc-len 250 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza \
  --o-base-transition-stats base-transition-stats.qza

qiime tools view denoising-stats.qzv
```

#### Summarize feature table, sequences, and denoising statistics
```bash
qiime feature-table summarize \
  --i-table table.qza \
  --o-feature-frequencies feature-frequencies.qza \
  --o-sample-frequencies sample-frequencies.qza \
  --o-summary table-summary.qzv

qiime tools view table-summary.qzv
```

```bash
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

qiime tools view rep-seqs.qzv
```

#### Assign taxonomy based on sequence classification
```bash
# Obtain a classifer for the common 515F/806R V4 region
wget \
  -O "gg-13-8-99-515-806-nb-classifier.qza" \
  "https://data.qiime2.org/classifiers/sklearn-1.4.2/greengenes/gg-13-8-99-515-806-nb-classifier.qza"

# Assign taxonomy to the sequences
qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
```

#### Generate plain text taxonomy table
```bash
# View the taxonomy table
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime tools view taxonomy.qzv
```

```bash
# Make a taxa composition plot
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file manifest.tsv \
  --o-visualization taxa-bar-plots.qzv

qiime tools view taxa-bar-plots.qzv
```

```bash
# Export the taxonomy table to plain text
qiime tools export \
  --input-path taxonomy.qza \
  --output-path exported-taxonomy

qiime tools export \
  --input-path table.qza \
  --output-path exported-table
```

#### Merging sequences sharing taxonomic label into rows
```bash
# Make a readable taxonomy table
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime tools view taxonomy.qzv

# Combine sequences with taxonomy
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --o-visualization rep-seqs-with-taxonomy.qzv
  
qiime tools view rep-seqs-with-taxonomy.qzv
```

[!NOTE]
The following visualization is helpful.
```bash
qiime tools view rep-seqs-with-taxonomy.qzv
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

qiime tools view genus-table-summary.qzv
```

```bash
# Export the genus table to a plain text file
qiime tools export \
  --input-path genus-table.qza \
  --output-path exported-genus-table
biom convert \
  -i exported-genus-table/feature-table.biom \
  -o genus-table.tsv \
  --to-tsv
```

#### Generate visualization of taxa in sample
```python
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------
# Read sample name from manifest.tsv
# ---------------------------
manifest = pd.read_csv("manifest.tsv", sep="\t")
sample_name = manifest.loc[0, "sample-id"]

# ---------------------------
# Read genus table
# Skip the first comment line: "# Constructed from biom file"
# ---------------------------
df = pd.read_csv("genus-table.tsv", sep="\t", skiprows=1)

# Clean column names just in case
df.columns = df.columns.str.strip()

# Rename the taxonomy column to something easier to work with
df = df.rename(columns={"#OTU ID": "Taxon"})

# The abundance column should match the sample name from manifest.tsv
if sample_name not in df.columns:
    raise ValueError(
        f"Sample '{sample_name}' was not found in genus-table.tsv. "
        f"Available columns: {list(df.columns)}"
    )

# ---------------------------
# Extract genus from taxonomy string
# ---------------------------
def extract_genus(taxon):
    taxon = str(taxon).strip()

    # Handle fully unassigned rows
    if taxon.lower().startswith("unassigned"):
        return "Unassigned"

    parts = [p.strip() for p in taxon.split(";")]

    for part in parts:
        if part.startswith("g__"):
            genus = part.replace("g__", "").strip()
            return genus if genus else "Unassigned"

    return "Unassigned"

df["Genus"] = df["Taxon"].apply(extract_genus)

# ---------------------------
# Aggregate counts at genus level
# ---------------------------
genus_counts = df.groupby("Genus")[sample_name].sum()
genus_counts = genus_counts[genus_counts > 0]

# Convert to relative abundance
rel = genus_counts / genus_counts.sum()

# Group small taxa into "Other"
threshold = 0.03  # 3%
major = rel[rel >= threshold].copy()
minor = rel[rel < threshold].sum()

if minor > 0:
    major.loc["Other"] = minor

major = major.sort_values(ascending=False)

# ---------------------------
# Plot donut chart
# ---------------------------
fig, ax = plt.subplots(figsize=(8, 8))

wedges, texts, autotexts = ax.pie(
    major,
    labels=None,
    autopct=lambda p: f"{p:.1f}%" if p >= 3 else "",
    startangle=90,
    wedgeprops={"width": 0.4, "edgecolor": "white", "linewidth": 1}
)

# Donut hole
centre_circle = plt.Circle((0, 0), 0.60, fc="white")
ax.add_artist(centre_circle)

# Title in center
ax.text(0, 0, sample_name, ha="center", va="center", fontsize=12)

# Legend
ax.legend(
    wedges,
    major.index,
    title="Genus",
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    frameon=False
)

ax.set_title(f"Genus composition of {sample_name}")
plt.tight_layout()
plt.show()
```
