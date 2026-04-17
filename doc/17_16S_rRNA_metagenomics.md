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
LOCATION="FossilLake"   # CHANGE AS NEEDED!
SAMPLE="30_1302373217_FL"        # CHANGE AS NEEDED!

# Make sample-specific folder and operate only in the folder
mkdir -p 16S_rRNA_seq_${SAMPLE}
mv 16S_rRNA_seq_${SAMPLE}.assembled.fastq.gz 16S_rRNA_seq_${SAMPLE}
cd 16S_rRNA_seq_${SAMPLE}

# Create a TSV-formatted manifest file
printf "sample-id\tabsolute-filepath\n" > manifest.tsv
printf "%s\t%s\n" "$LOCATION" "$(realpath 16S_rRNA_seq_${SAMPLE}.assembled.fastq.gz)" >> manifest.tsv
```

#### Importing FASTQ, assessing read quality, denoising reads
```bash
# Import FASTQ data into a QIIME 2 artifact
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path demux-single-end.qza \
  --input-format SingleEndFastqManifestPhred33V2

# Summarize demultiplexed sequence quality
qiime demux summarize \
  --i-data demux-single-end.qza \
  --o-visualization demux-summary.qzv

# Denoise sequences using DADA2 (single-end)
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux-single-end.qza \
  --p-trunc-len 250 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza \
  --o-base-transition-stats base-transition-stats.qza
```

#### Summarize feature table
```bash
# Summarize feature table and generate frequency outputs
qiime feature-table summarize \
  --i-table table.qza \
  --o-feature-frequencies feature-frequencies.qza \
  --o-sample-frequencies sample-frequencies.qza \
  --o-summary table-summary.qzv

# Tabulate representative sequences into a visualization
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

# Tabulate denoising statistics into a visualization
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv
```

#### Assign taxonomy based on sequence classification

##### OPTION 1: Classification using the Silva 138 database
```bash
# Download pre-trained SILVA classifier for your compatible version/region
wget \
-O "silva-138-99-515-806-nb-classifier.qza" \
"https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza"

# Assign taxonomy using the pre-trained classifier
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
```

##### Classification using BLAST
Relies on (and thus requires downloading of):
- Silva 138 SSURef NR99 515F/806R region sequences
- Silva 138 SSURef NR99 515F/806R region taxonomy
Both available from https://docs.qiime2.org/2024.10/data-resources/

```bash
# Assign taxonomy using BLAST consensus against SILVA reference database
qiime feature-classifier classify-consensus-blast \
  --i-query rep-seqs.qza \
  --i-reference-reads silva-138-99-seqs-515-806.qza \
  --i-reference-taxonomy silva-138-99-tax-515-806.qza \
  --p-perc-identity 0.97 \
  --p-query-cov 0.8 \
  --o-classification taxonomy.qza
```

##### OPTION 3: Classification using the Greengene database
```bash
# Obtain a classifer for the common 515F/806R V4 region
wget \
  -O "gg-13-8-99-515-806-nb-classifier.qza" \
  "https://data.qiime2.org/classifiers/sklearn-1.4.2/greengenes/gg-13-8-99-515-806-nb-classifier.qza"

# Assign taxonomy using the pre-trained classifier
qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
```

#### Generate taxonomy table
```bash
# Tabulate taxonomy into a visualization artifact
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# Generate taxa composition bar plot
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file manifest.tsv \
  --o-visualization taxa-bar-plots.qzv

# Export taxonomy artifact to plain text format
qiime tools export \
  --input-path taxonomy.qza \
  --output-path exported-taxonomy

# Export feature table artifact to plain text format
qiime tools export \
  --input-path table.qza \
  --output-path exported-table
```

#### Combine sequences with taxonomy
```bash
# Combine representative sequences with taxonomy into a visualization
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --o-visualization rep-seqs-with-taxonomy.qzv
```

#### Make genus-level abundance table
```bash
# Collapse feature table to genus level (taxonomic level 6)
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table genus-table.qza

# Summarize the genus-level feature table
qiime feature-table summarize \
  --i-table genus-table.qza \
  --o-feature-frequencies genus-feature-frequencies.qza \
  --o-sample-frequencies genus-sample-frequencies.qza \
  --o-summary genus-table-summary.qzv

# Export genus-level feature table to BIOM format
qiime tools export \
  --input-path genus-table.qza \
  --output-path exported-genus-table

# Convert BIOM table to TSV format
biom convert \
  -i exported-genus-table/feature-table.biom \
  -o genus-table.tsv \
  --to-tsv
```

#### Generate donut graph of the identified taxa
The script prints the lowest assigned rank when no genus information is present.
Note: Sequences that are "unassigned" are not included in the visualization!

##### `make_donut_graph.py`
```python
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys
import os

# ---------------------------
# Parse command-line arguments
# ---------------------------
parser = argparse.ArgumentParser(
    description="Generate a donut chart from a QIIME2 taxonomy table with fallback to lowest assigned rank."
)

parser.add_argument(
    "--manifest",
    required=True,
    help="Path to manifest.tsv (must contain 'sample-id' column)"
)

parser.add_argument(
    "--table",
    required=True,
    help="Path to taxonomy table (e.g., genus-table.tsv)"
)

parser.add_argument(
    "--threshold",
    required=True,
    type=float,
    help="Relative abundance threshold as a decimal fraction (e.g., 0.03 for 3%)"
)

args = parser.parse_args()

manifest_path = args.manifest
table_path = args.table
threshold = args.threshold

# ---------------------------
# Validate input files
# ---------------------------
if not os.path.exists(manifest_path):
    sys.exit(f"Error: manifest file not found: {manifest_path}")

if not os.path.exists(table_path):
    sys.exit(f"Error: taxonomy table not found: {table_path}")

if not 0 <= threshold <= 1:
    sys.exit("Error: --threshold must be between 0 and 1 (e.g., 0.03 for 3%).")

# ---------------------------
# Read sample name from manifest.tsv
# ---------------------------
manifest = pd.read_csv(manifest_path, sep="\t")
manifest.columns = manifest.columns.str.strip()

if "sample-id" not in manifest.columns:
    sys.exit(
        f"Error: 'sample-id' column not found in {manifest_path}. "
        f"Available columns: {list(manifest.columns)}"
    )

sample_name = str(manifest.loc[0, "sample-id"]).strip()

# ---------------------------
# Read taxonomy table
# ---------------------------
df = pd.read_csv(table_path, sep="\t", skiprows=1)
df.columns = df.columns.str.strip()

if "#OTU ID" not in df.columns:
    sys.exit(
        f"Error: '#OTU ID' column not found in {table_path}. "
        f"Available columns: {list(df.columns)}"
    )

df = df.rename(columns={"#OTU ID": "Taxon"})

if sample_name not in df.columns:
    sys.exit(
        f"Error: Sample '{sample_name}' not found in {table_path}. "
        f"Available columns: {list(df.columns)}"
    )

df[sample_name] = pd.to_numeric(df[sample_name], errors="coerce").fillna(0)

# ---------------------------
# Taxonomy parsing logic
# ---------------------------
RANK_NAMES = {
    "k": "Kingdom",
    "p": "Phylum",
    "c": "Class",
    "o": "Order",
    "f": "Family",
    "g": "Genus",
    "s": "Species",
}

def extract_lowest_assigned_taxon(taxon):
    taxon = str(taxon).strip()

    if taxon.lower().startswith("unassigned"):
        return "Unassigned"

    parts = [p.strip() for p in taxon.split(";")]

    assigned = []
    genus_value = None

    for part in parts:
        if "__" not in part:
            continue

        prefix, value = part.split("__", 1)
        prefix = prefix.strip()
        value = value.strip()

        if not value or value == "_":
            continue

        assigned.append((prefix, value))

        if prefix == "g":
            genus_value = value

    if genus_value:
        return genus_value

    if assigned:
        lowest_prefix, lowest_value = assigned[-1]
        rank_name = RANK_NAMES.get(lowest_prefix, lowest_prefix)
        return f"{lowest_value} ({rank_name})"

    return "Unassigned"

df["AssignedTaxon"] = df["Taxon"].apply(extract_lowest_assigned_taxon)

# ---------------------------
# Aggregate counts
# ---------------------------
taxon_counts = df.groupby("AssignedTaxon")[sample_name].sum()
taxon_counts = taxon_counts[taxon_counts > 0]

# Remove unassigned taxa
taxon_counts = taxon_counts[taxon_counts.index != "Unassigned"]

if taxon_counts.empty:
    sys.exit("Error: No positive counts found after removing unassigned taxa.")

rel = taxon_counts / taxon_counts.sum()

# ---------------------------
# Group minor taxa
# ---------------------------
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
    autopct=lambda p: f"{p:.1f}%" if p >= threshold * 100 else "",
    startangle=90,
    wedgeprops={"width": 0.4, "edgecolor": "white", "linewidth": 1}
)

centre_circle = plt.Circle((0, 0), 0.60, fc="white")
ax.add_artist(centre_circle)

ax.text(0, 0, sample_name, ha="center", va="center", fontsize=12)

ax.legend(
    wedges,
    major.index,
    title="Assigned taxon",
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    frameon=False
)

ax.set_title(f"Taxonomic composition of {sample_name}")
plt.tight_layout()

# Save figures
plt.savefig(f"{sample_name}_donut.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{sample_name}_donut.svg", bbox_inches="tight")

# ---------------------------
# Print summary table
# ---------------------------
summary = pd.DataFrame({
    "Count": taxon_counts,
    "RelativeAbundance": rel
}).sort_values("RelativeAbundance", ascending=False)

print(summary)
```

##### Execute `make_donut_graph.py`
```bash
python make_donut_graph.py --manifest manifest.tsv --table genus-table.tsv --threshold 0.01
```