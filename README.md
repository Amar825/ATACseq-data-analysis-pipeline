# ATAC-Seq Analysis: Yeast Hybridization Study

This repository documents a full end-to-end ATAC-seq analysis performed on yeast (**Saccharomyces cerevisiae** and its hybrid with **S. uvarum**). The analysis includes data download, preprocessing, alignment, peak calling, differential accessibility analysis, and visualization — all done using bash, R, and open-source bioinformatics tools.

---
## Objectives

This project was designed as a learning workflow with real biological data, aiming to:

- Understand chromatin accessibility changes in a yeast hybridization model
- Practice reproducible ATAC-seq data processing
- Explore downstream analysis with `DiffBind` and `MACS2`

---

## Experimental Design

- Dataset: SRA Project **PRJNA576599**
- 3 replicates of parental *S. cerevisiae*
- 3 replicates of hybrid *S. cerevisiae x S. uvarum*
- Sequencing: Illumina HiSeq 2500, paired-end 50 bp reads
- Reads/sample: ~60 million (high coverage)

---

## Directory Structure

```bash
raw_data/               # Raw FASTQ files
raw_data/trimming/      # Trimmed reads
ref_gen/                # Reference genomes + annotations
mapping/                # BAM files, stats
fastqc_initial/         # Raw read QC
fastqc_trimming/        # Trimmed read QC
diffbind/               # R-based peak comparison and visualization
```

## Setup: Computing environment
Using Conda for reproducibility. Create environments:
```bash
  # Core ATAC-seq environment
conda create -n atacseq
conda activate atacseq

conda install -c bioconda macs2 fastqc multiqc bwa samtools bedtools picard igv

# DiffBind environment
conda create -n diffbind
conda activate diffbind
conda install -c bioconda bioconductor-diffbind
```
### 1. Data Retrieval
Raw ATAC-seq data for this project was downloaded from the SRA using [SRA Explorer](https://sra-explorer.info/) with the project ID `PRJNA576599`.
To download the data, use the following script (modified from SRA Explorer):

```bash
bash raw_data/renamed_sra_explorer.sh
```
This will retrieve all six paired-end FASTQ files and place them in the raw_data/ directory.
Once downloaded, generate a sample list with:
```bash
  ls SRR*gz | cut -f 1 -d "_" | sort | uniq > sample_ids.txt
```
This list will be used for looping through all samples in later steps (trimming, mapping, etc.).
