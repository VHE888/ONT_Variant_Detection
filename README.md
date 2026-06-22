## Human Whole Genome Variation Detection (Oxford Nanopore)

## Platform

**Oxford Nanopore PromethION P2 Solo**

- Used for human whole genome sequencing (WGS)
- No onboard computing functionality

---

## 1. Data Generation

### Basecalling

Tools:

- MinKNOW
- Dorado

Settings:

- Basecalling: ON
- Modified Bases Detection: ON (for epigenetic information)

### Input Data

Generated files:

- POD5
- FASTQ
- BAM (contains sequence information and methylation information)

---

## 2. Data Analysis

### Workflow

**EPI2ME: wf-human-variation**

Required inputs:

- Reference genome (GRCh38, FASTA)
- BAM file

The workflow supports **real-time analysis**. As BAM files are generated during sequencing, the pipeline can begin analysis immediately.

### Analysis Options

- SNVs (Single Nucleotide Variants)
- SVs (Structural Variants)
- STRs (Short Tandem Repeats / Repeat Expansions)
- CNVs (Copy Number Variations)
- Methylation / Epigenetic Analysis
- Phasing (Haplotype Phasing)

### Output Files

- HTML report
- CRAM
- bedMethyl
- VCF

### Recommended Requirements

- CPUs: 32
- Memory: 128 GB

---

## 3. Reports and Visualization

### Alignment Reports

Metrics include:

- Total reads
- Read N50
- Mean coverage
- Mapping rate

### Variant Reports

Reports for:

- SNVs
- SVs
- STRs
- CNVs
- Methylation
- Phasing

### Visualization

**IGV (Integrative Genomics Viewer)**

Files used for visualization:

- `.cram`
- `.bw`
- `.vcf`
