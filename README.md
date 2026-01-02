# GDC-Data-Preprocessing-Pipeline
This script downloads open-access RNA-Seq gene expression data from the NCI Genomic Data Commons (GDC) [https://gdc.cancer.gov](https://gdc.cancer.gov) for selected cancer primary sites and generates expression matrices (TPM, FPKM, and FPKM-UQ) for tumor and normal samples separately.

Data are queried from the GDC Files and Data APIs, downloaded in TSV format, parsed, filtered, and finally written to CSV matrices.

---

## Features

* Queries RNA-Seq Gene Expression Quantification metadata from GDC
* Downloads expression files in parallel using `ThreadPoolExecutor`
* Automatically handles gzip-compressed or plain TSV files
* Filters out non-coding `N_` gene entries
* Extracts required expression fields:

  * `tpm_unstranded`
  * `fpkm_unstranded`
  * `fpkm_uq_unstranded`
* Builds expression matrices where:

  * rows = genes
  * columns = samples
* Saves outputs for each cancer site and sample type:

  * `tumor_tpm.csv`
  * `tumor_fpkm.csv`
  * `tumor_fpkm_uq.csv`
  * `normal_tpm.csv`
  * `normal_fpkm.csv`
  * `normal_fpkm_uq.csv`

---

## Output Structure

Output files are written under:

```
tcga_csvs/<PRIMARY_SITE>/
```

Example:

```
tcga_csvs/Breast/
 ├── tumor_tpm.csv
 ├── tumor_fpkm.csv
 ├── tumor_fpkm_uq.csv
 ├── normal_tpm.csv
 ├── normal_fpkm.csv
 └── normal_fpkm_uq.csv
```

Each matrix contains:

| gene_id | gene_name | Sample1 | Sample2 | … |

---

## Requirements

Install dependencies:

```bash
pip install requests pandas
```

Python ≥ 3.8 is recommended.

---

## How to Run

Run the script directly:

```bash
python script.py
```

By default it processes multiple GDC primary sites defined in:

```python
generate_expression_csvs([
    "Adrenal Gland", "Bladder", "Bone Marrow and Blood", "Brain",
    "Breast", "Cervix", "Colorectal", "Esophagus", "Eye",
    "Head and Neck", "Kidney", "Liver", "Lung", "Lymph Nodes",
    "Ovary", "Pancreas", "Pleura", "Prostate", "Rectum",
    "Skin", "Soft Tissue", "Stomach", "Testis",
    "Thymus", "Thyroid", "Uterus"
])
```

To process a smaller subset, edit the list, e.g.:

```python
generate_expression_csvs(["Breast", "Lung"])
```

---

## Key Functions

### `get_rna_seq_metadata(primary_site)`

Queries the GDC Files API and returns metadata for RNA-Seq STAR-Counts gene expression files for the given primary site.

### `process_file(file_id, sample_id, file_name, data_endpoint)`

Downloads and parses a single expression file, returning TPM/FPKM expression vectors.

### `generate_expression_csvs(cancer_site_list)`

Orchestrates metadata retrieval, parallel downloading, matrix assembly, and CSV export.


