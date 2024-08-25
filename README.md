# TCGA + GTEX Data Processing Pipeline

This repository contains a set of scripts to process and combine TCGA (The Cancer Genome Atlas) and GTEX (Genotype-Tissue Expression) datasets. The goal is to prepare these datasets for downstream analyses, including gene expression analysis and clinical data integration.

## Features

- **Data loading and merging**: Efficiently handles and merges large TCGA and GTEX gene expression datasets.
- **Custom data segmentation**: Automatically separates and organizes data by study type and detailed categories.
- **Flexible directory management**: Automatically creates and navigates directory structures for data and results storage.

## Prerequisites

Ensure you have the following packages installed:

- `dplyr`
- `stringr`
- `data.table`
- `librarian` (for managing R package dependencies)

These can be installed using the following R command:

```R
install.packages(c("dplyr", "stringr", "data.table"))
```

You will also need the following files (available from TCGA and GTEX databases):

- `TcgaTargetGtex_gene_expected_count.gz`: This file contains gene expression data.
- `TcgaTargetGTEX_phenotype.txt.gz`: This file contains the phenotype and clinical data.
- `detailed_category.csv`: A reference file that maps the primary site to detailed categories.
- `TCGA_catagory.csv` and `tcga_gtex_catagory.csv`: These files categorize the TCGA and combined TCGA-GTEX datasets by site and study.

## Usage

### 1. Initial Setup

```R
# Clear the environment
rm(list = ls())

# Set up the main working directories
root_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
data_dir <- shelf_dir("data", root_dir)
results_dir <- shelf_dir("results", root_dir)
ref_dir <- shelf_dir("reference", root_dir)
```

### 2. Data Download and Preprocessing

The script checks for existing processed data files (`gtex_tcga_expr_ph.Rdata`). If these files are not present, it will load and preprocess the raw expression and phenotype data.

```R
# Load and preprocess data if not already processed
if (!("gtex_tcga_expr_ph.Rdata" %in% list.files())) {
    TcgaTargetGtex_gene_expected_count_file <- "TcgaTargetGtex_gene_expected_count.gz"
    counts <- data.table::fread(TcgaTargetGtex_gene_expected_count_file, data.table = FALSE)
    ph <- data.table::fread("TcgaTargetGTEX_phenotype.txt.gz", data.table = FALSE)

    # Additional steps to filter and process the data are performed here
    save(expr, ph, file = "gtex_tcga_expr_ph.Rdata")
}
```

### 3. Separating Data by Source

After preprocessing, the script separates the data into GTEX and TCGA datasets and further categorizes them by primary sites for both datasets.

```R
# Separate GTEX data
ph_gtex <- ph[ph$`_study` == "GTEX", ]
# Separate TCGA data
ph_TCGA <- ph[ph$`_study` == "TCGA", ]
```

### 4. Further Data Segmentation

For each major location in the dataset, the script will save the corresponding expression and phenotypic data to separate files. Specifically, the correspondence between the main locations in the GTEX data and their actual location identifiers comes from `detailed_category.csv`, while the correspondence between the projects in the TCGA data and their actual location identifiers comes from `TCGA_category`.

```R
# Save segmented data by primary site for GTEX and TCGA
for (primary_site in primary_sites) {
    # Save GTEX data
    primary_site_dir <- shelf_dir(primary_site, gtex_dir)
    save(selected_ph, selected_expr, file = file.path(primary_site_dir, paste0(primary_site, "_expr_ph.Rdata")))

    # Save TCGA data
    catagory_dir <- shelf_dir(catagory_name, tcga_dir)
    save(selected_ph, selected_expr, file = file.path(catagory_dir, paste0(catagory_name, "_expr_ph.Rdata")))
}
```

### 5. Combined TCGA and GTEX Data

The correspondence between the GTEX and TCGA data is derived from `tcga_gtex_category.csv`. Additionally, it includes the actual number of cancer samples, adjacent cancer samples, and normal samples.

```R
# Combine TCGA and GTEX data based on categories
for (i in 1:nrow(tcga_gtex_catagory)) {
    # Process and save the combined data
    catagory_dir <- shelf_dir(catagory_name_brief, tcga_gtex_dir)
    save(selected_ph, selected_expr, file = file.path(catagory_dir, paste0(catagory_name, "_expr_ph.Rdata")))
}
```

## File Structure

- **data/**: Raw and processed data files.
- **results/**: Results from data processing.
- **reference/**: Reference files such as category mappings for TCGA and GTEX datasets.

## License

This project is licensed under the MIT License.
