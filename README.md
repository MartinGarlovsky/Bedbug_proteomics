# Analysis and code for bedbug male reproductive proteomics

This repository contains the analysis scripts used in the study: **"Evolution and adaptations of the seminal proteome in an insect with traumatic insemination."**

------------------------------------------------------------------------

## Repository Structure

### `code/` contains scripts used for genomic and bioinformatic analyses.

### `analysis/male_rt/` contains R scripts used to perform statistical analyses.

These scripts require raw data available from the **Open Access
Repository of Saxon Universities ([OPARA](https://opara.zih.tu-dresden.de/home))**: [DOI: https://doi.org/10.25532/OPARA-1157](https://doi.org/10.25532/OPARA-1157)

------------------------------------------------------------------------

## Analysis Workflow

Run the scripts in the `analysis/` directory in the following order. 

Users should download the raw data and update local file paths in the scripts accordingly. Outputs from each step should be saved to an `output/` directory for use
in downstream analyses.

1.  **male_rt_MSstats.R**\
    Main script for performing differential abundance analysis and downstream analyses characterising male reproductive proteins.

2.  **SLAP_tree.R**\
    Create bootstrapped phylogenetic tree for S-Lap protein sequences. Reads in the MAFFT alignment that was run locally with: `mafft --auto M17LAP3_SLAP12_GrSm_orthologs.fasta > M17SLAP3_SLAP12_GrSm_aligned.fasta`

5.  **topGO_results.R**\
    Performs Gene Ontology enrichment analyses. Run after creating outputs from `male_rt_MSstats.R`

------------------------------------------------------------------------
