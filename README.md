# Analysis and code for bedbug male reproductive proteomics

This repository contains the analysis scripts used in the study: **"Evolution and adaptations of the seminal proteome in an insect with traumatic insemination."**

------------------------------------------------------------------------

## Repository Structure

These scripts require raw data available from the **Open Access
Repository of Saxon Universities ([OPARA](https://opara.zih.tu-dresden.de/home))**: [DOI: https://doi.org/10.25532/OPARA-1157](https://doi.org/10.25532/OPARA-1157)

------------------------------------------------------------------------

### `code/` contains scripts used for genomic and bioinformatic analyses.

1.  **cimex2dmel_RBB_plus.py**\
    Finds reciprocal best BLASTp hits from input `.fasta` files containing protein sequences. Points to `.fasta` files containing *Cimex lectularius* S-Lap orthologs and *Drosophila melanogaster* S-Laps. Run using `python3` with: `python3 cimex2dmel_RBB_plus.py` 

2.  **region_constrained_evidenced.py**\
    Aligns sequences and searches for M17 leucylaminopeptidase catalytic domain based on known reference position from grammy-smith (UniProt: Q9V3D8) or LAP3 (P00727). The only data needed to run this analysis is a single .fasta file containing protein sequences. Run using `python3` with: `python region_constrained_evidenced.py M17LAP3_SLAP12_GrSm_orthologs.fasta`

3.  **zinc_chimerax_coord.py**\
    Measures zinc ion coordination for proteins complexed with two zinc ions. Points to directory containing collected `.cif` files output from AlphaFold3. Run using ChimeraX in headless mode with: `ChimeraX --nogui --script zinc_chimerax_coord.py`


### `analysis/male_rt/` contains R scripts used to perform statistical analyses.

Users should download the raw data and update local file paths in the scripts accordingly. Outputs from each step should be saved to an `output/` directory for use
in downstream analyses.

1.  **male_rt_MSstats.R**\
    Main script for performing differential abundance analysis and downstream analyses characterising male reproductive proteins.

2.  **SLAP_tree.R**\
    Create bootstrapped phylogenetic tree for S-Lap protein sequences. Reads in the MAFFT alignment that was run locally with: `mafft --auto M17LAP3_SLAP12_GrSm_orthologs.fasta > M17SLAP3_SLAP12_GrSm_aligned.fasta`

5.  **topGO_results.R**\
    Performs Gene Ontology enrichment analyses. Run after creating outputs from `male_rt_MSstats.R`

------------------------------------------------------------------------
