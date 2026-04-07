# `analysis/male_rt/` 

Users should download the raw data and update local file paths in the scripts accordingly. Outputs from each step should be saved to an `output/` directory for use
in downstream analyses.

1.  **male_rt_MSstats.R**\
    Main script for performing differential abundance analysis and downstream analyses characterising male reproductive proteins.

2.  **SLAP_tree.R**\
    Create bootstrapped phylogenetic tree for S-Lap protein sequences. Reads in the MAFFT alignment that was run locally with: `mafft --auto M17LAP3_SLAP12_GrSm_orthologs.fasta > M17SLAP3_SLAP12_GrSm_aligned.fasta`

5.  **topGO_results.R**\
    Performs Gene Ontology enrichment analyses. Run after creating outputs from `male_rt_MSstats.R`
