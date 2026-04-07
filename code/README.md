# Code

1.  **cimex2dmel_RBB_plus.py**\
    Finds reciprocal best BLASTp hits from input `.fasta` files containing protein sequences. Points to `.fasta` files containing *Cimex lectularius* S-Lap orthologs and *Drosophila melanogaster* S-Laps. Run using `python3` with: `python3 cimex2dmel_RBB_plus.py` 

2.  **region_constrained_evidenced.py**\
    Aligns sequences and searches for M17 leucylaminopeptidase catalytic domain based on known reference position from grammy-smith (UniProt: Q9V3D8) or LAP3 (P00727). The only data needed to run this analysis is a single .fasta file containing protein sequences. Run using `python3` with: `python region_constrained_evidenced.py M17LAP3_SLAP12_GrSm_orthologs.fasta`

3.  **zinc_chimerax_coord.py**\
    Measures zinc ion coordination for proteins complexed with two zinc ions. Points to directory containing collected `.cif` files output from AlphaFold3. Run using ChimeraX in headless mode with: `ChimeraX --nogui --script zinc_chimerax_coord.py`
