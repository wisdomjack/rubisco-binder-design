

## Phylogenetic Analysis of Rubisco (rbcL)

A maximum-likelihood phylogeny of the Rubisco large subunit (rbcL) was reconstructed to examine evolutionary relationships among representative monocot and dicot plant species. Protein sequences were aligned using MAFFT (L-INS-i algorithm; --localpair --maxiterate 1000), which incorporates local pairwise alignment information and iterative refinement to maximize alignment accuracy for moderately divergent sequences.

Phylogenetic inference was performed using IQ-TREE 3 with automated model selection via ModelFinder. The best-fitting substitution model according to the Bayesian Information Criterion (BIC) was JTTDCMUT+I+G4, indicating a JTT-derived matrix with empirical corrections, a proportion of invariant sites, and gamma-distributed rate heterogeneity. Branch support was assessed using both ultrafast bootstrap approximation (UFBoot, 1000 replicates) and SH-like approximate likelihood ratio test (aLRT, 1000 replicates), providing complementary measures of node confidence.

The resulting tree was rooted using Selaginella moellendorffii, a lycophyte that serves as an appropriate outgroup to angiosperms. Nodes with low statistical support (UFBoot < 50) were collapsed to reduce overinterpretation of poorly resolved relationships. To improve visualization of deep and shallow divergences within a single framework, branch lengths were log-transformed prior to rendering. The tree was visualized in a circular layout using the ggtree framework, with clade membership annotated for interpretability.

The topology reveals a clear and well-supported separation between monocot and dicot lineages, consistent with established angiosperm phylogeny. Monocot species, primarily grasses and related taxa, cluster distinctly from dicot species, which include a diverse set of eudicot crops and model organisms. The high proportion of conserved sites in the alignment (~74%) reflects the strong evolutionary constraint on Rubisco, a central enzyme in carbon fixation. Notably, several closely related species exhibit identical sequences, which were retained to preserve biological context rather than artificially reducing redundancy.


<img width="3400" height="3400" alt="rubisco_phylogeny_circular" src="https://github.com/user-attachments/assets/bea4bd3e-93f0-407a-84b7-341cae650230" />

