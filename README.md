
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![AI](https://img.shields.io/badge/AI-Protein%20Design-purple)





![phylogenetic_binder_comp_design_figures 001](https://github.com/user-attachments/assets/c6864535-ba8b-43b4-8a54-de53a3b2b2ed)

**Firgure 1 Assembly pathway of plant Rubisco.**
RbcL monomers dimerize and further assemble into an RbcL₈ core, which is then capped by eight RbcS subunits to form the mature L₈S₈ holoenzyme. Top and side views of the holoenzyme are shown (right). RbcL, green; RbcS, lavender.


## Phylogenetic Analysis of Rubisco (rbcL)

A maximum-likelihood phylogeny of the Rubisco large subunit (rbcL) was reconstructed to examine evolutionary relationships among representative monocot and dicot plant species. Protein sequences were aligned using MAFFT (L-INS-i algorithm; --localpair --maxiterate 1000), which incorporates local pairwise alignment information and iterative refinement to maximize alignment accuracy for moderately divergent sequences.

Phylogenetic inference was performed using IQ-TREE 3 with automated model selection via ModelFinder. The best-fitting substitution model according to the Bayesian Information Criterion (BIC) was JTTDCMUT+I+G4, indicating a JTT-derived matrix with empirical corrections, a proportion of invariant sites, and gamma-distributed rate heterogeneity. Branch support was assessed using both ultrafast bootstrap approximation (UFBoot, 1000 replicates) and SH-like approximate likelihood ratio test (aLRT, 1000 replicates), providing complementary measures of node confidence.

The resulting tree was rooted using Selaginella moellendorffii, a lycophyte that serves as an appropriate outgroup to angiosperms. Nodes with low statistical support (UFBoot < 50) were collapsed to reduce overinterpretation of poorly resolved relationships. To improve visualization of deep and shallow divergences within a single framework, branch lengths were log-transformed prior to rendering. The tree was visualized in a circular layout using the ggtree framework, with clade membership annotated for interpretability.

The topology reveals a clear and well-supported separation between monocot and dicot lineages, consistent with established angiosperm phylogeny. Monocot species, primarily grasses and related taxa, cluster distinctly from dicot species, which include a diverse set of eudicot crops and model organisms. The high proportion of conserved sites in the alignment (~74%) reflects the strong evolutionary constraint on Rubisco, a central enzyme in carbon fixation. Notably, several closely related species exhibit identical sequences, which were retained to preserve biological context rather than artificially reducing redundancy.


<img width="3400" height="3400" alt="rubisco_phylogeny_circular" src="https://github.com/user-attachments/assets/bea4bd3e-93f0-407a-84b7-341cae650230" />


**Firgure 2 Phylogenetic tree of angiosperm Rubisco RbcL subunit.**



# Evolutionary Conservation Analysis of Plant Rubisco

## Overview

To guide minibinder design toward the RbcL–RbcS interface pocket, we characterized the evolutionary conservation of both Rubisco subunits across 40 representative land plant species. Rather than treating conservation as a binary property, we quantified per-position chemical variability using a group-corrected Shannon entropy metric, computed over a chemical-property consensus sequence that was subsequently used as the structural template for AlphaFold 3 folding.

---

## Sequence Dataset

Protein sequences for the Rubisco large subunit (RbcL) and small subunit (RbcS) were compiled from 40 representative land plant species spanning monocots and dicots, including major crop species (*Oryza sativa*, *Zea mays*, *Triticum aestivum*, *Glycine max*, *Solanum lycopersicum*) and model organisms (*Arabidopsis thaliana*, *Nicotiana tabacum*, *Spinacia oleracea*). Sequences were sourced from UniProtKB and curated manually. Entries containing ambiguous amino acid characters (X, B, Z, U, O) were excluded, as were sequences below minimum length thresholds (400 aa for RbcL, 80 aa for RbcS), to ensure alignment quality and biological relevance.

---

## Multiple Sequence Alignment

Sequences were aligned independently for each subunit using MAFFT v7 with the L-INS-i algorithm (`--localpair --maxiterate 1000`), which incorporates local pairwise alignment information and iterative refinement. This algorithm is recommended for moderately divergent protein sequences where local structural homology is expected to dominate alignment accuracy. The resulting alignments comprised 40 sequences each, spanning 488 columns for RbcL and 135 columns for RbcS.

---

## Chemical Consensus Sequence

A chemical-property consensus sequence was derived from each alignment to serve as a structurally representative template for folding. For each alignment column, residues were first classified into six chemical property groups: hydrophobic (A, V, I, L, M), aromatic (F, W, Y), polar (S, T, N, Q), positively charged (K, R, H), negatively charged (D, E), and special (G, P, C). The dominant chemical group was identified by majority vote across non-gap residues, and the most frequent amino acid within that dominant group was selected as the consensus residue. This two-step procedure ensures that the consensus reflects the chemical character of the position rather than being biased by splits among chemically similar residues — for example, a column that is 60% hydrophobic (split among A, V, I, L) correctly yields a hydrophobic consensus rather than defaulting to the single most frequent individual residue.

Alignment columns in which more than 50% of sequences contained a gap were excluded entirely from the consensus, as these represent lineage-specific insertions absent from the majority of sequences and therefore inappropriate for a consensus structural template. This filtering step removed 11 columns from the RbcL alignment and 9 columns from the RbcS alignment, yielding consensus sequences of 477 aa (RbcL) and 126 aa (RbcS). These sequences were submitted to AlphaFold 3 as eight copies each (L₈S₈ stoichiometry) to generate the consensus structural model.

---

## Structure Prediction

The chemical consensus sequences for RbcL (477 aa) and RbcS (126 aa) were folded as a complete L₈S₈ holoenzyme using AlphaFold 3. Five models were generated; all achieved high confidence scores (pTM = 0.92, ipTM = 0.92, ranking score = 0.94–0.95), indicating robust prediction of both intra- and inter-chain contacts. Model 0 (ranking score = 0.95) was selected for all subsequent analyses and visualization. The predicted structure recapitulates the canonical Rubisco architecture, with eight RbcL subunits forming the catalytic core and eight RbcS subunits capping the top and bottom faces.

| Model | pTM | ipTM | Ranking score |
|-------|-----|------|---------------|
| 0 | 0.92 | 0.92 | 0.95 |
| 1 | 0.92 | 0.92 | 0.95 |
| 2 | 0.92 | 0.92 | 0.95 |
| 3 | 0.93 | 0.92 | 0.94 |
| 4 | 0.92 | 0.92 | 0.94 |

---

## Conservation Quantification

Per-position evolutionary conservation was quantified using Shannon entropy calculated over chemical group frequencies rather than individual amino acid frequencies. For each alignment column retained in the consensus (i.e., gap-dominant columns excluded), the entropy $H$ was computed as:

$$H = -\sum_{g} p_g \log_2 p_g$$

where $p_g$ is the frequency of chemical group $g$ among non-gap residues at that position. Computing entropy over the six chemical groups rather than the 20 individual amino acids provides a more biologically meaningful measure of functional constraint: positions that tolerate conservative substitutions within a chemical class (e.g., K↔R, or V↔L) are correctly scored as conserved, whereas standard per-residue entropy would overestimate their variability. The theoretical maximum entropy under this scheme is $\log_2(6) \approx 2.585$ bits, corresponding to a completely uniform distribution across all six chemical groups.

Observed maximum entropy was 2.03 bits for RbcL and 2.06 bits for RbcS. To enable direct visual comparison across both subunits — which is essential when examining the interface pocket formed jointly by RbcL and RbcS residues — entropy values were normalized to a single shared scale defined by the global observed maximum of **2.06 bits**. Normalization to the observed rather than theoretical maximum ensures that the full color range is utilized for the actual variability present in this dataset of 40 plant species, without implying that any position approaches complete chemical disorder.

---

## Structural Visualization

Conservation was mapped onto the consensus structural model using ChimeraX. Each residue was assigned a color by linearly interpolating through a nine-bin pastel color scale (matching the binning convention of ConSurf) ranging from mint green (fully conserved, entropy = 0 bits) through sky blue, periwinkle, and lavender, to hot pink (maximally variable, entropy = 2.06 bits). The scale was applied identically to both RbcL and RbcS chains, so that color values are directly comparable across the subunit interface. Surfaces were rendered in flat graphics mode with silhouette edges, and the final figure was exported at 4000 × 4000 px with 3× supersampling.

---

## Results and Discussion

The surface conservation map reveals a striking pattern consistent with the functional architecture of Rubisco. The large subunit (RbcL) core is predominantly conserved (mint green), reflecting the strong evolutionary constraint imposed by the catalytic mechanism of CO₂ fixation, which is essentially invariant across land plants. Variable positions (pink to hot pink) are enriched at the periphery of the holoenzyme and at inter-subunit contact surfaces, suggesting that these regions have tolerated substitution without compromising catalytic function.

The binder pocket — formed at the interface between one RbcS subunit and two adjacent RbcL subunits — displays an intermediate-to-low entropy profile, with the pocket floor contributed primarily by conserved RbcL residues and the pocket rim showing moderate variability from RbcS. This conservation pattern is favorable for minibinder design for two reasons. First, a conserved binding surface increases the likelihood that a binder designed against the consensus sequence will retain affinity across diverse plant Rubisco orthologs, a desirable property for a biosensor intended for broad deployment. Second, the presence of some variability at the RbcS-contributed rim provides chemical diversity that can potentially be exploited to achieve ortholog selectivity if required.

The use of a chemical-property consensus sequence as the folding template, rather than a single representative sequence, is deliberate: it biases the structural model toward the chemically central solution across the phylogenetic sampling, smoothing over idiosyncratic substitutions in any individual species. The high AlphaFold 3 confidence scores (ipTM = 0.92) across all five predicted models indicate that the consensus sequence folds robustly into the canonical Rubisco structure, validating this approach.

---







<img width="1920" height="1080" alt="conservation 001" src="https://github.com/user-attachments/assets/fad83e48-cfff-4e86-b138-31b57bddc360" />

**Figure 3. Evolutionary conservation of the Rubisco holoenzyme surface.**
*(Left)* Surface representation of the consensus L₈S₈ Rubisco holoenzyme showing the RbcL–RbcS–RbcL trimer interface. Subunits are colored by identity: RbcL, green; RbcS, lavender. The binder pocket at the intersubunit interface is indicated.
*(Right)* The same view colored by group Shannon entropy. Color scale ranges from mint green (entropy = 0 bits, fully conserved) to hot pink (entropy = 2.06 bits, maximally variable), normalized to the global observed maximum across both subunits. Entropy was computed over six chemical property groups (hydrophobic, aromatic, polar, positive, negative, special) from a 40-species alignment of land plant Rubisco sequences. The pocket floor is predominantly conserved, consistent with its suitability as a minibinder target.
