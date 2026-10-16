# Future Work

## Overview

The current pipeline establishes a reproducible workflow for filtering natural products, comparing them with the FKBP12-facing region of rapamycin, characterizing chemical space, docking a stratified candidate set, and modeling Vina scores.

The next phase should focus on determining whether the highest-ranked candidates remain plausible after more rigorous structural, chemical, and biological evaluation.

## 1. Candidate Identity and Quality Control

Before additional computation, the highest-ranked structures should be checked against their source records.

Recommended checks include:

- canonical SMILES and InChIKey verification;
- stereochemistry and salt-form review;
- tautomer and protonation-state enumeration;
- removal of duplicates and disconnected structures;
- identification of reactive or unstable functional groups;
- confirmation of natural-product provenance;
- commercial availability or synthetic feasibility.

This step is especially important for candidates represented only by database identifiers or incomplete names.

## 2. Docking Validation

The docking protocol should be benchmarked before interpreting score differences among candidates.

Recommended validation includes:

- redocking rapamycin into the FKBP12 pocket;
- comparing the predicted pose with the crystallographic pose;
- measuring heavy-atom root-mean-square deviation;
- repeating docking with multiple random seeds;
- testing exhaustiveness values such as 8 (retest), 16, and 32;
- evaluating sensitivity to docking-box placement and dimensions;
- inspecting hydrogen bonds, hydrophobic contacts, and steric clashes;
- including reference compounds and property-matched decoys.

Candidate ranking should consider pose plausibility and reproducibility in addition to the best Vina score.

## 3. Receptor Ensembles

The current workflow uses a single FKBP12 structure. Protein flexibility could be addressed by docking candidates against an ensemble assembled from:

- additional experimental FKBP12 structures;
- representative conformations from molecular dynamics;
- alternative side-chain conformations near the binding pocket.

Consensus across receptor states would provide stronger evidence than a favorable result from one rigid structure.

## 4. Ternary-Complex Evaluation

FKBP12 binding alone does not establish rapamycin-like activity. The most important mechanistic extension is to evaluate whether an FKBP12-bound candidate presents a compatible surface to the FRB domain of mTOR.

Possible approaches include:

- structural alignment with the FKBP12–rapamycin–FRB complex;
- steric and interface analysis after placing FRB near the docked complex;
- protein–protein interface scoring;
- restrained docking of FRB to candidate-bound FKBP12;
- molecular dynamics of the proposed ternary complex.

This stage could distinguish potential rapalogs from compounds that merely occupy the FKBP12 pocket.

## 5. Molecular Dynamics

Molecular dynamics could be applied to a small number of prioritized complexes to examine whether docked poses remain stable under an explicit-solvent model.

Potential analyses include:

- ligand and binding-site root-mean-square deviation;
- ligand–residue contact occupancy;
- hydrogen-bond persistence;
- conformational changes in FKBP12;
- ligand displacement from the pocket;
- comparison across independent simulation replicas.

Molecular dynamics would test pose stability but would not, by itself, establish biological activity.

## 6. Binding-Energy Refinement

MM-GBSA (Molecular Mechanics Generalized Born Surface Area) or MM-PBSA (Molecular Mechanics Poisson-Bolzman Surface Area) calculations could be used as approximate secondary ranking methods after molecular dynamics, helping identify a single core scaffold. Their results should be interpreted as model-dependent estimates rather than experimentally validated free energies.

Alchemical free-energy methods, including free-energy perturbation or thermodynamic integration, may become useful later for comparing closely related compounds with a shared binding mode. They are not well suited to the initial comparison of chemically diverse natural products without an established congeneric series.

## 7. Improved Machine Learning

The present linear models provide interpretable baselines but do not explicitly encode three-dimensional ligand–protein interactions.

Potential extensions include:

- nonlinear tree-based models;
- graph neural networks for molecular representation;
- molecular-language-model embeddings;
- protein–ligand interaction fingerprints;
- three-dimensional geometric neural networks;
- uncertainty estimation and applicability-domain analysis;
- nested cross-validation and external validation.

More complex models should be compared against the existing baselines and justified by improved held-out performance and interpretability.

## 8. ADMET and Chemical-Liability Assessment

High-ranking candidates should be evaluated for absorption, distribution, metabolism, excretion, toxicity (ADMET), and chemical liabilities.

Relevant assessments include:

- solubility and permeability;
- metabolic stability;
- cytochrome P450 interactions;
- P-glycoprotein interactions;
- hERG liability;
- mutagenicity and hepatotoxicity predictions;
- pan-assay interference and aggregation risk;
- reactive functional groups and nonspecific covalent chemistry.

These predictions should be treated as prioritization evidence rather than automatic pass/fail decisions.

## 9. Experimental Validation

Computational ranking should ultimately be tested experimentally.

A staged validation strategy could include:

1. confirmation of compound identity and purity;
2. direct FKBP12-binding assays;
3. competition with rapamycin;
4. structural or biophysical confirmation of binding;
5. assessment of FKBP12–compound–FRB ternary-complex formation;
6. measurement of downstream mTOR complex 1 signaling;
7. evaluation of autophagy-related responses;
8. cytotoxicity and off-target testing.

Agreement between independent assay types would provide stronger evidence than repeating the same measurement principle.

## 10. Reproducibility and Deployment

The software workflow can be strengthened by:

- adding unit and integration tests;
- recording software and database versions;
- saving run-specific configuration snapshots;
- recording random seeds and file checksums;
- publishing compact example data;
- distributing large checkpoints through versioned release assets;
- adding continuous-integration checks;
- packaging reusable pipeline stages as command-line functions.

These improvements would make the repository easier to evaluate, reproduce, and adapt to other structure-based natural-product discovery projects.