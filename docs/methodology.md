# Rapalog Discover Pipeline Methodology

## Overview
Rapamycin is an FDA-approved pharmaceutical used clinically as an immunosuppressant following organ transplantation. In addition to its clinical applications, rapamycin has been shown to extend lifespan and improve healthspan in multiple model organisms. These effects are believed to arise primarily through inhibition of the mechanistic Target of Rapamycin Complex 1 (mTORC1), a central regulator of cellular growth, metabolism, and autophagy.

The biological activity of rapamycin involves a two-step binding process. First, rapamycin binds intracellularly to the FKBP12 protein. The resulting FKBP12–rapamycin complex subsequently binds the FKBP12-Rapamycin Binding (FRB) domain of mTORC1, leading to suppression of mTORC1 signaling.

Despite its therapeutic potential, rapamycin exhibits known side effects including immunosuppression, insulin resistance, and increased risk of metabolic dysfunction. Consequently, significant effort has been devoted to identifying rapamycin analogs ("rapalogs") that preserve beneficial modulation of mTOR signaling while reducing undesirable side effects.

Computational approaches provide a scalable method for prioritizing candidate molecules from extremely large chemical libraries prior to experimental validation. The objective of this project was to identify candidate natural products with structural and functional characteristics consistent with FKBP12 binding by combining cheminformatics, machine learning, unsupervised and supervised learning, linear regression, and molecular docking methods.

The workflow screened approximately 740,000 natural products from the COCONUT database and progressively narrowed this search space through molecular filtering, similarity analysis, clustering, predictive modeling, and docking-based prioritization.

## Data Sources
### FKBP12–Rapamycin Complex    
The experimentally determined FKBP12–rapamycin structure was obtained from the Protein Data Bank (PDB). Structural analysis was performed to identify rapamycin atoms involved in interactions with FKBP12 and to extract a binding-region fragment representing the pharmacophore most relevant to FKBP12 recognition.

### COCONUT Natural Product Database
Candidate molecules were obtained from the Collection of Open Natural Products (COCONUT) database containing approximately 740,000 naturally occurring compounds represented as SMILES strings.

SMILES structures were converted to RDKit molecular objects and subjected to quality-control filtering. Molecules containing invalid structures, duplicate entries, or missing descriptor information were removed. Additional filtering based on molecular weight, hydrogen-bonding characteristics, lipophilicity, ring count, and topological polar surface area reduced the dataset to approximately 97,000 candidate molecules.

## FKBP12 Binding Fragment Extraction
### Objective
The FKBP12–rapamycin complex was analyzed to identify rapamycin atoms participating in direct interactions with the FKBP12 binding pocket. The full rapamycin molecule contains structural features not directly involved in FKBP12 recognition. To focus similarity analysis on the pharmacologically relevant region of the molecule, a binding-region fragment was extracted from the FKBP12-rapamycin crystal structure.

An interaction map was constructed from contacting atoms identified within the crystal structure. These interaction-associated atoms were mapped to RDKit molecular indices and used to extract a rapamycin binding fragment representing the structural region most directly involved in FKBP12 recognition.

The extracted fragment was sanitized and validated using RDKit prior to downstream similarity analysis.

### Procedure
1. Load FKBP12-rapamycin crystal structure from the Protein Data Bank.
2. Identify rapamycin atoms participating in FKBP12 interactions.
3. Construct an interaction map based on spatial proximity.
4. Classify interaction types.
5. Map interaction-associated atoms to RDKit molecular indices.
6. Extract the corresponding rapamycin fragment.
7. Sanitize and validate the resulting molecular structure.

### Rationale
Using the FKBP12-interacting region rather than the entire rapamycin molecule focuses similarity calculations on structural elements most likely to contribute to FKBP12 binding and may improve identification of candidate rapalogs possessing alternative molecular scaffolds.

## Candidate Molecule Processing
### SMILES Processing
Candidate molecules were obtained from the Collection of Open Natural Products (COCONUT) database represented as SMILES (Simplified Molecular Input Line Entry System) strings.  RDKit served as the primary cheminformatics framework for molecular representation, descriptor generation, fingerprinting, and downstream analysis. 

SMILES strings were converted into RDKit molecular objects to enable downstream cheminformatics analyses. Molecules that failed RDKit parsing or sanitization were excluded from further processing. Duplicate molecular entries and structures with incomplete descriptor information were also removed.

The resulting RDKit molecule objects served as the canonical molecular representation used throughout the remainder of the analysis pipeline.

### Molecular Filtering
The initial COCONUT dataset contains a broad range of natural products, many of which fall outside physicochemical property ranges typically associated with drug-like molecules.

To focus the search space on compounds with rapamycin-like characteristics, candidate molecules were filtered using descriptor-based criteria chosen to identify compounds possessing rapamycin-like physicochemical characteristics:

  - Molecular Weight: 350–700 Da
  - Hydrogen Bond Acceptors: 3–8 
  - Hydrogen Bond Donors: 1–4
  - LogP: 2–5
  - Ring Count: 2–5
  - TPSA: <150 Å²

The filtering strategy was designed to remove compounds that were unlikely to function as viable FKBP12-binding candidates while preserving sufficient chemical diversity for subsequent similarity analysis.

Application of these filters reduced the candidate pool from approximately 740,000 molecules to approximately 97,000 molecules.

### Descriptor Generation
Molecular descriptors were calculated using RDKit for each candidate molecule.

Descriptors included:

  - Molecular Weight (MW)
  - Hydrogen Bond Acceptors (HBA)
  - Hydrogen Bond Donors (HBD)
  - LogP
  - Ring Count
  - Aromatic Ring Count
  - Topological Polar Surface Area (TPSA)
  - Amide Count

These descriptors served two purposes. First, they were used to characterize the physicochemical properties of candidate molecules and support clustering analyses. Second, they were used as input features for supervised regression models predicting molecular similarity and binding affinity.


## Molecular Representation, Fingerprinting, and Similarity Analysis
### Morgan Fingerprints
Candidate molecules and the extracted rapamycin binding fragment were represented using Morgan circular fingerprints generated with RDKit.

Morgan fingerprints encode local molecular neighborhoods around each atom and provide a compact representation of molecular structure suitable for similarity calculations and virtual screening workflows.

Fingerprint generation was evaluated using multiple hyperparameter settings, including radii of 1 and 2 and a fingerprint length of 2048 bits. Radii of 1 and 2 were evaluated to compare local versus broader neighborhood representations of molecular structure. Fingerprint length of 2048 bits was chosen to balance computational efficiency which might be compromised with 4096 bits, and hash table collisions which might be excessive with only 1024 bits. 

The extracted rapamycin binding fragment was converted into a reference fingerprint against which all candidate molecules were compared.

### Tanimoto Similarity
Tanimoto similarity was used as the primary measure of structural similarity between candidate molecules and the rapamycin binding fragment.

The Tanimoto coefficient measures the fraction of fingerprint features shared between two molecules relative to the total number of unique features present in either molecule.

Because Tanimoto similarity is symmetric and widely used throughout cheminformatics and virtual screening applications, it provides a useful baseline for evaluating structural similarity.

All candidate molecules were scored against the rapamycin binding fragment and ranked according to Tanimoto similarity.

### Tversky Similarity
Tversky similarity was additionally evaluated to allow asymmetric weighting of molecular features.

Unlike Tanimoto similarity, Tversky similarity can independently control penalties associated with missing desired features and additional superfluous features. This property is particularly useful when comparing candidate molecules of varying size against a smaller pharmacophore-like reference structure such as the rapamycin binding fragment.

Two Tversky configurations were evaluated:

  - **Symmetric** Tversky similarity, which penalizes missing and additional features equally (α = 0.5, β = 0.5)
  - **Feature-focused** Tversky similarity, which prioritizes retention of rapamycin-like binding features while reducing penalties for additional molecular complexity (α = 0.1, β = 0.9).

Candidate molecules were ranked according to both Tanimoto and Tversky similarity scores. Correlation analyses were performed to evaluate whether different similarity metrics produced consistent candidate prioritizations across the filtered data set. 

### Candidate Ranking
Structural similarity was used as a prioritization strategy for downstream analysis and docking studies rather than a direct predictor of biological activity. 

Following fingerprint generation and similarity calculation, candidate molecules were ranked according to their similarity to the extracted rapamycin binding fragment. Rankings were generated using Tanimoto similarity as the primary metric and compared against rankings produced by both **symmetric** and **feature-focused** Tversky similarity.

The objective of candidate ranking was not to identify exact structural replicas of rapamycin, but rather to prioritize compounds retaining key FKBP12-binding characteristics while potentially exhibiting alternative molecular scaffolds. This approach preserves structural diversity while enriching the candidate pool for molecules likely to share functionally relevant features.

Similarity rankings were subsequently used in three ways:

1. Selection of high-priority candidate molecules for molecular docking.
2. Investigation of structure-function relationships through clustering analyses.
3. Construction of supervised regression models relating molecular descriptors to similarity and binding affinity.

Because similarity calculations are computationally inexpensive relative to molecular docking, candidate ranking served as an efficient first-stage screening strategy for reducing a large chemical search space to a tractable set of molecules for more computationally intensive analyses.

Rather than relying solely on similarity ranking, subsequent analyses evaluated whether molecular descriptors and similarity measures could predict FKBP12 binding affinity, enabling assessment of the relationship between structural similarity and functional activity.

### Similarity Method Determination
Tanimoto and Tversky (**symmetric** and **feature-focused**) similarity metrics were compared across the full dataset and a high-similarity scoring subset candidates to ensure both statistical robustness and enrichment for chemically relevant structures. While the metrics were highly correlated across all molecules (r = 0.97), this correlation decreased substantially among top 5000 candidates (r = 0.79).  Histogram analysis of top candidates shows a rightward shift (greater similarity generally) in Tversky **feature-focused** similarity scores compared to Tversky **symmetric** and Tanimoto methods, consistent with reduced penalty for additional molecular features and improved sensitivity to fragment containment.

This supports the use of Tversky **feature-focused** similarity for identifying pharmacophore-relevant structures in subsequent analyses.

## Unsupervised Learning
### Feature Engineering
The molecular descriptors generated during candidate processing were assembled into a feature matrix for downstream unsupervised learning analyses. Features were standardized using z-score normalization prior to dimensionality reduction and clustering analyses to ensure comparable scaling across descriptors with different numerical ranges.

The resulting feature matrix served as the basis for principal component analysis (PCA), K-Means clustering, and UMAP visualization.

### Dimensionality Reduction
PCA was applied to the standardized chemical feature matrix to identify the dominant sources of variance within the candidate molecule population.  The first two principal components captured approximately 57% of total feature variance.

PCA provided a lower-dimensional representation of chemical feature space, facilitating visualization, interpretation of feature relationships, and subsequent clustering analyses. Principal component loadings were examined to identify the molecular descriptors contributing most strongly to each component.

PCA projections were also used to visualize cluster structure and investigate relationships between chemical properties and molecular similarity.

### Clustering
K-Means clustering was applied to the principal-component representation of the standardized chemical feature matrix to identify groups of chemically similar compounds.

The optimal number of clusters was evaluated using elbow plots, silhouette scores, and silhouette distributions. Silhouette analysis indicated that a two-cluster solution provided the strongest balance between cluster cohesion and separation and was therefore selected for primary downstream analyses.  The four-cluster solution was additionally explored to investigate potential additional chemical substructure within the candidate population. 

Cluster assignments were subsequently used to characterize chemical feature patterns, evaluate relationships between molecular similarity and chemistry, and guide stratified sampling for molecular docking studies.

### UMAP
Uniform Manifold Approximation and Projection (UMAP) was applied to visualize candidate molecules within a nonlinear low-dimensional embedding of chemical feature space. Unlike PCA, which captures linear variance, UMAP preserves local neighborhood structure and can reveal complex relationships among chemically related compounds. UMAP visualizations were used to evaluate cluster separation, identify potential subclusters, and compare cluster structure against PCA-based representations.


## Predicting Similarity from Chemical Features
### Dataset Construction
A supervised regression dataset was constructed using molecular descriptors as input features and molecular similarity scores as target values. Similarity scores were calculated using the rapamycin binding fragment as the reference structure. Descriptor features were standardized prior to model training, and training and testing datasets were generated to evaluate predictive performance.

### Similarity Regression Modeling
Multiple regression approaches were evaluated, including baseline linear regression, then including Elastic Net, Ridge, and Lasso regularization. These models were used to assess the extent to which conventional molecular descriptors could predict similarity to the rapamycin binding fragment. Specificially, Elastic Net was included to combine L1 (Lasso) and L2 (Ridge) regularization, enabling feature selection while maintaining coefficient stability in the presence of correlated molecular descriptors.  Model coefficients, residual analysis, and diagnostic plots were examined to identify descriptors most strongly associated with rapamycin-like structural characteristics and to evaluate predictive performance. This stage evaluated whether inexpensive molecular descriptors could serve as surrogate predictors of fingerprint-based similarity scores, potentially reducing the need for large-scale similarity calculations in future screening workflows.

## Molecular Docking
### Candidate Selection
Molecular docking was computationally impractical for the full filtered candidate set. A stratified sampling strategy was therefore employed to select 1,000 candidate molecules for docking, with 971 molecules successfully completing ligand preparation and docking workflows. Sampling was performed proportionally across the two major chemical clusters identified during unsupervised learning (~60/40 distribution) and across similarity-score quartiles within each cluster. This approach preserved both chemical diversity and representation across the similarity-ranking spectrum while maintaining computational tractability.

### Docking Box Construction
The FKBP12 binding site was extracted from the experimentally determined FKBP12–rapamycin crystal structure. Contacting atoms associated with rapamycin binding were used to define the docking region. Docking box dimensions and coordinates were validated numerically and visually to ensure adequate coverage of the FKBP12 binding pocket.

### AutoDock Vina
AutoDock Vina was used to estimate binding affinity between candidate molecules and the FKBP12 binding pocket. Ligands were prepared through hydrogen addition, cleanup, and conversion to PDBQT format. The FKBP12 receptor was similarly prepared for docking. Docking parameters exhaustiveness=8 and n_poses=9 were used to balance screening throughput and computational cost, enabling evaluation of 971 candidate molecules on workstation-class hardware. Docking scores were used as computational estimates of binding affinity and do not substitute for experimental validation.

## Predicting Binding Affinity from Chemical Features + Similarity
### Dataset Construction
A second supervised regression dataset was constructed using molecular descriptors and similarity scores as predictor variables and AutoDock Vina binding affinity scores as target values. The resulting dataset represented the relationship between molecular structure, rapamycin binding fragment similarity, and predicted FKBP12 binding affinity.

### Regression Models
Multiple regression approaches were evaluated, including baseline linear regression, Elastic Net regression, Ridge regression, and Lasso regression. Again, Elastic Net was included to combine L1 (Lasso) and L2 (Ridge) regularization, enabling feature selection while maintaining coefficient stability in the presence of correlated molecular descriptors.  Regularized models were included to mitigate overfitting, improve coefficient stability, and identify the most informative molecular descriptors associated with binding affinity.

### Model Evaluation
Model performance was evaluated using multiple regression error metrics and diagnostic visualizations. Predicted versus observed binding affinity values, residual analyses, and coefficient comparisons were used to assess model quality and interpret relationships between molecular descriptors, similarity measures, and predicted binding affinity. This stage evaluated whether molecular descriptors and similarity measures could serve as surrogate predictors for computationally expensive docking calculations.

## Workflow Summary
1. Acquire FKBP12–rapamycin structure and COCONUT natural product database.
2. Extract rapamycin binding fragment from the FKBP12-Rapamycin interaction region.
3. Filter and standardize candidate molecules.
4. Generate molecular descriptors and Morgan fingerprints.
5. Compute Tanimoto and Tversky similarity scores.
6. Perform PCA, K-Means, and UMAP clustering/visualization.
7. Train supervised models to predict rapamycin similarity from chemical features.
8. Perform docking of a stratified subset of candidate molecules using AutoDock Vina.
9. Train regression models to predict binding affinity from chemical descriptors and similarity measures.
10. Prioritize candidate rapalogs for future experimental validation.


## Computational Environment
The analysis was performed in Python using:

  - RDKit
  - NumPy
  - pandas
  - scikit-learn
  - UMAP-learn
  - Matplotlib
  - AutoDock Vina
  - Open Babel

All workflows were implemented using reproducible notebook-based and Python module-based pipelines.


## Reproducibility
All data processing, feature generation, similarity calculations, clustering analyses, regression modeling, docking workflows, and visualization steps were implemented as reproducible Python workflows. Repository organization separates reusable computational modules from exploratory notebooks to facilitate reproducibility and future extension.


## Study Limitations
Several limitations should be noted. Binding affinity estimates were derived from AutoDock Vina and therefore represent computational predictions rather than experimentally measured binding. Docking was performed on a stratified subset of candidate molecules due to computational constraints, and alternative docking protocols, molecular dynamics simulations, and experimental validation may yield additional insights. Similarity-based screening may also fail to identify structurally distinct compounds exhibiting comparable functional activity.

