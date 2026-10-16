# Biological Interpretation

## Scientific Question

Rapamycin first binds FKBP12 through one region of its structure and then presents a separate surface to the FKBP–rapamycin-binding (FRB) domain of mTOR complex 1. This project asks whether natural products contain chemical features resembling the FKBP12-facing portion of rapamycin and whether selected compounds are predicted to bind favorably within the FKBP12 pocket.

The present pipeline addresses only the first part of that mechanism. It does not establish that a candidate can recruit the FRB domain, inhibit mTOR complex 1, promote autophagy, or reproduce rapamycin’s biological effects.

## Why Molecular Similarity Is Useful

Molecular similarity provides an efficient first-pass method for searching a large chemical collection. Morgan fingerprints encode local atomic environments, while Tanimoto and Tversky measures quantify the overlap between candidate and reference fingerprints.

The fragment-directed Tversky score is particularly useful here because it emphasizes whether a candidate contains features resembling the selected FKBP12-facing region of rapamycin without requiring the candidate to resemble the entire rapamycin molecule.

Similarity therefore provides a rational way to enrich the candidate set before more computationally expensive docking.

## Why Similarity Does Not Fully Predict Binding

Fingerprint similarity and molecular docking measure different properties.

Similarity is based primarily on two-dimensional molecular substructures. Docking additionally considers a proposed three-dimensional pose within the protein pocket, including molecular shape, steric compatibility, hydrogen bonding, hydrophobic contacts, and conformational flexibility.

This distinction is visible in the results. The highest-ranked docking candidates had fragment-directed similarity scores ranging from approximately 0.19 to 0.33 rather than uniformly possessing the highest similarity values in the screened collection. Compounds with moderate similarity could therefore receive favorable Vina scores if their three-dimensional structures were compatible with the FKBP12 pocket.

Conversely, a compound with high fingerprint similarity may dock less favorably because of its conformation, substituent orientation, protonation state, steric clashes, or other features not represented adequately by the fingerprint.

## Interpretation of the Regression Results

The similarity regression achieved a held-out test R² of 0.534. This indicates that the selected physicochemical descriptors captured a meaningful but incomplete portion of the structural information represented by the fragment-directed fingerprint score.

The Vina-score regression achieved a lower test R² of approximately 0.30. Ordinary linear regression, Lasso, Ridge, and Elastic Net produced similar results, and regularization did not materially improve prediction.

The lower docking-model performance is scientifically plausible. Molecular weight, cLogP, topological polar surface area, hydrogen-bond counts, ring counts, amide count, and fragment similarity summarize broad molecular properties, but they do not explicitly represent:

- three-dimensional ligand geometry;
- the orientation of chemical groups within the binding pocket;
- specific ligand–residue contacts;
- conformational strain;
- ligand and receptor flexibility;
- tautomeric and protonation states;
- water-mediated interactions;
- uncertainty in the Vina scoring function.

The regression models should therefore be viewed as interpretable surrogate models rather than replacements for docking.

## Chemical-Space Clusters

The two-cluster solution separated the filtered natural products partly by size and polarity. Cluster 0 contained compounds with lower average molecular weight, topological polar surface area, hydrogen-bonding capacity, and amide content. Cluster 1 contained larger and more polar compounds on average.

Stratifying the docking sample across clusters and similarity quartiles reduced the risk of evaluating only one narrow chemical class. This was important because favorable FKBP12 docking was not restricted to compounds with the highest fragment-similarity scores.

The clusters are descriptive groupings rather than experimentally validated biological classes. They should not be interpreted as distinct mechanisms of action.

## Interpretation of the Top Candidates

NPD14003, STL518675, CNP0253821.1, CNP0243466.1, CNP0342192.2, and Pentacecilide B were among the strongest-scoring candidates in the production docking run.

These compounds are computational leads, not confirmed FKBP12 binders. Before biological conclusions are drawn, each candidate should be evaluated for:

- chemical identity and stereochemistry;
- source-database reliability;
- protonation and tautomeric state;
- docking-pose plausibility;
- repeated docking and scoring stability;
- problematic or reactive functional groups;
- physicochemical and ADMET liabilities;
- availability or feasibility of synthesis;
- experimental FKBP12 binding;
- ability to support or disrupt interactions involving the FRB domain.

## Relationship to Rapamycin-Like Activity

Binding FKBP12 is necessary for the canonical rapamycin mechanism, but it is not sufficient. A candidate could bind FKBP12 without creating the molecular surface required for FRB recruitment. It could also bind in an orientation that blocks ternary-complex formation or produces an entirely different biological effect.

Accordingly, the current candidates should be described as potential FKBP12-binding compounds rather than rapalogs. Demonstrating rapalog-like activity would require evidence of ternary-complex formation and downstream modulation of mTOR complex 1.

## Biological Conclusions

The study supports three main conclusions:

1. Fragment-directed similarity can reduce a large natural-product collection to a chemically relevant candidate space.
2. Similarity and docking provide complementary information; neither should be used alone to infer biological activity.
3. Simple molecular descriptors provide only limited prediction of docking scores, reinforcing the importance of three-dimensional structural analysis and experimental validation.

The pipeline is therefore best understood as a reproducible hypothesis-generation system for prioritizing compounds for more rigorous computational and experimental study.