# Analysis scripts for the ACKR1 study

This repository contains the analysis scripts for the manuscript and revision supplementary figures, including Supplementary Figures 7 and 8. Analyses start at Figure 2 because Figure 1 is a study workflow diagram. The graphical abstract is not an R-generated analysis figure.

**R version used for syntax checking:** R 4.6.1. 

## Figure-to-script map

| Figure | Analysis | Script |
|---|---|---|
| Figure 2A-H | Internal single-cell overview, endothelial subtype validation, differential analysis, GO and KEGG | `01_Internal_SingleCell.R` |
| Figure 3A-J | Internal single-cell mechanism, metabolism, immune and clinical-association analyses | `01_Internal_SingleCell.R` |
| Figure 4A-G | GSE163154 differential expression, heatmap and core-gene/pathway expression | `03_External_Transcriptome_Proteome.R` |
| Figure 5A-F | GSE163154 correlation analyses | `03_External_Transcriptome_Proteome.R` |
| Figure 6A-J | GSE163154 correlation/HIF1A analyses, internal proteomic validation and PPI | `03_External_Transcriptome_Proteome.R` |
| Figure 7A-F, H-L | LASSO, clinical model, ROC, DCA, risk score, nomogram, calibration and performance metrics | `03_External_Transcriptome_Proteome.R` |
| Figure 7G | Fully nested LOOCV validation | `05_Supplementary_Analysis.R` |
| Figure 8A-J | GSE43292 and GSE28829 external validation | `03_External_Transcriptome_Proteome.R` |
| Figure 9 | Spatial transcriptomics colocalization and proximity analysis | `04_Spatial_Transcriptomics.R` |
| Supplementary Figure 1A-J | GSE253903 external single-cell validation | `02_External_SingleCell.R` |
| Supplementary Figure 2A-H | GSE163154 GO, KEGG, HALLMARK, glycolysis and hypoxia enrichment | `03_External_Transcriptome_Proteome.R` |
| Supplementary Figure 3A-J | GSE100927 carotid and femoral validation | `03_External_Transcriptome_Proteome.R` |
| Supplementary Figure 4A-E | GSE100927 infra-popliteal validation | `03_External_Transcriptome_Proteome.R` |
| Supplementary Figure 5A-F | ACKR1/glycolysis redraw and internal single-cell CellChat | `05_Supplementary_Analysis.R` |
| Supplementary Figure 6A-H | Revision model stability, composition adjustment, CD8 proportion, key-gene and bulk-signature validation | `05_Supplementary_Analysis.R` |
| Supplementary Figure 7A-C | GSE43292 within-patient paired comparisons and conditional logistic regression | `05_Supplementary_Analysis.R` |
| Supplementary Figure 7D-E | Internal proteomic key proteins and KEGG pathways before and after clinical adjustment | `05_Supplementary_Analysis.R` |
| Supplementary Figure 7F | GSE163154 nested LOOCV calibration and bootstrap confidence intervals | `05_Supplementary_Analysis.R` |
| Supplementary Figure 7G | Proteome-wide effect sizes before and after clinical adjustment | `05_Supplementary_Analysis.R` |
| Supplementary Figure 8A-B | GSE163154 model features and CD8 cytotoxicity-score genes | `05_Supplementary_Analysis.R` |

## Running the analyses

Replace each `PATH_TO_DATA/...` placeholder with the local path to the corresponding input data, and align local sample labels with the `XXXX`-style placeholders while preserving pairing and group membership. Install the packages required by each script and run sections in their listed order; later sections may depend on earlier outputs. Record the R and package versions using `sessionInfo()`.

Internal single-cell, spatial, proteomic and clinical inputs are not included. Public GEO accession numbers are retained for reproducibility. Chinese input filenames, worksheet names and column names are retained as data-reading keys and must match the local input files. Some sections install packages, download public data or query STRING/KEGG; results depend on the input files, package and database versions, and service availability.

## Privacy and responsible use

Internal sample and patient labels have been replaced with `XXXX`-style placeholders, and local data paths with `PATH_TO_DATA/...`. These replacements do not anonymize runtime inputs or outputs. Users are responsible for obtaining required data-access and ethics approvals and reviewing all outputs before sharing. Do not publish patient records, identifier mappings, credentials, private datasets, unredacted logs or R session history, including through GitHub issues or pull requests. 

These scripts support research and peer review and have not been validated for clinical use. Cite the associated study and the datasets and software used, and record the repository commit or release. Code licensing is governed by the repository's `LICENSE` file; third-party software and datasets remain subject to their respective terms. A code license does not authorize access to or redistribution of patient data.
