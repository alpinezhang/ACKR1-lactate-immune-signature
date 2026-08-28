# Reorganized analysis scripts (revision figure order)

R version used for syntax checking: **R 4.6.1**.

Before running the scripts, replace each `PATH_TO_DATA/...` placeholder with the directory containing the corresponding input data. For internal single-cell and spatial data, make the local sample-folder names agree with the anonymized `Sample_XXXX##` mapping (or edit that mapping to match your locally stored files). The scripts are provided in the figure order shown below; outputs are written to the subdirectories specified in each script.

Analysis code starts at Figure 2 because Figure 1 is a study workflow diagram. The graphical abstract is not an R-generated analysis figure.

| Current figure | Content | Script |
|---|---|---|
| Figure 2A-H | Internal single-cell overview, endothelial subtype validation, differential analysis, GO and KEGG | `01_Internal_SingleCell.R` |
| Figure 3A-J | Internal single-cell mechanism, metabolism, immune and clinical-association analyses | `01_Internal_SingleCell.R` |
| Figure 4A-G | GSE163154 differential expression, heatmap and core-gene/pathway expression | `03_External_Transcriptome_Proteome.R` |
| Figure 5A-F | GSE163154 correlation analyses | `03_External_Transcriptome_Proteome.R` |
| Figure 6A-J | GSE163154 correlation/HIF1A analyses plus internal proteomic validation and PPI | `03_External_Transcriptome_Proteome.R` |
| Figure 7A-F, H-L | LASSO, clinical model, ROC, DCA, risk score, nomogram, calibration and performance metrics | `03_External_Transcriptome_Proteome.R` |
| Figure 7G | Fully nested LOOCV validation | `05_Revision_Supplementary.R` |
| Figure 8A-J | GSE43292 and GSE28829 external validation | `03_External_Transcriptome_Proteome.R` |
| Figure 9 | Spatial transcriptomics colocalization and proximity analysis | `04_Spatial_Transcriptomics.R` |
| Supplementary Figure 1A-J | GSE253903 external single-cell validation | `02_External_SingleCell.R` |
| Supplementary Figure 2A-H | GSE163154 GO, KEGG, HALLMARK, glycolysis and hypoxia enrichment | `03_External_Transcriptome_Proteome.R` |
| Supplementary Figure 3A-J | GSE100927 carotid and femoral validation | `03_External_Transcriptome_Proteome.R` |
| Supplementary Figure 4A-E | GSE100927 infra-popliteal validation | `03_External_Transcriptome_Proteome.R` |
| Supplementary Figure 5A-F | ACKR1/glycolysis redraw and internal single-cell CellChat | `05_Revision_Supplementary.R` |
| Supplementary Figure 6A-H | Revision model stability, composition adjustment, CD8 proportion, key-gene and bulk-signature validation | `05_Revision_Supplementary.R` |

## Privacy redaction

Internal sample and patient identifiers were replaced with `XXXX`-style placeholders while preserving the original sample pairing and group structure. Spatial-transcriptomics sample labels were also replaced with anonymized placeholders. Local data paths were replaced with `PATH_TO_DATA` placeholders, so the input directories must be set locally before execution. Public GEO accession numbers (for example, GSE and GSM identifiers) were retained because they identify publicly available datasets required for reproducibility. The optional STRING network request uses a public endpoint and contains no embedded credential. No API tokens, database credentials, or server IP addresses are included.
