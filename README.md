# SC_ST_CTB_figures_scripts
A Spatiotemporal Single-Cell Atlas of Cutaneous Tuberculosis Infection and Treatment Response

## Repository Contents
| Directory | Description |
| :--- | :--- |
| **`deconvolution_RCTD/`** | Cell-type deconvolution of spatial transcriptomics data using **RCTD (Robust Cell Type Decomposition)**. |
| **`decoupleR_score/`** | Calculation of pathway activity scores at single-cell/spot resolution using the **decoupleR** tool. |
| **`pyscenic/`** | Inference of gene regulatory networks and transcription factor activity via the **pySCENIC** workflow. |
| **`single_cell_cellcrosstalk/`** | Analysis of cell-cell communication from single-cell data. |
| **`spatial_distance_calculation/`** | Quantification of spatial distances between different cell types or features within tissue sections. |
| **`spatial_niche_analysis/`** | Identification and characterization of spatial cellular niches (microenvironments). |
| **`spatial_NMF_analysis/`** | Dimensionality reduction and pattern discovery in spatial data using **Non-negative Matrix Factorization (NMF)**. |

## Data Availability
The raw and processed sequencing data generated in this study have been deposited in the **Gene Expression Omnibus (GEO)** under accession number **[GSExxxxxx](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSExxxxxx)** (Link is a placeholder; replace `GSExxxxxx` with the actual number).

*Please replace `GSExxxxxx` with the actual accession number once available.*

## Getting Started
*   **Dependencies**: Most analyses are implemented in R or Python. Key packages are listed within each module's directory (e.g., in a `requirements.txt` or `environment.yml` file).
*   **Usage**: Navigate to the relevant directory for a specific analysis. Each typically contains a main script (e.g., `run_analysis.R`) and a dedicated README with detailed instructions.
*   **General Workflow**: The analysis is modular. You can start with deconvolution (`deconvolution_RCTD/`) to obtain cell type proportions, then proceed to downstream analyses like niche identification or cell-cell communication.

## Citation
If you use the code or data from this repository, please cite our accompanying publication:
> *[Publication Title].* Authors. Journal (Year). DOI: [DOI-Link]

*(Citation placeholder; update upon publication)*
