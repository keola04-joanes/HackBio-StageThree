# Final Report — Project 1: SARS-CoV-2 Infection Dynamics (scRNA-seq)

# 

This report summarizes the complete analysis workflow performed on the single-cell RNA sequencing dataset of SARS-CoV-2–infected human bronchial epithelial cells (GSE166766). The project includes data loading, quality control, integration, clustering, cell type annotation, and pseudotime trajectory inference.

* * *

# 1\. Data Acquisition & Organization

# 

The dataset consisted of four 10x Genomics MTX-format directories:

*   **mock/** (uninfected control)
    
*   **1dpi/** (1 day post-infection)
    
*   **2dpi/** (2 days post-infection)
    
*   **3dpi/** (3 days post-infection)
    

Each folder contained:

    matrix.mtx.gz
    barcodes.tsv.gz
    features.tsv.gz
    

These were loaded using `scanpy.read_10x_mtx()` with gene symbols as variable names.

Each dataset was labeled with a `condition` field: _mock_, _1dpi_, _2dpi_, _3dpi_.

* * *

# 2\. Quality Control (QC)

# 

To evaluate cell quality, the following QC metrics were computed:

*   **n\_genes\_by\_counts** (cell complexity)
    
*   **total\_counts** (library size)
    
*   **pct\_counts\_MT** (mitochondrial %)
    
*   **pct\_counts\_RIBO** (ribosomal %)
    
*   **pct\_counts\_HB** (hemoglobin genes)
    

Filtering criteria applied:

*   Cells with **<200 genes** removed.
    
*   Cells with **\>20% mitochondrial RNA** removed.
    

QC was performed separately for mock, 1dpi, 2dpi, and 3dpi before merging.

* * *

# 3\. Dataset Integration

# 

Filtered datasets were concatenated into a single AnnData object using `anndata.concat()`, with a `batch` label indicating sample origin.

* * *

# 4\. Normalization & Feature Selection

# 

Processing steps:

1.  **Total-count normalization** to 10,000 counts per cell
    
2.  **Log1p transformation**
    
3.  **Highly Variable Gene (HVG) selection** using Seurat v3 flavor (3000 HVGs)
    
4.  **Scaling** with regression of unwanted variation
    
5.  **Principal Component Analysis (PCA)**
    

* * *

# 5\. Nearest-Neighbor Graph, UMAP & Clustering

# 

*   Built with **20 neighbors** and **40 PCs**
    
*   UMAP computed for visualization
    
*   Leiden clustering performed with **resolution = 0.5**
    

The UMAP revealed distinct clusters corresponding to epithelial and immune cell identities.

* * *

# 6\. Cell Type Identification

# 

Using known marker genes:

*   **Basal**: KRT5, TP63
    
*   **Ciliated**: FOXJ1, TUBB4B
    
*   **Club**: SCGB1A1
    
*   **Goblet**: MUC5AC, MUC5B
    
*   **Ionocyte**: FOXI1, CFTR
    
*   **Immune**: PTPRC, LYZ
    

A dot plot was generated to inspect marker expression per cluster.  
Clusters were manually annotated and stored as `adata.obs['cell_type']`.

* * *

# 10\. Pseudotime Trajectory Inference (DPT)

# 

Pseudotime was computed using Diffusion Pseudotime (DPT):

1.  **Root selection**: Basal cells chosen as the starting point.
    
2.  **Diffusion map computed** using `sc.tl.diffmap()`.
    
3.  **DPT applied** with `sc.tl.dpt()`.
    
4.  UMAP colored by pseudotime to visualize lineage progression.
    

**Interpretation:**  
Basal → differentiating → ciliated/secretory → infected/stressed transitions were visible along pseudotime, reflecting epithelial progression and infection-induced remodeling.

* * *

#   

# 

*   Biological Interpretation & Results

Below are the biological answers and interpretations based on the scRNA-seq analysis of SARS-CoV-2 infected human bronchial epithelial cells across mock, 1 dpi, 2 dpi, and 3 dpi.

* * *

## 1\. Cell types identified at different stages of infection

Across all four samples, the following major airway epithelial and immune cell types were identified:

*   **Basal cells** – stem-like progenitor cells that maintain the epithelium.
    
*   **Ciliated cells** – fully differentiated epithelial cells responsible for mucociliary clearance.
    
*   **Secretory / Club cells** – protective secretory cells producing antimicrobial and surfactant proteins.
    
*   **Goblet cells** – mucus-producing epithelial cells.
    
*   **Ionocytes** – rare CFTR-high epithelial cells.
    
*   **Differentiating / Intermediate epithelial cells** – transitional states between basal and mature lineages.
    
*   **Immune cells** – including infiltrating macrophage-like or leukocyte populations.
    

**In the mock sample**, all major cell populations appear at baseline proportions with no infection-induced transcriptional signatures.

**At 1 dpi**, early infection responses begin to appear, with initial activation in ciliated and secretory cells.

**At 2 dpi**, interferon-stimulated genes and viral response pathways increase across infected clusters, and stress-like epithelial states become more apparent.

**At 3 dpi**, infected and damaged epithelial states dominate, with increased representation of stressed, inflamed, or dedifferentiated epithelial cells.

* * *

## 2\. Why these cell types correlate with COVID-19 infection

Ciliated and secretory/club cells show the strongest correlation with SARS-CoV-2 infection because:

*   They express **ACE2** and **TMPRSS2**, the viral entry receptor and protease required for infection.
    
*   They are **luminally exposed**, making them the first cells to contact inhaled virus.
    
*   Their **transcriptional and metabolic states** support viral entry and replication.
    
*   Infection of these cells impairs mucociliary clearance and disrupts epithelial barrier function.
    

Thus, cell-type vulnerability results from a combination of receptor availability, anatomical position, and biological function.

* * *

## 3\. Is ACE2 a good marker for tracking COVID-19 infection rate?

**No, ACE2 alone is not a reliable infection-rate marker** in single-cell RNA-seq.

Reasons:

*   ACE2 marks **susceptibility**, not active infection.
    
*   It is **lowly expressed** and prone to dropout in scRNA-seq.
    
*   ACE2 expression levels do **not correlate linearly** with viral load.
    
*   Viral transcripts or interferon-stimulated gene (ISG) signatures are **far better indicators** of infection intensity.
    

ACE2 should be interpreted as a receptor-associated vulnerability marker rather than a direct infection readout.

* * *

## 4\. Difference between ENO2 and ACE2 as biomarkers

*   **ACE2** is the **entry receptor** for SARS-CoV-2. It identifies which cells can be infected.
    
*   **ENO2** is a **metabolic/stress-related enzyme** not involved in viral entry. In infection datasets, ENO2 reflects **transcriptional reprogramming**, metabolic stress, or infected-cell state changes.
    

Thus, ACE2 indicates **which cells can be infected**, whereas ENO2 indicates **how cells respond after infection**.

* * *

## 5\. Which cell cluster has the highest ACE2 expression at 3 dpi & biological meaning

At 3 dpi, a single cluster shows the highest ACE2 abundance. This cluster:

*   Represents the **most ACE2-enriched epithelial population** at this late infection stage.
    
*   Indicates the population most permissive or still vulnerable to continued infection.
    
*   Shows localized ACE2 intensity on UMAP and clear expression differences in violin plots.
    

**Biological interpretation:**  
Cells in this ACE2-high cluster likely represent the epithelial subtype most actively interacting with the virus at this stage. High ACE2 expression suggests ongoing susceptibility, possible reinfection cycles, and significant epithelial remodeling or stress within this population.

This aligns with known COVID-19 pathology where specific epithelial subtypes (often ciliated or secretory) remain major viral targets as infection progresses.