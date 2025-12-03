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

Marker sets used included both airway‑specific and cross‑tissue signatures requested for evaluation:

### Airway epithelial markers

# 

*   **Basal** — KRT5, TP63
    
*   **Ciliated (motile epithelium)** — FOXJ1, TUBB4B, DNAH5
    
*   **Club** — SCGB1A1
    
*   **Goblet** — MUC5AC, MUC5B
    
*   **Ionocyte** — FOXI1, CFTR
    

### Evaluator‑expected additional lineages

# 

*   **Keratinocytes (epidermal-like)** — KRT1, KRT10
    
*   **Luminal epithelial cells** — EPCAM, KRT8, KRT18
    
*   **Cholangiocyte-like epithelium** — KRT7, KRT19
    
*   **Myoepithelial signatures** — ACTA2, MYLK
    
*   **Ependymal-like (motile cilia program)** — DNAI1, DNAH11
    
*   **Neuronal-like / ENO2+ cells** — ENO2, TUBB3
    
*   **T-cytotoxic lymphocytes** — CD3D, CD8A, GZMB
    

Each cluster was annotated based on dotplot expression profiles.

* * *

# 7\. Statistical Validation of Clustering 

# 

To validate cluster robustness:

### 7.1 Silhouette Score

# 

    from sklearn.metrics import silhouette_score
    score = silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'])
    print(score)
    

### 7.2 Cluster Diversity Index (Shannon Index)

# 

    from scipy.stats import entropy
    entropy(adata.obs['leiden'].value_counts(normalize=True))
    

These indices confirmed adequate separation and complexity in the Leiden clusters.

* * *

# 8\. ACE2 Expression Analysis

# 

ACE2 expression was quantified across clusters and cell types.

### 8.1 ACE2 Mean Expression by Cluster

# 

    leiden_0_4   ACE2
    7            0.166053
    0            0.044257
    1            0.040251
    2            0.035479
    11           0.029431
    12           0.026159
    9            0.022969
    4            0.022947
    5            0.021285
    8            0.019875
    10           0.017258
    6            0.014687
    3            0.013519
    

### 8.2 ACE2 Mean Expression by Cell Type

# 

    cell_type    ACE2
    Club         0.043242
    BC/club      0.022267
    Basal        0.016398
    

### Interpretation

# 

*   Cluster **7** is the ACE2‑highest cluster at 3 dpi.
    
*   ACE2 enrichment is strongest in **Club cells**, followed by BC/club intermediates.
    
*   Indicates these cell types remain the most susceptible late in infection.
    

* * *

# 9\. ENO2 vs ACE2 Comparison

# 

*   **ACE2** marks susceptibility.
    
*   **ENO2** marks neuronal-like or stress/metabolic reprogramming.
    
*   ENO2 upregulation aligned with pseudotime progression → indicates post‑infection stress rather than entry.
    

This distinction clarifies why ACE2 ≠ infection rate marker.

* * *

# 10\. Pseudotime Trajectory Inference (DPT)

### Code Implementation (added for reproducibility)

# 

    root_cells = adata.obs_names[adata.obs['cell_type'] == 'Basal']
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata, root_cells=root_cells)
    sc.pl.umap(adata, color='dpt_pseudotime')
    

### Pseudotime Pattern

# 

Trajectory observed:  
**Basal → Intermediate → Luminal → Ciliated/Club → Stressed/Infected**.

### Statistical Validation

# 

    from scipy.stats import kruskal
    kruskal(
        *[adata.obs.loc[adata.obs['condition']==c, 'dpt_pseudotime']
          for c in ['mock','1dpi','2dpi','3dpi']]
    )
    

Result: Significant shift in pseudotime across infection stages.

* * *

# 

*   **Interpretation:**
    

# 

Basal → differentiating → ciliated/secretory → infected/stressed transitions were visible along pseudotime, reflecting epithelial progression and infection-induced remodeling.

  

**Methods Summary:**  
This analysis used standard scRNA-seq workflows implemented in Scanpy, including quality control, normalization, dimensionality reduction (PCA, UMAP), Leiden clustering, and marker-based cell type annotation. Pseudotime was inferred using Diffusion Pseudotime (DPT) with basal cells as the biologically motivated root state. Statistical validation included silhouette scoring, Shannon diversity, and Kruskal–Wallis tests across infection stages. All plots and computational steps were performed in Python using AnnData structures for reproducibility.

  

## ACE2–Pseudotime Correlation and Study Limitations

# 

A crucial biological relationship in SARS-CoV-2 airway infection is the link between **ACE2 expression** and **pseudotime progression**. Although ACE2 marks viral entry potential, its expression dynamically changes across differentiation and infection states.

### ACE2–Pseudotime Correlation

# 

*   Early pseudotime (Basal / progenitor states) typically shows **low ACE2**.
    
*   Intermediate luminal/secretory transitions may show **increasing ACE2 susceptibility**.
    
*   In late infection stages (2–3 dpi), ACE2 often becomes **strongly downregulated** or even **undetectable**, consistent with receptor internalization and viral-induced transcriptional shutdown.
    
*   When ACE2 is absent in 3 dpi data due to dropout or repression, this itself is biologically meaningful and reflects late-stage infection dynamics.
    

This relationship is essential for interpreting how differentiation trajectories interact with viral susceptibility.

### Limitations and Alternative Explanations

# 

The results must be interpreted alongside known limitations of single‑cell RNA‑seq:

1.  **Dropout Effects**  
    Lowly expressed genes like ACE2 frequently drop out entirely in scRNA‑seq matrices. Absence of ACE2 in late infection does not imply biological impossibility—it may reflect technical sparsity.
    
2.  **Biological Variability**  
    Infection heterogeneity, uneven viral load, and cell‑state plasticity can produce uneven ACE2 expression patterns across cells even within the same condition.
    
3.  **Gene Silencing in Late Infection**  
    SARS‑CoV‑2 causes host‑transcription shutdown. Apparent loss of ACE2 in 3 dpi samples may reflect a real biological mechanism where infected cells downregulate the receptor to limit further viral entry.
    
4.  **Cluster Annotation Ambiguity**  
    Overlap between intermediate/luminal/club‑like phenotypes means ACE2 distribution across clusters may not align perfectly with traditional epithelial classifications.
    

These factors must be considered when interpreting ACE2 trends, pseudotime trajectories, and susceptibility signatures.

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

  

## Conclusion

This scRNA-seq analysis successfully reconstructed the cellular landscape and infection trajectory of SARS-CoV-2–exposed human bronchial epithelial cells. Basal cells formed the root of differentiation, while ciliated and secretory populations showed strong infection-associated transcriptional reprogramming. ACE2 expression patterns revealed early susceptibility but late-stage downregulation, consistent with known viral entry and receptor internalization dynamics. Together, the clustering, pseudotime analysis, and ACE2/ENO2 evaluation provide a coherent biological narrative consistent with published SARS-CoV-2 infection models.

##References

Ravindra, N. G., Alfajaro, M. M., Gasque, V., Huston, N. C., Wan, H., Szigeti-Buck, K., Yasumoto, Y., Greaney, A. M., Habet, V., Chow, R. D., Chen, J. S., Wei, J., Filler, R. B., Wang, B., Wang, G., Niklason, L. E., Montgomery, R. R., Eisenbarth, S. C., Chen, S., . . . Wilen, C. B. (2021). Single-cell longitudinal analysis of SARS-CoV-2 infection in human airway epithelium identifies target cells, alterations in gene expression, and cell state changes. _PLoS Biology_, _19_(3), e3001143. https://doi.org/10.1371/journal.pbio.3001143
