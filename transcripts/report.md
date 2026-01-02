# Session 3_Analysis of transcriptomes (W dataset)

<details>
<summary><strong>🔴 Gene/isoform co-expression network construction and module identification in watered (W) RNA-seq samples using WGCNA</strong></summary>


## Goal
Construct a gene/isoform co-expression network using RNA-seq expression data from watered (W) samples and identify co-expression modules using WGCNA.

## Data and scripts
**Input files**
- `TPM_counts_Drought_W_dataset.csv` (expression matrix in TPM)
- `TRAITS_W.csv` (phenotypic traits for W samples)

**Script**
- `Practical_WGCNA_W_dataset.Rmd`

## Workflow summary (what was done)

### 1) Data import and formatting
The TPM expression table was loaded in R and converted to an expression matrix suitable for WGCNA:
- The first column contains isoform IDs (`target_id`).
- Columns correspond to samples.
- The dataset was transposed so that **rows = samples** and **columns = isoforms** (object `datExprW`).

### 2) Filtering genes/samples (quality check)
To remove problematic features, WGCNA’s quality-control function was applied:
- `goodSamplesGenes(datExprW)` checks missing values and removes genes/samples with excessive missingness.
- If any genes/samples failed QC, the matrix was subset to keep only valid ones.

### 3) Outlier detection and sample removal
Outlier samples were detected by hierarchical clustering:
- A sample dendrogram was built with `hclust(dist(datExprW), method="average")`.
- A cut height was applied (`cutHeight = 300000`) to define the main sample cluster.
- Only samples in the main cluster were retained (`keepSamplesW`), producing a filtered expression matrix saved as `datExpr_W.RData`.

### 4) Soft-threshold power selection
To build a scale-free-like network, a range of candidate powers was tested using:
- `pickSoftThreshold(datExprW, powerVector=powers, networkType="unsigned")`
The final power was chosen based on the scale-free topology fit criterion (R² threshold shown in the plot).

### 5) Network construction (adjacency and TOM)
Using the selected power:
- Adjacency was computed: `adjacency(datExprW, power = chosen_power)`
- Topological Overlap Matrix (TOM) was computed: `TOMsimilarity(adjacency)`
- A TOM-based dissimilarity was derived (`1 - TOM`) and used to cluster isoforms.

### 6) Module identification
Co-expression modules were detected using dynamic tree cutting:
- `cutreeDynamic(...)` assigned each isoform to a module.
- Modules were translated into colors (`labels2colors`) to label module membership.
A dendrogram with module colors was generated to visualize module structure.

### 7) Module eigengenes and module merging
For each module, an eigengene (ME) was calculated:
- `moduleEigengenes(datExprW, colors=modulecolors_W)`
Highly similar modules were merged:
- Eigengene dissimilarity was computed and clustered.
- Modules were merged using `mergeCloseModules(..., cutHeight = 0.25)`.
This produced the final merged module set and updated eigengenes.

### 8) Outputs produced
Key outputs for reporting and downstream analysis:
- Filtered expression matrix: `datExpr_W.RData`
- Network and modules: `net_W.RData`
- Per-module gene/isoform lists: `Genes_per_module_W.tsv`
- Module–trait association heatmap: `heatmap_traits_modules_W.pdf`

## References
- CIHEAM Zaragoza bioinformatics materials: https://eead-csic-compbio.github.io/bioinformatics  
- Practical WGCNA script: https://github.com/eead-csic-compbio/bioinformatics/blob/main/coexp/Practical_WGCNA_W_dataset.Rmd
- OpenAI ChatGPT – used for language refinement and clarification  
of BLAST output redirection during Session 3


</details>
