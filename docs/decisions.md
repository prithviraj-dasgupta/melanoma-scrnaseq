# Decisions

Method and scope decisions for melanoma-scrnaseq, with the reason for each. Each entry says
where it is (or will be) implemented. Entries marked *provisional* are re-checked once the data
file has been inspected (Phase 1).

## 2026-11-24: Plan

### D1. Dataset: Tirosh et al. 2016 melanoma (GEO GSE72056)
- Human metastatic melanoma, Smart-seq2 (full-length, plate-based) single-cell RNA-seq,
  about 4,600 cells from 19 patients.
- Chosen because every cell carries **author labels**: malignant / non-malignant /
  unresolved, plus a non-malignant cell type (T, B, macrophage, endothelial, CAF, NK). That
  lets us check two of our own steps (marker annotation and CNV-based malignancy calls)
  against an outside reference. It is one file and small enough to run on a laptop.
- Alternative considered: GSE148673 (10x Genomics, triple-negative breast cancer, raw UMI
  counts). Rejected because its tumour labels come from another CNV tool (copyKAT), so
  validation would be weaker, though it would have shown standard raw-count QC.
- *Where:* `scripts/01_download.py`, `src/melscrna/load.py` (Phase 1).

### D2. Engineering standards
- Same standards as the `cancer-fingerprints` repo: `uv` with a lock file, `ruff` lint and
  format, `pytest`, GitHub Actions CI on Python 3.11 and 3.12.
- Raw data is downloaded by a script and never committed. Tests run on tiny synthetic data
  so CI needs no download; tests that need the internet are marked `network` and skipped
  in CI.
- *Where:* `pyproject.toml`, `.github/workflows/ci.yml`.

### D3. The data are log2(TPM/10 + 1), not raw counts (confirmed in D6)
- GSE72056 is shared as processed expression values, not read or UMI counts.
- Consequences:
  - **No count-depth normalisation** (`normalize_total` + `log1p`). The values are already
    length- and depth-normalised and log-transformed.
  - **No scrublet doublet detection.** Scrublet simulates doublets by adding raw count
    profiles together, which is not valid on log TPM. Smart-seq2 cells are sorted one per
    well by FACS, so doublets are rare by design. This is stated as a limitation rather
    than hidden.
  - QC uses genes detected per cell. Mitochondrial share cannot be used: the file has no
    `MT-` genes (D6). The authors deposited only cells that already passed
    their own QC, so few cells are expected to be removed.
- *Where:* `src/melscrna/qc.py` (Phase 2).

### D4. Agreement with author malignancy labels is not independent ground truth
- Tirosh et al. defined malignant cells partly from inferred CNV. Agreement between our
  infercnvpy calls and their labels shows that two CNV-based approaches agree; it does not
  prove either is correct.
- *Where:* README limitations, Phase 5.

### D5. Biological question: MITF-high vs AXL-high malignant states
- The non-malignant cells are immune and stromal cells, not normal melanocytes. A
  malignant-vs-non-malignant differential expression test would mostly rediscover
  lineage genes (melanocyte vs immune), which is trivial.
- Instead we ask a modest replication question: do malignant cells show MITF-high and
  AXL-high programmes, are they anti-correlated, and how does the mix vary between patients?
  AXL-high states have been linked to resistance to BRAF/MEK inhibitors.
- *Where:* Phase 6.

## 2026-11-24: Phase 1, what the file actually contains

### D6. Data check before any analysis
Run with `scripts/02_inspect_data.py` on the downloaded file
(`GSE72056_melanoma_single_cell_revised_v2.txt.gz`, 75.0 MB,
MD5 `9c05cb22103d01d3086a2a952e97e96b`).
- **4,645 cells x 23,686 genes, 19 tumours.** Two gene symbols appear twice; they are kept
  and made unique (`GENE-1`) rather than silently dropped.
- **The values are log2(TPM/10 + 1), as stated.** None are whole numbers, so they are not
  counts. Undoing the transform and summing each cell gives a median of 999,999 TPM (range
  951,669 to 1,000,160). TPM sums to one million per cell by definition, so this confirms
  the scale. This confirms D3: no count normalisation and no scrublet.
- **No mitochondrial (`MT-`) genes are in the file.** The authors removed them before
  sharing, so the usual "% mitochondrial reads" QC (a sign of dying or broken cells) is not
  possible. This is a limitation, reported in the README.
- **Genes detected per cell:** median 4,058, minimum 1,370, maximum 13,276. Every cell
  already has many detected genes, consistent with the authors having filtered low-quality
  cells before deposit.
- **Author labels:** 1,257 malignant, 3,256 non-malignant, 132 unresolved. Of the
  non-malignant cells, 2,887 have a cell type (T 2,068; B 515; macrophage 126; endothelial
  65; CAF 61; NK 52) and 369 have none. Unresolved cells are kept through clustering and
  left out of label-agreement statistics.
- *Where:* `src/melscrna/load.py`, `scripts/02_inspect_data.py`,
  `results/tables/cells_per_tumour.tsv`.

## 2026-11-24: Phase 2, quality control

Results in `results/tables/qc_summary.tsv` and `results/figures/qc_distributions.png`.
4,645 cells in, **4,619 cells x 22,285 genes out.**

### D7. Cell filter on genes detected: median +/- 5 MADs (log scale) -> 0 cells removed
- A cell with very few detected genes is usually an empty or broken well. One with very many
  can be a doublet. We use a data-driven cut-off (median +/- 5 median absolute deviations,
  on log1p values) rather than a fixed number, because "low" depends on the technology and
  the tissue.
- Bounds were 1,068 to 15,405 genes. No cell falls outside them, which is expected because
  the authors filtered cells before deposit (D6).
- **Why 5 MADs and not the stricter 3:** malignant cells detect more genes than immune cells
  (median 5,523 vs 3,626), because they are larger and transcriptionally more active. A
  3-MAD cut-off (1,822 to 9,035) would remove 75 cells, 31 of them malignant, mostly from
  the high end. That is real biology, not poor quality, so a tight upper bound would bias
  the dataset against tumour cells.
- *Where:* `src/melscrna/qc.py::run_qc`, `scripts/03_qc.py`.

### D8. Doublet check from marker co-expression -> 26 cells removed
- Scrublet cannot be used (D3). Instead we use a biological rule: no single cell is both an
  immune cell (PTPRC/CD45) and a melanocyte-lineage cell (PMEL, MLANA or TYR). Cells with
  both above 7 log2(TPM/10+1) (about 1,270 TPM, i.e. strong expression of both) are removed.
- 26 cells were removed. The authors had labelled all of them non-malignant (25 T cells,
  1 NK cell).
- We cannot tell whether these are true doublets (a T cell stuck to a melanoma cell) or
  heavy contamination with tumour RNA. Either way their profile is mixed, so they are removed.
- **Reported, not filtered:** at a lower level (any marker > 4, about 150 TPM), 308 of 2,068
  T cells contain melanocyte transcripts. This graded, low-level signal looks like ambient
  or contaminating tumour RNA rather than doublets, so those cells are kept. It is worth
  remembering when interpreting T-cell expression.
- The threshold of 7 is a judgement call and is stated as such.
- *Where:* `src/melscrna/qc.py::flag_lineage_conflict`.

### D9. Gene filter: detected in at least 3 cells -> 1,401 genes removed
- Genes seen in fewer than 3 cells cannot help define cell populations and add noise.
  840 genes were never detected at all.
- *Where:* `src/melscrna/qc.py::run_qc`.

### D10. Ribosomal share reported, not filtered
- Ribosomal protein genes make up a median of 17% of TPM (13% in malignant cells, 18% in
  non-malignant). There is no agreed cut-off for this on TPM data and it reflects cell
  type, so it is plotted but not used to remove cells.

## 2026-11-24: Phase 3, embedding and clustering

Results in `results/tables/cluster_composition.tsv`, `results/figures/umap_overview.png` and
`results/figures/pca_elbow.png`.

### D11. 2,000 highly variable genes, Seurat (dispersion) method
- Most genes vary little between cells and only add noise. We keep the 2,000 genes whose
  variance is highest relative to their mean.
- The `seurat` flavour works on log-normalised values. The `seurat_v3` flavour, often used
  for 10x data, models raw counts and so is not valid here (D3).
- *Where:* `src/melscrna/cluster.py::embed_and_cluster`.

### D12. 30 principal components, 15 neighbours, Leiden resolution 0.5 -> 20 clusters
- Genes are scaled (mean 0, variance 1, clipped at 10) before PCA, so highly expressed genes
  do not dominate. PC1 explains 4.2% of variance and the curve flattens after about 20 PCs
  (`pca_elbow.png`); we use 30 to be safe, which together explain 22.1%.
- Leiden finds communities in the 15-nearest-neighbour graph. Resolution controls how many:
  0.5 gives 20 clusters, which is enough to separate the main cell types without splitting
  them into many small groups that are hard to label. Random seeds are fixed (0).
- UMAP is used only to draw the data. Distances between far-apart UMAP islands are not
  meaningful, so no conclusions are drawn from them.
- *Where:* `src/melscrna/cluster.py`, `scripts/04_cluster.py`.

### D13. No batch correction across patients
- Malignant cells cluster by patient: in the 7 clusters that are at least 80% malignant, a
  median of 99% of cells come from a single tumour. In the 12 clusters that are at most 20%
  malignant, the median is only 28%; immune and stromal cells from different patients mix.
- This matches Tirosh et al.: each patient's tumour has its own copy-number changes and
  expression programme, while normal cell types look alike across people. If it were a
  technical batch effect, non-malignant cells from the same patient would also group
  together, and they do not.
- Batch correction (e.g. Harmony) would remove this real biology, so it is not used.
- One cluster (9) is mixed: 55 malignant and 54 non-malignant cells, mostly from Mel59.
  It is examined during annotation (Phase 4).
- *Where:* `scripts/04_cluster.py`, README results.

## 2026-11-24: Phase 4, marker-based annotation

Results in `results/tables/cluster_annotation.tsv` (label, scores and top genes per cluster),
`results/tables/annotation_vs_author.tsv`, `results/figures/marker_dotplot.png` and
`results/figures/umap_annotation.png`.

### D14. Annotate clusters, not single cells, from short canonical marker lists
- Single-cell data are sparse: a gene can read zero in a cell that expresses it. Averaging
  over a cluster is more robust than labelling each cell, at the cost of missing small
  populations inside a cluster (see D16).
- Each cell type has 3 to 5 textbook markers (`src/melscrna/annotate.py::MARKERS`), e.g.
  CD3D/CD3E/CD2 for T cells, CD19/MS4A1/CD79A for B cells, PECAM1/VWF/CDH5 for endothelium,
  COL1A1/COL1A2/COL3A1/DCN for fibroblasts, and PMEL/MLANA/TYR/MITF/SOX10 for the melanocyte
  lineage.
- Rules: the top-scoring type wins; if no list averages 2.5 or more the cluster is
  `Unclear`; if a second, unrelated lineage scores at least 2 and at least half the top
  score, it is `Mixed`. Plasma/B and NK/T pairs share genes for biological reasons and are
  not called mixed.
- Every label is checked against the cluster's own top genes from a Wilcoxon test (one
  cluster vs the rest), which do not depend on our marker choice. For example, cluster 1's
  top genes are CCL5, NKG7, CD8A and PRF1 (cytotoxic CD8 T cells), and cluster 7's are
  TYROBP, AIF1, LYZ and CD68 (macrophages).
- **"Melanocytic" is not "malignant".** Normal melanocytes express the same genes, and some
  melanoma cells switch them off (D16). The malignant call is made from copy number in
  Phase 5.

### D15. Marker fixes found by looking at the data
- The first plasma-cell markers (MZB1, IGJ/JCHAIN, DERL3) were also high in cluster 10,
  whose top genes (IL3RA, PLD4, IRF7, IRF8) mark **plasmacytoid dendritic cells (pDCs)**,
  which also express those genes. Plasma cells are now marked by SDC1 (CD138), PRDM1
  (BLIMP1) and TNFRSF17 (BCMA), and a pDC list (LILRA4, IL3RA, CLEC4C, PLD4) was added.
  This gives 19 pDCs and 25 plasma cells, neither of which is a category in the author
  labels.

### D16. Result: 93.7% agreement with author cell types, with honest gaps
- Of 2,861 cells with an author cell type, our cluster labels agree for **2,681 (93.7%)**:
  T cells 98.1%, B cells 90.9% (456 in the B-cell cluster and 12 in the plasma-cell
  cluster; most of the other 47 sit in T-cell clusters 13 and 8), macrophages 100%,
  endothelial 100%.
- **NK cells: 0 of 51.** They sit inside the large cytotoxic T-cell cluster (cluster 1), and
  a cluster-level method cannot pull out 51 cells from 1,045. Sub-clustering could, but is
  out of scope; this is reported as a limitation.
- **CAFs: 18 of 61.** The other 43 are in cluster 9, which is labelled `Mixed: CAF /
  Melanocytic` because it combines 54 non-malignant cells (mostly CAFs) and 55 malignant
  cells, mostly from Mel59.
- **165 author-malignant cells are in `Unclear` clusters (3 and 15).** Cluster 3 (Mel81)
  has low MITF and SOX10 and expresses collagens and FN1; cluster 15 (Mel82) keeps SOX10
  and NCAM1 but loses MITF and PMEL. These look like de-differentiated, MITF-low melanoma
  cells, which marker-based annotation cannot call. This motivates the copy-number step
  (Phase 5) and the MITF/AXL question (Phase 6).
- Overall, 1,078 of 1,257 author-malignant cells (85.8%) are in Melanocytic or mixed
  Melanocytic clusters.

## 2026-11-24: Phase 5, copy-number inference and malignant calls

Results in `results/tables/cnv_vs_author.tsv`, `cnv_sensitivity_per_tumour.tsv`,
`cnv_by_cell_type.tsv`, `cnv_sensitivity_analysis.tsv`, and `results/figures/cnv_heatmap.png`
and `cnv_signal_vs_correlation.png`.

### D17. Gene positions from GENCODE v19 (hg19), parsed by our own code
- infercnvpy needs each gene's chromosome and position. GENCODE v19 matches the genome
  build of the 2016 data. 20,652 of 22,285 genes (92.7%) were placed; the rest (mostly
  old or non-standard symbols) are skipped.
- chrY is left out, so the CNV profile cannot pick up patient sex.
- We parse the GTF ourselves (`read_gene_positions`, ~30 lines, tested) instead of adding
  the `gtfparse` dependency.

### D18. Reference = half of the immune and stromal cells, one baseline per cell type
- The reference comes from **our** marker-based labels (Phase 4), not the author labels,
  so the comparison with the authors is not fed their answer.
- Cells in T, B, plasma, macrophage, pDC, endothelial and CAF clusters are split 50/50
  within each type (seed 0): 1,610 reference cells and 1,610 held-out normal cells. The
  1,399 cells in Melanocytic, Mixed and Unclear clusters are the query.
- Each normal cell type is its own reference category. With a single pooled reference
  (mostly T cells), genes that fibroblasts express from neighbouring loci look like copy
  gains, and 17 of 32 CAFs were called malignant. With per-type references, 0 are.
- *Where:* `src/melscrna/cnv.py::split_reference`, `scripts/06_cnv.py`.

### D19. HLA genes excluded
- The MHC region on chr6p21 is a dense cluster of genes expressed far more in immune cells.
  Against an immune reference it showed up as a "deletion" in every non-immune cell (a
  vertical stripe on chr6 in the first heatmap). The 25 HLA genes are removed before
  inference. The pattern also covers immunoglobulin and T-cell receptor segment genes, but
  none are in this file.

### D20. Malignant call: two CNV measures, thresholds from held-out normal cells
- Following Tirosh et al., each cell gets a **CNV signal** (mean squared CNV value) and a
  **CNV correlation** with its own tumour's consensus profile (mean of the 10% of that
  tumour's cells with the highest signal, chosen without labels).
- A cell is malignant if either measure is above the 99th percentile of the **held-out**
  normal cells (signal > 0.0102 or correlation > 0.563). The reference cells cannot be used
  for this because their profiles are flat by construction, which would set the bar too
  low.
- 23 of 1,610 held-out normal cells (1.4%) pass, close to the 1-2% expected from two 99th
  percentile cut-offs.
- Window: the infercnvpy default of 100 genes.

### D21. Result: 91.5% agreement, high specificity, moderate sensitivity
- Against the author labels (132 unresolved cells excluded): **agreement 91.5%,
  sensitivity 73.3% (922 of 1,257), specificity 98.6% (3,184 of 3,230).**
- Sensitivity varies by tumour: at least 85% in 7 of the 15 tumours with author-malignant
  cells, but 60% in Mel79 (279 of 468) and 0 of 3 in Mel75. The heatmap shows that Mel79's
  missed cells do carry its CNV pattern (losses on chr10, 13 and 19), just more weakly
  than the threshold. The cut-off is deliberately strict.
- The Unclear clusters from Phase 4 are mostly malignant by copy number: 107 of 155 cells
  in cluster 3 and 18 of 36 in cluster 15. So the MITF-low cells that markers could not
  place are tumour cells.
- **Sensitivity analysis** (`cnv_sensitivity_analysis.tsv`): a pooled reference raises
  sensitivity to 89.8% (100-gene window) or 94.6% (150), but calls 26 to 27 CAF,
  endothelial or macrophage cells malignant. A 150-gene window with per-type references
  gives 82.3%. We keep the stricter main setting because the next question (AXL-high
  malignant cells) would be biased by fibroblasts called malignant, since fibroblasts
  express AXL. This choice was made for that reason, not to maximise agreement.
- Reminder (D4): the author labels were themselves partly CNV-based, so agreement means
  two CNV methods agree, not that either is ground truth.

## 2026-11-24: Phase 6, MITF-high vs AXL-high malignant states

Results in `results/tables/programme_genes.tsv`, `programme_scores.tsv`,
`programmes_per_tumour.tsv` and `results/figures/mitf_axl_programmes.png`.

### D22. Which cells: CNV-called malignant cells outside immune/stromal clusters
- 941 cells from 14 tumours. Using our CNV calls (not the author labels) keeps the analysis
  self-contained. Excluding the few immune or stromal cells that passed the CNV threshold
  stops, for example, fibroblasts (which express AXL) from inflating the AXL state.
- Genes detected in more than 10% of these cells are used (12,567).

### D23. Programmes learned on one half of the cells, tested on the other
- A single gene is noisy in one cell, so each state is a programme: the 50 genes most
  correlated with MITF or with AXL, computed **within tumours** (each tumour's mean
  removed), so that genes are not picked just because some tumours are high in both.
- Genes are chosen on a random half of the malignant cells in each tumour (seed 0). The
  correlation tests use only the other half. Choosing genes and testing them on the same
  cells would build in part of the answer.
- The programmes look biologically sensible without being told what to find. The MITF
  programme includes MLANA, PMEL, TYR, GPR143, S100B and ERBB3 (melanocyte
  differentiation). The AXL programme includes NGFR, TNC, LIF and EMP1 (genes linked to
  the invasive, undifferentiated melanoma state). The two share no genes.
- A cell's programme score is the mean z-score of the programme genes. "AXL-dominant"
  means AXL score above MITF score. Both are relative to the average malignant cell in
  this dataset, not an absolute biological state.

### D24. Result: the states differ mainly between patients; the within-tumour link is weak
- **Within tumours, the MITF and AXL programmes are only weakly anti-correlated:** Spearman
  rho = -0.10 (p = 0.04) across 473 held-out cells. On its own, only Mel79 (rho -0.34,
  140 cells) and Mel78 (-0.30, 57 cells) show a clear negative correlation. The other 6
  tumours that could be tested show no significant correlation in either direction.
- **The big differences are between tumours.** In Mel81, Mel82 and Mel53, 99-100% of
  malignant cells are AXL-dominant. In Mel78, Mel88 and Mel80 only 7-22% are. Mel81 and
  Mel82 are the de-differentiated cells that marker-based annotation could not label
  (clusters 3 and 15).
- Across the 11 tumours with at least 10 malignant cells, median MITF and AXL scores are
  negatively but not significantly correlated (rho -0.30, p = 0.37). Eleven tumours is too
  few to say more.
- **What this does and does not show.** It is consistent with Tirosh et al. that MITF-high
  and AXL-high states exist and that tumours differ in their mix. It gives only weak
  support for a strong anti-correlation within each tumour, and with Smart-seq2 data from
  up to 279 malignant cells per tumour we cannot say whether that is biology or limited
  power. There is no treatment information here, so nothing can be said about drug
  resistance itself.

### Note on reproducibility across machines
- The numbers in this log come from runs on macOS (Apple Silicon). Re-running on Linux gave
  the same QC, clusters and headline CNV metrics, but a few cells near a threshold can
  flip (e.g. 278 vs 279 Mel79 cells called malignant), because floating-point results
  differ slightly between platforms. UMAP layouts also look rotated between machines.
