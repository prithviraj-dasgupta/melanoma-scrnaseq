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
