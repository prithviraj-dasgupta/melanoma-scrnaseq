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

### D3. The data are log2(TPM/10 + 1), not raw counts (*provisional*)
- GSE72056 is shared as processed expression values, not read or UMI counts.
- Consequences:
  - **No count-depth normalisation** (`normalize_total` + `log1p`). The values are already
    length- and depth-normalised and log-transformed.
  - **No scrublet doublet detection.** Scrublet simulates doublets by adding raw count
    profiles together, which is not valid on log TPM. Smart-seq2 cells are sorted one per
    well by FACS, so doublets are rare by design. This is stated as a limitation rather
    than hidden.
  - QC uses genes detected per cell and, if mitochondrial genes are present in the file,
    mitochondrial expression share. The authors deposited only cells that already passed
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
