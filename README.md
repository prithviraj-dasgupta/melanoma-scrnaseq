# melanoma-scrnaseq

**Work in progress.** Separating malignant from non-malignant cells in human melanoma
single-cell RNA-seq (Tirosh et al. 2016, GEO [GSE72056](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056))
using Scanpy, marker-based cell-type annotation and inferred copy-number variation (infercnvpy).

The full README (question, methods, results with real numbers, limitations and how to
reproduce) will be written when the analysis is finished. Method choices, and the reasons for
them, are recorded as they are made in [docs/decisions.md](docs/decisions.md).

## Setup

```bash
uv sync
uv run pytest
```
