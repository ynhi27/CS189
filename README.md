# 🧬 Systematic Benchmarking of Differential Gene Expression Methods for ER Subtype Classification in Breast Cancer with Literature-Aware Candidate Prioritization

<div align="center">

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?style=flat-square&logo=python&logoColor=white)
![R](https://img.shields.io/badge/R-4.3%2B-276DC3?style=flat-square&logo=r&logoColor=white)
![Platform](https://img.shields.io/badge/Platform-Google%20Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white)
![Dataset](https://img.shields.io/badge/Dataset-TCGA--BRCA-00897B?style=flat-square)
![Status](https://img.shields.io/badge/Status-Complete-27AE60?style=flat-square)

**Systematic benchmarking of DESeq2, edgeR, and limma-voom for estrogen receptor subtype classification in breast cancer, with automated literature-aware candidate prioritization.**

[Overview](#-overview) · [Results](#-key-results) · [Pipeline](#-pipeline) · [Quick Start](#-quick-start) · [Outputs](#-output-files) · [Team](#-team)

</div>

---

## Overview

This project benchmarks three widely used RNA-seq differential expression methods — **DESeq2**, **edgeR**, and **limma-voom** — on the TCGA-BRCA cohort for ER+/ER− subtype classification.

Rather than the standard direct ER+ vs ER− comparison, we use a **three-group design** (each subtype independently vs. matched Normal tissue) and the **Li et al. four-class dysregulation framework** to capture the full landscape of subtype-specific gene expression changes.

A **novel contribution** of this work is an automated PubMed literature confidence scoring system, augmented by an **Observe–Decide–Act agentic refinement loop** using the Anthropic Claude API, which identifies reproducible but understudied gene candidates.

---

## Key Results

| Metric | Value |
|---|---|
| Samples | 100 ER+ · 100 ER− · ≤100 Normal (TCGA-BRCA) |
| Genes after filtering | **18,627** |
| DESeq2 DEGs | **11,945** |
| edgeR DEGs | **9,120** |
| limma-voom DEGs | 0 ⚠️ (edgeR v4 API issue — see [Limitations](#⚠️-known-limitations)) |
| Jaccard similarity (DESeq2 vs edgeR) | **0.38** |
| log2FC Pearson r (all pairs) | **> 0.80** |
| DESeq2 CV (robustness) | **1.92%** |
| DESeq2 Jaccard stability | **0.876 ± 0.012** |

> All known ER markers (ESR1, GATA3, TFF1, XBP1, CDH2) recovered in both methods.  
> **Primary finding:** Class 2 Tier 3 consensus genes — statistically reproducible, not yet studied in ER context.

---

## Pipeline

```
TCGA-BRCA (GDC API)
        │
        ▼
 Sample Selection          100 ER+ · 100 ER− · N Normal  (seed = 42)
        │
        ▼
  CPM Filtering            CPM ≥ 1 in ≥ 10% samples → 18,627 genes
        │
        ├──────────────────────────────────────┐
        ▼                                      ▼
   DESeq2                                   edgeR
  (pydeseq2)                            (rpy2 + R)
  ER+ vs Normal                        ER+ vs Normal
  ER− vs Normal                        ER− vs Normal
        │                                      │
        └───────────────┬──────────────────────┘
                        ▼
              Consensus Gene Set
           (DESeq2 ∩ edgeR fallback)
                        │
          ┌─────────────┼──────────────┐
          ▼             ▼              ▼
    Concordance    Four-Class     Robustness
    Analysis      Classification  (80% × 10)
   Jaccard, r    Class 1/2/3/4    CV, Jaccard
          │             │              │
          └─────────────┼──────────────┘
                        ▼
             Biological Validation
            (Markers + KEGG Enrichment)
                        │
                        ▼
          Literature Confidence Scoring
          PubMed + MyGene.info aliases
          Tier 1 (ER-validated)
          Tier 2 (BC literature)
          Tier 3 (understudied)
                        │
                        ▼
           Agentic Refinement Loop
         Observe → Decide → Act
       (Anthropic Claude API + alias expansion)
                        │
                        ▼
         final_ranked_candidates.csv
```

---

## Quick Start

### Prerequisites

This notebook runs on **Google Colab** with a Google Drive mount for caching.

### 1. Clone the repository

```bash
git clone https://github.com/ynhi27/CS189.git
cd <your-repo>
```

### 2. Open in Google Colab

Upload `Group8_Final_project_v4.ipynb` to Google Colab, or click the badge below:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/<your-username>/<your-repo>/blob/main/Group8_Final_project_v4.ipynb)

### 3. Install dependencies

Run **Cell 3** in the notebook. It installs all required packages automatically:

```python
# Python packages (installed via pip)
pydeseq2, requests, matplotlib, seaborn, scipy,
statsmodels, matplotlib-venn, gseapy, anthropic

# R packages (installed via rpy2/BiocManager)
edgeR, limma
```

### 4. Set up API keys

The agentic refinement loop (Section 3.9) requires an Anthropic API key.
Add it to your Colab secrets or set it as an environment variable:

```python
import os
os.environ["ANTHROPIC_API_KEY"] = "your-api-key-here"
```

> The NCBI Entrez API is used for PubMed queries. No key is required for low-volume use,
> but setting `Entrez.email` in the notebook is recommended.

### 5. Run all cells

Use **Runtime → Run all**. The notebook includes auto-caching:
if intermediate outputs already exist in Google Drive, they are loaded automatically.

---

## Repository Structure

```
.
├── Group8_Final_project_v4.ipynb   # Main analysis notebook
├── README.md                       # This file
└── TCGA_BRCA_Project/              # Generated after running notebook (output folder)
    ├── sample_metadata.csv         # Sample IDs, ER status, condition labels
    ├── count_filtered.parquet      # Filtered count matrix (18,627 × N samples)
    ├── gene_class.csv              # Four-class dysregulation labels per gene
    ├── gene_literature_scores.csv  # PubMed tier assignments
    └── final_ranked_candidates.csv # Ranked Class 2 Tier 3 novel candidates
```

All outputs are cached to **Google Drive** (`/content/drive/MyDrive/TCGA_BRCA_Project/`).
Re-running the notebook after a cache hit skips expensive API calls automatically.

---

## Output Files

| File | Description |
|---|---|
| `sample_metadata.csv` | Sample IDs, ER status (IHC), condition label (ER_pos / ER_neg / Normal) |
| `count_filtered.parquet` | CPM-filtered raw count matrix |
| `gene_class.csv` | Per-gene four-class label (Class 1/2/3/4) + method-specific significance flags |
| `gene_literature_scores.csv` | PubMed hit counts (ER-specific, BC-general), Tier assignment, alias list |
| `final_ranked_candidates.csv` | Class 2 · Tier 3 candidates ranked by \|log2FC\|, post-agentic refinement |

---

## Methods Summary

<details>
<summary><strong>DEG Analysis — click to expand</strong></summary>

| Method | Implementation | Normalization | Test | Threshold |
|---|---|---|---|---|
| DESeq2 | `pydeseq2` | Median-of-ratios | Wald test | padj < 0.05, \|LFC\| ≥ 1 |
| edgeR | `rpy2` + R `edgeR` | TMM | Quasi-likelihood F-test | FDR < 0.05, \|LFC\| ≥ 1 |
| limma-voom | `rpy2` + R `limma` | TMM → voom | eBayes (trend + robust) | ⚠️ Failed (see Limitations) |

Both contrasts run: **ER+ vs Normal** and **ER− vs Normal**.

</details>

<details>
<summary><strong>Four-Class Dysregulation Framework (Li et al.) — click to expand</strong></summary>

| Class | Definition | Biological Meaning |
|---|---|---|
| Class 1 | Same direction in both ER+ and ER− vs Normal | Shared tumor program |
| **Class 2** | **Opposite directions** | **ER-specific divergence — primary candidates** |
| Class 3 | Significant in ER+ only | ER+ specific program |
| Class 4 | Significant in ER− only | ER− specific program |

</details>

<details>
<summary><strong>Literature Confidence Tiers — click to expand</strong></summary>

| Tier | Criteria | Interpretation |
|---|---|---|
| Tier 1 | ≥ 3 ER-specific PubMed hits | Well-validated in ER context |
| Tier 2 | ≥ 3 breast cancer hits, < 3 ER-specific | Breast cancer literature, not ER-specific |
| **Tier 3** | **< 3 hits in either category** | **Reproducible but understudied — novel candidates** |

Alias resolution via **MyGene.info** (symbol, alias, other_names, name).
Agentic refinement re-queries with expanded aliases for all Class 2 · Tier 3 genes.

</details>

---

## Known Limitations

- **limma-voom failure:** The `calcNormFactors` function was renamed to `normLibSizes` in edgeR v4, breaking the rpy2 bridge. limma-voom returned 0 DEGs. All consensus analyses use the **DESeq2 ∩ edgeR fallback set**. Fix: implement custom TMM normalization or downgrade edgeR.
- **Single dataset:** Results are specific to TCGA-BRCA. Generalizability to other cohorts is not established.
- **Adjacent-normal tissue:** TCGA Solid Tissue Normal samples are adjacent-to-tumor, which may carry residual tumor-field effects.
- **Robustness:** Subsample stability testing was performed for DESeq2 only.

---

## Team

| Member | Sections |
|---|---|
| Y Nhi Tran | Introduction · Related Work · Conclusion · References |
| Yasmeen Soe | Abstract · Methodology · Limitations & Future Work |
| Katrina Huynh | Results & Experiments · Discussion |

**Group 8** — Final Project, Computational Biology / Bioinformatics

---

## References

Key citations (full list in the paper):

- Love et al. (2014). **DESeq2.** *Genome Biology*, 15(12):550.
- Robinson et al. (2010). **edgeR.** *Bioinformatics*, 26(1):139–140.
- Law et al. (2014). **limma-voom.** *Genome Biology*, 15(2):R29.
- Li et al. (2013). Four-class dysregulation framework. *(verify citation)*
- Perou et al. (2000). Molecular portraits of breast tumours. *Nature*, 406:747–752.
- Soneson & Delorenzi (2013). DEG method comparison. *BMC Bioinformatics*, 14:91.
- Schurch et al. (2016). Biological replicates in RNA-seq. *RNA*, 22(6):839–851.

---

## License

This project is for academic purposes only. Data obtained from TCGA via the GDC API is subject to the [TCGA Data Use Agreement](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga/using-tcga/using-tcga-data).
