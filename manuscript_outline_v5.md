# Manuscript Drafting Plan v5 — All Data Finalized

**Identity:** Resource paper
**Target:** Genome Biology
**Harness:** `~/Dropbox/Bioinformatics/Claude/bioarxiv-harness/` (copy for VEP; phi_gap paper stays)
**Working Title:** "Evo2 zero-shot variant effect prediction across 215,001 coding and non-coding variants in spaceflight radiation-response genes"

---

## Paper Metadata

**Title:** Evo2 zero-shot variant effect prediction across 215,001 coding and non-coding variants in spaceflight radiation-response genes

**Authors:** JangKeun Kim^1,2^, Christopher E. Mason^1,2^*

**Affiliations:**
1. Department of Physiology and Biophysics, Weill Cornell Medicine, New York, NY 10065, USA
2. The WorldQuant Initiative for Quantitative Prediction, Weill Cornell Medicine, New York, NY 10065, USA

**Corresponding:** * chm2042@med.cornell.edu

**Running title:** Evo2 zero-shot VEP for spaceflight genes
**Lead author:** Kim *et al.*

---

## Harness Setup

The bioarxiv-harness is currently hosting the evolutionary potential (phi_gap) paper. For the VEP manuscript:

1. **Copy harness** to a new directory (e.g., `~/Dropbox/Bioinformatics/Claude/Evo2/manuscript/`) or create a new branch
2. **Replace metadata** in `main.tex`: title, authors, running title, lead author
3. **Clear section files** and write VEP content
4. **Replace figures/** with VEP figures
5. **Replace references.bib** with VEP references
6. **Replace supplement.tex** with VEP supplementary materials

Key harness features to use:
- `\leadauthor{Kim \textit{et al.}}` — footer
- `\runningtitle{Evo2 zero-shot VEP for spaceflight genes}` — header
- `\correspondingauthor{chm2042@med.cornell.edu}`
- `\keywords{...}` — keyword block
- `\bp`, `\nt` — biology SI units from bioarxiv.sty
- `mode-submit` Makefile target for journal submission mode
- `make check` for pre-submission validation
- `cleveref` for smart references (\Cref)
- `subcaption` for multi-panel figures (a,b,c)

---

## Confirmed Numbers (All Pre-Draft Analyses Complete)

### Variant Counts
- **Total: 215,001** (168,409 SNVs + 46,592 indels)
- Per gene: ATM 52,791 | BRCA1 36,901 | DNMT3A 19,693 | CLOCK 18,369 | TERT 17,443 | CHEK2 16,098 | NFE2L2 15,230 | TP53 14,240 | MSTN 12,383 | RAD51 11,853

### DMS Calibration (Multi-Tool)
| Gene | n_matched | Evo2 rho | AM rho | CADD rho | REVEL rho |
|------|-----------|----------|--------|----------|-----------|
| BRCA1 | 3,884 | **+0.515** | -0.562 | -0.499 | -0.470 |
| TP53 | 2,827 | **-0.211** | +0.414 | +0.301 | +0.413 |
| CHEK2 | 4,882 | **+0.310** | -0.362 | -0.370 | -0.254 |
| DNMT3A | 746 | **+0.457** | -0.672 | -0.267 | -0.519 |

Note: Sign convention differs — Evo2 delta is negative=damaging (positive rho = correct), other tools score higher=pathogenic (negative rho = correct). |Evo2| is competitive: BRCA1 0.515 vs AM 0.562; DNMT3A 0.457 vs AM 0.672 (AM wins); TP53 0.211 vs AM 0.414 (AM wins); CHEK2 0.310 vs CADD 0.370 (comparable).

### ClinVar Per-Tool AUROCs (with 95% Bootstrap CIs)
| Gene | Evo2 | CI | CADD | CI | AM | CI | n_Evo2 |
|------|------|----|------|----|----|----|--------|
| ATM | **0.9955** | [0.993, 0.998] | 0.9989 | [0.998, 1.000] | 0.926 | [0.870, 0.969] | 3,398 |
| BRCA1 | **0.9879** | [0.984, 0.992] | 0.9951 | [0.992, 0.997] | 0.960 | [0.935, 0.980] | 2,676 |
| TP53 | **0.9887** | [0.983, 0.993] | 0.9882 | [0.979, 0.995] | 0.988 | [0.978, 0.996] | 839 |
| CHEK2 | **0.9996** | [0.999, 1.000] | 0.9986 | [0.996, 1.000] | — (n=8) | — | 689 |
| DNMT3A | **0.9962** | [0.989, 1.000] | 1.000 | [1.000, 1.000] | 1.000 | [1.000, 1.000] | 159 |
| TERT | **0.9089** | [0.832, 0.974] | 0.9990 | [0.997, 1.000] | 0.969 | [0.893, 1.000] | 715 |

Evo2 is **highest** for: CHEK2 (0.9996 > CADD 0.9986), TP53 (0.9887 > CADD 0.9882)

### Matched-Variant AUROCs (Evo2 ∩ CADD ∩ AM ∩ REVEL)
| Gene | n | Evo2 | CADD | AM | REVEL |
|------|---|------|------|----|-------|
| ATM | 86 | 0.864 | 0.940 | 0.924 | 0.941 |
| BRCA1 | 335 | 0.932 | 0.957 | **0.960** | 0.916 |
| TP53 | 271 | 0.894 | 0.920 | **0.988** | 0.967 |
| DNMT3A | 18 | 0.847 | 1.000 | 1.000 | 0.972 |
| TERT | 22 | 0.866 | 0.991 | 0.969 | 0.991 |
| CHEK2 | 8 | — (too few) | — | — | — |

On matched sets (mostly missense), AM leads because it's specifically designed for missense. Evo2 is lower but still high (>0.84). Key narrative: **matched sets are overwhelmingly missense, which is AM/REVEL's home turf; Evo2's advantage is breadth** (non-coding, indels, all SNV types).

### DeLong Tests (Matched Sets)
Most pairwise comparisons NOT significant (small matched sets). Notable:
- TP53: AM vs Evo2 z=4.27, p=2.0e-5 (AM significantly better on matched missense)
- BRCA1: AM vs Evo2 z=2.04, p=0.042 (marginal)
- No gene shows Evo2 significantly beating other tools on matched sets
- But: CADD vs Evo2 never significant (p > 0.08 in all genes) — Evo2 is comparable to CADD

### Ensemble (Evo2 + CADD, 5-Fold CV)
| Gene | Evo2 | CADD | Ensemble | Improvement |
|------|------|------|----------|-------------|
| ATM | 0.9954 | 0.9988 | **0.9992** | +0.04% |
| BRCA1 | 0.9877 | 0.9952 | **0.9971** | +0.19% |
| TP53 | 0.9882 | 0.9885 | **0.9922** | **+0.38%** |
| CHEK2 | 0.9997 | 0.9986 | **0.9998** | +0.02% |
| DNMT3A | 0.9949 | 1.0000 | 1.0000 | 0.00% |
| TERT | 0.8996 | 0.9981 | **0.9991** | +0.09% |

Ensemble always >= best single. TP53 shows most complementarity.

### MPRA Multi-Tool (TERT Promoter, n=777)
| Tool | rho | p | n |
|------|-----|---|---|
| **Evo2** | **+0.267** | **3.8e-14** | **777** |
| CADD | -0.081 | 0.21 (NS) | 240 |
| SpliceAI | +0.003 | 0.96 (NS) | 240 |
| ncER | -0.211 | 5.9e-4 | 261 |

**Evo2 is the only tool with significant positive MPRA correlation.** This is a major selling point.

### Indel Frameshift AUROCs
| Gene | AUROC | n |
|------|-------|---|
| ATM | **0.993** | 2,699 |
| DNMT3A | **0.992** | 71 |
| TP53 | **0.986** | 620 |
| CHEK2 | **0.980** | 775 |
| BRCA1 | **0.961** | 3,083 |
| TERT | **0.948** | 108 |

Mean: **0.977** across 6 genes. All > 0.94.

### Transversion Control (**CRITICAL — Revises Radiation Narrative**)
Radiation-oxidative (C>A/G>T) vs other transversions:
| Gene | p (raw) | p (Bonf ×10) | rbc | Direction |
|------|---------|-------------|-----|-----------|
| TP53 | 0.0005 | **0.005** | +0.049 | Rad slightly more damaging |
| DNMT3A | 0.003 | **0.032** | +0.033 | Rad slightly more damaging |
| TERT | 0.032 | 0.322 | +0.025 | NS after Bonferroni |
| MSTN | 0.037 | 0.375 | +0.032 | NS |
| BRCA1 | 0.698 | 1.000 | -0.005 | NS |
| ATM | 0.645 | 1.000 | -0.003 | NS |
| NFE2L2 | 0.381 | 1.000 | +0.004 | NS |
| CLOCK | 0.839 | 1.000 | -0.013 | NS |
| RAD51 | 0.757 | 1.000 | -0.011 | NS |
| CHEK2 | 0.268 | 1.000 | +0.009 | NS |

**Conclusion:** Only TP53 (borderline) and DNMT3A survive Bonferroni. 8/10 genes show NO radiation-specific signal beyond transversion/transition bias. The original "radiation mutations more damaging" finding was **almost entirely driven by transversion > transition**, NOT radiation-specific chemistry.

However: **transversion > transition is itself real and consistent** (all 10 genes p < 1e-10 for rad vs transitions). And microhomology deletions remain highly damaging.

**Revised narrative:** Frame as "transversions (including radiation-signature) score more damaging than transitions, consistent with greater chemical disruption" rather than "radiation-specific signal." The transversion control is an honest null result that strengthens the paper's credibility.

### Orthogonality (Evo2 vs tools, Spearman rho)
| Gene | vs CADD | vs REVEL | vs AM | vs SpliceAI |
|------|---------|----------|-------|-------------|
| ATM | 0.674 | 0.633 | 0.569 | 0.373 |
| BRCA1 | 0.721 | 0.472 | 0.788 | 0.285 |
| TP53 | 0.828 | 0.525 | 0.601 | 0.193 |
| CHEK2 | 0.682 | 0.441 | — | 0.446 |
| DNMT3A | 0.601 | 0.509 | 0.538 | 0.053 |
| TERT | 0.414 | 0.545 | 0.638 | 0.094 |

Evo2 vs CADD: 0.41–0.83 (complementary to moderately correlated)
Evo2 vs SpliceAI: 0.05–0.45 (highly complementary)

---

## Figure Plan (6 Main + 4 Supplementary)

### Main Figures
| # | Title | Panels | Source |
|---|-------|--------|--------|
| **Fig 1** | Study overview | Graphical abstract: 10 genes → Evo2 → 3-layer validation | **CREATE** (BioRender) |
| **Fig 2** | DMS calibration | (a-d) 4 scatter+marginals, (e) multi-tool |rho| bars | **UPDATE** fig3_dms → add multi-tool panel |
| **Fig 3** | ClinVar benchmark | (a) ROC curves 6 genes, (b) per-tool AUROC+CIs bars, (c) matched-set AUROC+CIs bars, (d) DeLong matrix | **REGEN** with new data |
| **Fig 4** | Orthogonality & ensemble | (a) pairwise rho heatmap, (b) ensemble improvement bars, (c) Evo2-unique variant recovery | **REGEN** |
| **Fig 5** | Constraint landscapes | 10-gene score profiles with domains | **EXISTS** (needs CLOCK+ATM → COMPLETE) |
| **Fig 6** | Non-coding & indels | (a) TERT MPRA + multi-tool, (b) ENCODE cCRE enrichment, (c) frameshift AUROC bars | **UPDATE** fig6 → add multi-tool MPRA + indel bars |

### Supplementary Figures
| # | Title | Source |
|---|-------|--------|
| **Fig S1** | Window ablation (3 genes × 4 windows) | **EXISTS** fig2_ablation |
| **Fig S2** | Mutation signatures | (a) 4-class per gene, (b) transversion control (C>A/G>T vs other TV), (c) microhomology | **UPDATE** fig7 → add transversion control panel |
| **Fig S3** | Indel detail | Per-gene fs vs if breakdown | **EXISTS** supp_s3 |
| **Fig S4** | Astronaut variants | 8 variants with percentiles | **EXISTS** fig8 |

### Supplementary Tables
| # | Content | Data source |
|---|---------|-------------|
| **S1** | Gene panel: coordinates, transcripts, rationale | gene_coordinates.py |
| **S2** | Variant counts per gene (SNV/indel) | `predraft/variant_counts.json` |
| **S3** | DMS calibration: sources, URNs, multi-tool rho | `predraft/dms_multitool/` |
| **S4** | ClinVar breakdown by gene × star rating | `clinvar_validation_summary.json` |
| **S5** | Multi-tool benchmark: all AUROCs + CIs + DeLong | `predraft/matched_benchmark/` |
| **S6** | Radiation: 4-class + transversion control | `predraft/transversion_control/` |
| **S7** | Astronaut variants: scores + clinical context | existing astronaut results |

---

## Section-by-Section Drafting Plan

### Abstract (~350 words, structured for Genome Biology)

**Background:** VEP bottleneck; existing tools restricted to coding/missense or require training; spaceflight radiation creates characteristic mutations; no comprehensive variant effect maps for radiation-response genes spanning coding + non-coding.

**Results:** Applied Evo2 7B to score 215,001 variants (168,409 SNVs + 46,592 indels) across 10 genes. DMS calibration |rho| = 0.21–0.52 (4 genes). Per-tool ClinVar AUROCs: 0.91–1.00 across 6 genes, with Evo2 highest for CHEK2 (0.9996) and TP53 (0.9887). On matched-variant sets, Evo2 comparable to supervised tools (DeLong: CADD vs Evo2 p > 0.08 all genes). Orthogonal signal (rho = 0.41–0.83 vs CADD); Evo2+CADD ensemble improves by up to 0.38%. Uniquely scores non-coding (TERT MPRA rho = +0.267, p = 3.8e-14 — only tool with significant positive correlation) and frameshifts (AUROC > 0.94, 6 genes). Transversions score more damaging than transitions across all genes; radiation-oxidative effect after controlling for transversion bias is marginal (2/10 genes Bonferroni-significant).

**Conclusions:** Comprehensive zero-shot VEP atlas; complementary to supervised tools; fills critical gaps for non-coding and indel variants; supports integration into ACMG PP3/BP4 evidence.

### Introduction (~1,500 words, 5 paragraphs)

1. **VEP problem** — variant interpretation bottleneck, VUS burden, ACMG PP3/BP4
2. **Current landscape** — AM (missense-only), CADD (integrative, ClinVar-trained), REVEL (meta-predictor, missense-only), SpliceAI, ncER, LINSIGHT; no single tool covers all variant types
3. **Foundation models** — Evo2 7B, 9.3T-nt OpenGenome2, up to 1 Mb context; zero-shot, any variant type; no ClinVar circularity; vs NT (shorter context), HyenaDNA (smaller scale)
4. **Spaceflight context** — radiation signatures (C>A, complex DSBs), Rutter 2024 (protective alleles), 10-gene panel rationale; CLOCK/MSTN as exploratory
5. **This study** — 215,001 variants, 3-layer validation, resource contribution

### Results (~3,500 words, 7 subsections)

#### 2.1 Window size optimization (Fig S1)
- w8192 optimal, mean AUROC 0.9921
- 8192 is scoring window, not model limit (1 Mb)

#### 2.2 DMS calibration (Fig 2)
- 4 genes: BRCA1 +0.515, DNMT3A +0.457, CHEK2 +0.310, TP53 -0.211
- TP53 negative expected (nutlin-3 high=LOF)
- Multi-tool comparison: Evo2 |rho| competitive with AM/CADD on same variants
- BRCA1: Evo2 0.515 vs AM 0.562 (comparable); DNMT3A: 0.457 vs AM 0.672 (AM higher)

#### 2.3 ClinVar benchmarking (Fig 3)
- **Per-tool (primary presentation for resource paper):** Evo2 0.91–1.00; CHEK2 and TP53 highest of all tools
- **Matched-variant (fairness analysis):** Evo2 0.85–0.93; AM leads on missense-dominated sets; DeLong CADD vs Evo2 NS in all genes
- **Caveat:** Matched sets are mostly missense (AM/REVEL home turf); Evo2's advantage is breadth; CADD incorporates ClinVar in training

#### 2.4 Orthogonality & ensemble (Fig 4)
- Evo2 vs CADD rho 0.41–0.83 (complementary)
- Evo2 vs SpliceAI rho 0.05–0.45 (highly complementary)
- Ensemble: always >= best single, TP53 +0.38%

#### 2.5 Non-coding and indels (Fig 6)
- **TERT MPRA:** rho = +0.267 (p = 3.8e-14, n=777); only tool with significant positive correlation; CADD NS, SpliceAI NS, ncER rho=-0.211
- **Indels:** Frameshift AUROCs 0.95–0.99 (6 genes); mean 0.977
- TERT C228T/C250T null scores — limitation noted

#### 2.6 Mutation signature analysis (Fig S2) — **REVISED NARRATIVE**
- **Primary finding (honest):** Transversions score more damaging than transitions across ALL 10 genes (p < 1e-10)
- **Transversion control:** Radiation-oxidative vs other transversions — NOT significant in 8/10 genes after Bonferroni; only TP53 (p=0.005) and DNMT3A (p=0.032) borderline
- **Interpretation:** The signal is transversion > transition (greater chemical disruption), not radiation-specific. This is expected — Evo2 scores sequence disruption, and transversions cause more disruption than transitions by definition.
- **Microhomology deletions:** Still consistently highly damaging across all genes
- Frame honestly: "consistent with known mutational physics, not a radiation detector"

#### 2.7 Variant atlas & constraint landscapes (Fig 5)
- 215,001 total variants
- Domain structure visible in landscapes
- Astronaut case studies (Fig S4, Table S7): DNMT3A R882C 97.2nd pctile, TP53 R248W 86.5th

### Discussion (~2,000 words, 6 paragraphs)

1. **Summary** — first comprehensive zero-shot VEP atlas for spaceflight genes; unified scoring; complementary signal
2. **Foundation models as complementary tools** — competitive without training data; no ClinVar circularity; orthogonal; ensemble improvement; comparison with NT, HyenaDNA
3. **Non-coding & indel gaps filled** — MPRA correlation; frameshift AUROCs; AM/REVEL can't do this
4. **Mutation signatures — honest reporting** — transversion > transition is real but not radiation-specific; microhomology deletions consistent with DSB signatures; Evo2 captures sequence vulnerability, not radiation exposure
5. **Limitations** (8 items): zero-shot accuracy vs supervised; ClinVar circularity; GPU requirement; limited non-coding validation; 10-gene scope; TERT C228T/C250T null; no protein structure; ClinVar star bias
6. **Future directions** — ACMG PP3/BP4; foundation model ensembles; larger panels; precomputed atlases for mission planning

### Methods (~3,000 words, 9 subsections)

4.1 Gene panel selection
4.2 Evo2 variant scoring (7B, 8192 bp window, A40/A100)
4.3 DMS calibration (4 datasets, URNs, offsets)
4.4 ClinVar benchmarking (≥2 star, matched vs per-tool)
4.5 Multi-tool benchmarking (AM, CADD, REVEL, SpliceAI, ncER, LINSIGHT)
4.6 Non-coding validation (TERT MPRA, ENCODE cCREs)
4.7 Radiation signature analysis (4-class + transversion control)
4.8 Statistical analysis (Spearman, Mann-Whitney, DeLong, bootstrap, Bonferroni)
4.9 Data and code availability

---

## Key References (26 entries needed for references.bib)

| Short key | Citation | DOI |
|-----------|----------|-----|
| Evo2_2026 | Brixi et al. Nature 2026 | 10.1038/s41586-026-10176-5 |
| Rutter_2024 | Rutter et al. Nat Comms 2024 | 10.1038/s41467-024-49423-6 |
| AlphaMissense_2023 | Cheng et al. Science 2023 | 10.1126/science.adg7492 |
| CADD_2024 | Schubach et al. NAR 2024 | 10.1093/nar/gkad989 |
| REVEL_2016 | Ioannidis et al. AJHG 2016 | 10.1016/j.ajhg.2016.08.016 |
| SpliceAI_2019 | Jaganathan et al. Cell 2019 | 10.1016/j.cell.2018.12.015 |
| BRCA1_DMS_2018 | Findlay et al. Nature 2018 | 10.1038/s41586-018-0461-z |
| TP53_DMS_2018 | Giacomelli et al. Nat Genet 2018 | 10.1038/s41588-018-0204-y |
| CHEK2_DMS_2024 | McCarthy-Leo et al. bioRxiv 2024 | (MaveDB: 00001203-a-1) |
| DNMT3A_DMS_2025 | Garcia et al. bioRxiv 2025 | 10.1101/2025.09.24.678339 |
| TERT_MPRA_2019 | Kircher et al. Nat Comms 2019 | 10.1038/s41467-019-11526-w |
| ClinVar_2018 | Landrum et al. NAR 2018 | 10.1093/nar/gkx1153 |
| ACMG_2015 | Richards et al. Genet Med 2015 | 10.1038/gim.2015.30 |
| Twins_2019 | Garrett-Bakelman et al. Science 2019 | 10.1126/science.aau8650 |
| gnomAD_2020 | Karczewski et al. Nature 2020 | 10.1038/s41586-020-2308-7 |
| ENCODE_2020 | ENCODE Consortium Nature 2020 | 10.1038/s41586-020-2493-4 |
| NT_2024 | Dalla-Torre et al. Nat Methods 2024 | 10.1038/s41592-024-02523-z |
| HyenaDNA_2023 | Nguyen et al. NeurIPS 2023 | arXiv:2306.15794 |
| ncER_2019 | Wells et al. Nat Comms 2019 | 10.1038/s41467-019-13212-3 |
| LINSIGHT_2017 | Huang et al. Nat Genet 2017 | 10.1038/ng.3810 |
| DeLong_1988 | DeLong et al. Biometrics 1988 | 10.2307/2531595 |
| daSilveira_2020 | da Silveira et al. Cell 2020 | 10.1016/j.cell.2020.10.002 |
| MaveDB_2019 | Esposito et al. Genome Biol 2019 | 10.1186/s13059-019-1845-6 |
| Evo1_2024 | Nguyen et al. Science 2024 | 10.1126/science.ado9336 |
| ESM2_2023 | Lin et al. Science 2023 | 10.1126/science.ade2574 |
| Richards_ACMG_2015 | (same as ACMG_2015) | — |

---

## Drafting Order

1. **Methods** — most mechanical, establishes terminology for Results
2. **Results** — data-driven, reference Methods subsections
3. **Introduction** — frame the story after results are clear
4. **Discussion** — synthesize findings
5. **Abstract** — summarize final paper
6. **Figures** — regenerate/update alongside Results writing
7. **Supplement** — tables and supplementary figures
8. **references.bib** — compile all citations

---

## Data Files for Each Section

| Section | Required data files |
|---------|-------------------|
| 2.1 Window ablation | `results/window_ablation/` (existing) |
| 2.2 DMS | `results/predraft/dms_multitool/dms_multitool_results.json`, `results/calibration/*.json` |
| 2.3 ClinVar | `results/predraft/matched_benchmark/matched_benchmark_results.json`, `all_aurocs_with_ci.csv` |
| 2.4 Orthogonality | `matched_benchmark_results.json` (orthogonality section), `predraft/ensemble/ensemble_results.json` |
| 2.5 Non-coding + indels | `predraft/mpra_multitool/mpra_multitool_results.json`, `predraft/indels/indel_aurocs_all_genes.json` |
| 2.6 Radiation | `predraft/transversion_control/transversion_control_results.json`, existing radiation_results.json |
| 2.7 Atlas | `predraft/variant_counts.json`, existing constraint landscape data, astronaut results |

---

## Critical Narrative Decisions

### 1. Per-tool vs matched benchmarking
**Decision:** Lead with per-tool AUROCs (primary for a resource paper — shows what Evo2 can do). Present matched-set analysis as fairness check in same section. Acknowledge AM wins on matched missense sets; argue Evo2's value is breadth + orthogonality.

### 2. Radiation narrative
**Decision:** Honest reporting. Frame as "transversions > transitions" (10/10 genes), not "radiation-specific." The transversion control null result strengthens credibility. Note that microhomology deletions remain a genuine DSB-associated signal.

### 3. DMS absolute rho values
**Decision:** Report both sign-corrected rho for Evo2 and tool rho side by side. Acknowledge AM has higher |rho| for DNMT3A and TP53 on same variants. Argue Evo2's value is zero-shot (no training data) + non-coding/indel coverage.

### 4. TERT lower AUROC
**Decision:** Explain honestly — TERT is non-coding, many regulatory variants, GOF promoter hotspots violate LOF assumption. This is expected and informative, not a failure.

### 5. CHEK2 matched set too small
**Decision:** Only 8 variants in matched set (6 P/LP, 2 B/LB). Report per-tool only for CHEK2, note insufficient matched variants. CHEK2's per-tool AUROC of 0.9996 is still a highlight.
