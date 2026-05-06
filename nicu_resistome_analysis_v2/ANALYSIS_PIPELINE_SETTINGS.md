# NICU Resistome Analysis v2 — Pipeline Settings

**Date:** 2026-05-05
**Companion to:** `RADAR_PIPELINE_PARAMETERS.md` (which documents the RADAR detection step that produces `radar_runs/unified/<sid>/`).
**Scope:** This document covers everything *downstream* of the RADAR per-sample outputs — the aggregator, the strict/permissive filtering, and the analysis pipeline (`scripts/python/0[1-9]_*.py`, `1[0-8]_*.py`).

---

## 1. Aggregator (`scripts/aggregate_radar_outputs.py`)

Runs against `radar_runs/unified/`, writes `data/`. Two behaviors worth pinning:

- **`gene_name` is aliased to `gene_family`.** Downstream scripts group by `gene_name` for differential abundance, longitudinal, and mixed-model analyses, so this places the analysis at the curated **family level** (e.g. all `blaCTX-M-*` alleles collapse to a single `blaCTX-M` row). Allele-level identifiers are preserved as `gene_symbol`. This avoids the ~50× multiple-testing inflation of allele-level tests and matches v1's analysis granularity.
- **`gene_family` and other annotation columns are NaN in the per-sample RADAR TSV** (regression from the simplified DIAMOND DB headers — see `RADAR_PIPELINE_PARAMETERS.md` §"Important schema note"). The aggregator restores them via metadata join. The relevant patches are documented inline; both `gene_family` and `gene_name` (alias) are 100% populated in `data/nicu_amr_only.tsv` after aggregation.

---

## 2. v1-comparable filter (recommended default for manuscript analyses)

The DIAMOND-based v2 RADAR detects ~4× more AMR rows per sample than v1's k-mer + extension `amr_detection`. To compare against v1 numbers (and to exclude the `FRAGMENT` confidence tier that has no v1 analog), apply this filter to `data/nicu_amr_only.tsv` before running `scripts/python/01_*.py` onward:

```python
strict = df[
    (df['detection_confidence'] == 'HIGH') &
    (df['percent_coverage'] >= 80) &
    (df['mean_identity']    >= 90)
]
```

Notes on the thresholds:
- **`detection_confidence == 'HIGH'`** is the primary gate. RADAR's HIGH ≈ v1's HIGH+MODERATE in practice — keeps ~38% of rows.
- **`percent_coverage >= 80`** matches v1's `--min-coverage 0.8` flag in `batch_process_nicu_samples.sh`. HIGH-confidence DIAMOND hits already pass this in nearly every case (only ~4k/122k HIGH rows fail), but it's worth keeping for consistency.
- **`mean_identity >= 90`** matches v1's `--min-identity 0.9`. Note that DIAMOND identity is on the **0–100 scale** (not 0–1 like the v1 flag), and it's protein-level rather than nucleotide.

After filtering, recompute `data/sample_summary.tsv` from the filtered long-format table by summing `total_reads` and `rpm` within each `(sample_name, gene_family)` group, then computing per-sample richness, Shannon, drug-class count.

The full apply-filter-and-recompute snippet is preserved in the project log under `logs/run_amr_*.summary.log`.

### Permissive mode (kept for reference)

The unfiltered v2 outputs (HIGH + MEDIUM + FRAGMENT, all detected hits) are preserved at:

- `data/nicu_amr_only_permissive.tsv`
- `data/sample_summary_permissive.tsv`

To switch back to permissive:

```bash
cp data/nicu_amr_only_permissive.tsv data/nicu_amr_only.tsv
cp data/sample_summary_permissive.tsv data/sample_summary.tsv
# then re-run scripts/python/0[1-9]_*.py
```

The strict view should be the manuscript default unless we want to claim greater detection sensitivity than v1.

---

## 3. Effect of the filter on top-line stats (2026-05-05 run, n=604 samples)

| Metric | v1 (k-mer) | v2 permissive | v2 strict |
|---|---|---|---|
| Gene families pass QC | n/a (alleles) | 473 | **250** |
| Median Axilla `total_amr_rpm` | 609 | 3,110 | 2,692 |
| Median Stool `total_amr_rpm` | 3,357 | 7,849 | 6,559 |
| Stool/Axilla median ratio | 5.5× | 2.5× | 2.4× |
| Mixed-model Stool coef | +0.615, p<10⁻²⁰ | +0.316, p=0.018 | +0.299, p=0.032 |
| Mixed-model Groin coef | +0.008, ns | -0.470, p=0.0005 | -0.530, p=0.0002 |
| **Groin × Wk3 interaction** | not in v1 | not significant | **+0.523, p=0.0075** |
| PCA PC1 variance | (hardcoded "30.8%" in v1 report) | 24.6% | **30.6%** |
| ICC | (n/a) | 0.092 | 0.110 |
| Spearman ρ abx-AMR (overall) | (n/a) | -0.020 | -0.004 |

The strict filter recovers PCA PC1 variance to ~v1 levels and surfaces a real **groin × Week 3 interaction** that the permissive run masked. The body-site magnitude difference (stool/axilla 2.4× vs v1's 5.5×) does **not** recover under filtering — that's a genuine alignment-engine difference, not a filtering artifact.

---

## 4. ARO and MEGARes annotations — present but unused

The aggregator joins `data/amr_protein_db/amr_hierarchy_annotated.tsv` and writes the following columns into `nicu_amr_only.tsv`. As of 2026-05-05 these are populated in the data but **not consumed by any of the analysis scripts** (`scripts/python/0*.py`, `1*.py`).

| Column | % populated | Example values |
|---|---|---|
| `aro_id` | 49.9% | `ARO:3003953`, `ARO:3000822` |
| `aro_match_method` | 49.9% | `gene_name`, `gene_family`, `gene_prefix` |
| `subclass` | 85.9% | `EFFLUX`, `QUINOLONE`, `MUPIROCIN` |
| `scope` | 85.9% | `core`, `plus` |
| `resistance_type` | 85.9% | `intrinsic`, `acquired` |
| `megares_type` | 85.6% | `Multi-compound`, `Drugs`, `Biocides` |
| `megares_class` | 85.6% | `Drug and biocide resistance`, `Mupirocin` |
| `megares_class_source` | 85.6% | `megares_lookup`, `ncbi_fallback` |
| `megares_mechanism` | 35.7% | `Drug and biocide MATE efflux pumps` |
| `megares_group` | 35.7% | `HMRM`, `PMRA`, `MUPA` |
| `megares_meg_id` | 35.7% | `MEG_3268`, `MEG_5797` |

ARO coverage is partial because CARD does not annotate every NCBI ReferenceGeneCatalog entry (heavy-metal, virulence, and some intrinsic-resistance accessions are skipped). MEGARes coverage at the type/class level is broad but the mechanism/group/meg_id rows reflect MEGARes proper, with the ~50% gap filled by NCBI fallback at the class-source level only.

**Status (2026-05-05): all three TODOs done.** The analysis pipeline now consumes ARO and MEGARes annotations via three new scripts. They use the strict-filtered `nicu_amr_only.tsv` and the `master_metadata_with_qc.tsv` produced by 02:

- **MEGARes-class differential abundance**: `scripts/python/04b_megares_class_differential.py`. Per-class linear mixed-effects models with fixed effects for BodySite, SampleCollectionWeek, PostNatalAbxCohort, MaternalAntibiotics, AnyMilk, GestationTime; SubjectID random intercept. FDR-corrected per covariate across classes. Outputs at `results/differential_megares/megares_class_clinical_effects_{long,wide,topN}.tsv` and a heatmap at `figures/differential_megares/megares_class_clinical_effects_heatmap.pdf`. (No Location-vs-Location framing — that contrast is not the headline for this manuscript.)
- **Intrinsic vs acquired split**: `scripts/python/09b_intrinsic_acquired_split.py`. Splits per-sample RPM totals by NCBI's `resistance_type`; emits `results/summary/intrinsic_vs_acquired_{per_sample,summary}.tsv` and `acquired_fraction_per_sample.tsv`, plus a 3-panel figure. The summary report (09) now reads these and adds a "Section 6: Intrinsic vs Acquired Resistance" block.
- **ARO-level panel**: `scripts/python/04c_aro_panel.py`. Per-ARO median RPM and prevalence stratified by SampleType × Week and by PostNatalAbxCohort, with `gene_family`, `subclass`, `resistance_type`, and `megares_class` retained for cross-reference to AMRFinderPlus / CARD. Outputs at `results/aro_panel/` and a top-30 heatmap at `figures/aro_panel/`.

The Location-vs-Location framing in 04 is preserved but should not be the manuscript's headline; the clinical contrasts in 04b are the right ones for this question (body site, postnatal age, infant antibiotic cohort, maternal antibiotics, feeding, gestational age).

---

## 4b. ZCH sample-name canonicalization (2026-05-05)

The pipeline previously had two collisions in ZCH naming that silently dropped all 240 ZJH AMR samples from analysis:

1. The metadata `sample_name` column was set from the unprefixed `Sample` column (`N01_1_2`) instead of the prefixed `SampleID` (`ZJH_N01_1_2`), so the metadata→AMR join couldn't match ZCH rows.
2. ZCH single-digit subjects were padded to two-digit (`N1` → `N01`) at the manifest stage, even though the raw FASTQs and clinical convention use the single-digit form.

Fixed by `scripts/canonicalize_zjh_names.py` (rewrites `data/nicu_amr_only*.tsv`, `data/sample_summary*.tsv`, `data/sample_manifest.csv`, and the source `metadata/AllNICUSampleKeyRevised20250320*.csv`) and `scripts/canonicalize_radar_dirs.sh` (renames the 48 `radar_runs/unified/ZJH_N0X_*/` dirs and their inner files). Originals preserved as `*.pre_canonicalize`. Companion fix in `01_create_master_metadata.py`: removed the leading-zero pad on the antibiotics `Subject` column.

`02_quality_control.py` was also fixed to group on `SubjectID` (already disambiguated as `N01..N68` for UCMC vs `ZJH_N1..ZJH_N59` for ZCH) instead of `PatientID` (which collides between the two NICUs for `N10..N59`).

Finally, `data/sample_summary.tsv` was augmented with zero-rows for the 54 radar-processed samples (43 ZCH + 11 UCMC) that had no strict-filter-passing AMR detections — those samples were sequenced and should remain in the cohort for body-site / postnatal-age / antibiotic comparisons rather than be dropped by selection bias.

Net effect: post-canonicalization the cohort is 642 samples (359 UCMC + 283 ZCH) instead of the prior buggy 347 (UCMC-only).

---

## 5. Hard-coded numbers — eliminated as of 2026-05-05

The original `09_create_summary_report.py` had hard-coded headline numbers (gene counts, mixed-model coefficients, ICC, PC1 variance, Spearman ρ) that did not update across runs. These are now read from source files:

| Field | Source |
|---|---|
| Gene families analyzed | `results/qc/genes_passing_qc.txt`, `results/qc/gene_prevalence_stats.tsv` |
| Mixed-model Location p, Stool/Groin coefs | `results/mixed_models/full_model_coefficients.tsv` |
| PCA PC1 variance | `results/exploratory/pca_variance_explained.tsv` |
| Overall Spearman ρ (abx vs AMR) | `results/correlations/total_abx_amr_correlation_overall.tsv` (newly emitted by `05_antibiotic_correlations.py`) |
| ICC | `results/mixed_models/variance_components.tsv` (newly emitted by `07_mixed_effects_models.py`) |
| Complete subjects per location | computed from `summary_by_location` |

The narrative section titles ("Significant Longitudinal Changes in Groin and Stool", "Body Site is the Primary Driver of Resistome Composition") and the Conclusions section remain as author-editable prose — these are interpretation, not numbers, and shouldn't be auto-generated.
