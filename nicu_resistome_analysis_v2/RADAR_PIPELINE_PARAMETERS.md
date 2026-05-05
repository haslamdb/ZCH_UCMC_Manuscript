# RADAR Pipeline Parameters — v2 NICU Resistome Re-run

**Run date**: 2026-05-05
**Cohort**: 642 NICU samples (UCMC + ZCH), reusing prior `nicu_resistome_analysis/` cohort
**Output**: `nicu_resistome_analysis_v2/radar_runs/unified/<sample_id>/`
**Pipeline binary**: `biogpu/build/runtime/kernels/unified_pipeline` (orchestrator that runs AMR via DIAMOND + FQ via clean_resistance_pipeline in parallel for each sample)

---

## Hardware

| GPU (PCI) | Device | VRAM | Role |
|-----------|--------|------|------|
| 0 (`01:00.0`) | RTX PRO 6000 Blackwell #1, sm_120 | 96 GB GDDR7 | RADAR worker chunk 0 (214 samples) |
| 1 (`C1:00.0`) | RTX A6000, sm_86 | 48 GB | RADAR worker chunk 1 (214 samples) |
| 2 (`E1:00.0`) | RTX PRO 6000 Blackwell #2, sm_120 | 96 GB GDDR7 | RADAR worker chunk 2 (214 samples) |

`CUDA_DEVICE_ORDER=PCI_BUS_ID` is set so the CLI device index matches `nvidia-smi` (without it, A6000 ends up as CUDA index 2, not 1).

CUDA 12.9, driver 590.48.01.

---

## Input Data

- **FASTQs**: `/bulkpool/sequence_data/mss_data/<sample_id>_R{1,2}.fastq.gz` (paired-end Illumina NextSeq, ~25–30M paired reads/sample for NICU samples)
- **Sample manifest**: `nicu_resistome_analysis_v2/data/sample_manifest.csv`
- **Sample ID rename rule**: ZCH samples in the prior metadata used a `ZJH_N0X_` prefix for subjects 1–9; FASTQs are named `ZJH_NX_` (no leading zero). The manifest builder strips the leading zero only for `ZJH_N0[1-9]_` to resolve all 642 samples.

---

## Databases

| Component | Path | Notes |
|-----------|------|-------|
| AMR protein bloom filter | `biogpu/data/amr_protein_db/protein_bloom.bin` | Pre-built; ~99% read rejection rate before DIAMOND |
| AMR DIAMOND DB | `biogpu/data/amr_protein_db/proteins.dmnd` | NCBI ReferenceGeneCatalog–derived, 9,538 proteins |
| AMR metadata JSON | `biogpu/data/amr_protein_db/protein_details.json` | id-keyed metadata for post-hoc annotation join |
| AMR hierarchy TSV | `biogpu/data/amr_protein_db/amr_hierarchy_by_genename.tsv` | NCBI hierarchy for mechanism / intrinsic flag (used by aggregator only) |
| FQ nucleotide DB | `biogpu/data/fq_resistance_db_v2/nucleotide` | **v2** (validated 2026-03-07: 92.8% sens / 81.1% prec on 49 E. coli isolates) |
| FQ protein DB | `biogpu/data/fq_resistance_db_v2/protein` | Companion to nucleotide DB |
| FQ resistance DB JSON | `biogpu/data/fq_resistance_db_v2/resistance_db.json` | Mutation positions, species mappings. Symlinked from `biogpu/data/resistance_db.json` so unified_pipeline picks it up via the default. |

The unified_pipeline shells out to `clean_resistance_pipeline` from cwd `/home/david/projects/biogpu`, which reads `data/resistance_db.json` (the symlink) by default. There is no `--fq-db-json` flag exposed at the unified_pipeline level.

---

## Per-Sample Invocation (one per chunk, one chunk per GPU)

```bash
cd /home/david/projects/biogpu

CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=<gpu_idx> \
  /home/david/projects/biogpu/build/runtime/kernels/unified_pipeline \
    --csv      <chunk_<gpu_idx>.csv> \
    --output-dir /home/david/projects/ZCH_UCMC_Manuscript/nicu_resistome_analysis_v2/radar_runs/unified \
    --reference-db   /home/david/projects/biogpu/data/fq_resistance_db_v2/nucleotide \
    --resistance-db  /home/david/projects/biogpu/data/fq_resistance_db_v2/protein \
    --protein-bloom  /home/david/projects/biogpu/data/amr_protein_db/protein_bloom.bin \
    --diamond-db     /home/david/projects/biogpu/data/amr_protein_db/proteins \
    --amr-metadata   /home/david/projects/biogpu/data/amr_protein_db/protein_details.json \
    --gpu-device     0
```

`<chunk_N.csv>` columns (round-robin assigned across GPUs by sample order in manifest):

```
SampleName,FilePath,R1 file,R2 file
N01_1_2,/bulkpool/sequence_data/mss_data,N01_1_2_R1.fastq.gz,N01_1_2_R2.fastq.gz
...
```

Each chunk has 214 samples; total 642.

---

## Detection Thresholds

Inherited from biogpu defaults (`unified_pipeline` does not surface tuning flags for these):

- **AMR bloom filter min k-mer hits**: 3 (default)
- **AMR DIAMOND**: default sensitivity, 128 threads
- **AMR coverage_calculator** (post-DIAMOND filter for "detected"):
  - `--min-coverage` 80 (% protein covered)
  - `--min-depth` 5 (mean depth)
  - `--all-features` enabled (full ML feature set in `_gene_stats.tsv`)
- **FQ resistance**: bloom enabled (50% pass rate observed on N01_1_2 smoke), Smith-Waterman enabled, `--min-allele-depth 5`, default mutation thresholds.

---

## Per-Sample Outputs

Per-sample directory `radar_runs/unified/<sid>/`:

| File | Contents |
|------|----------|
| `<sid>_gene_stats.tsv` | Per-gene AMR statistics (440 rows total, ~295 detected for a typical NICU sample) |
| `<sid>_amr_alignments.tsv` | Raw DIAMOND BLAST-X hits |
| `<sid>_detected_genes.txt` | Filtered AMR gene list (those passing coverage/depth thresholds) |
| `<sid>_sample_stats.txt` | Sample-level AMR summary (total reads, gene count, drug class abundances, Shannon diversity) |
| `<sid>_allele_frequencies.tsv` | FQ resistance per-position allele frequencies (~80 rows for a typical NICU sample) |
| `<sid>_protein_matches.tsv` | Raw FQ protein-level alignments |
| `<sid>_clinical_fq_report.{json,html,txt}` | Clinical interpretation of FQ findings |
| `<sid>.h5` | HDF5 with raw alignment tables (kept; not deleted) |
| `<sid>_amr_filtered_R{1,2}.fastq.gz` | Bloom-filtered reads passed to DIAMOND |

### Important schema note: `gene_stats.tsv` annotation columns are empty by design

The `gene_symbol`, `gene_family`, `drug_class`, `mechanism`, `gene_type` columns in `_gene_stats.tsv` come back empty from the `coverage_calculator` subprocess. This is because `coverage_calculator` parses DIAMOND alignment headers expecting an 11-part `gene_id|accession|x|y|gene_symbol|gene_family|mechanism|z|drug_class|type|description` format, while the actual DIAMOND DB headers use a compact 4-part `id|gene_symbol|drug_class|accession` form (e.g. `1011|merA|MERCURY|CAC69251.1`).

The information is fully recoverable post-hoc; the aggregator handles this:

1. Parse the pipe-delimited `protein_id` column → recover `id`, `gene_symbol`, `drug_class`, `accession`.
2. Join with `protein_details.json` (id-keyed, 9,538 entries) → recover `gene_family`.
3. Join with `amr_hierarchy_by_genename.tsv` → recover `resistance_mechanism`, `is_intrinsic`, NCBI subclass.

This is implemented in `nicu_resistome_analysis_v2/scripts/aggregate_radar_outputs.py`.

---

## Aggregation

Runs after the batch completes:

```bash
python /home/david/projects/ZCH_UCMC_Manuscript/nicu_resistome_analysis_v2/scripts/aggregate_radar_outputs.py
```

Produces (in `nicu_resistome_analysis_v2/data/`):

- `nicu_amr_only.tsv` — long-format `sample_name × gene` AMR matrix (matches prior schema with hierarchy columns added)
- `nicu_fq_allele_frequencies.csv` — concatenated FQ allele frequency table
- `nicu_fq_protein_matches.csv` — concatenated FQ protein matches
- `sample_summary.tsv` — per-sample summary stats

---

## Smoke Results (sample `N01_1_2`)

| Metric | Prior analysis | RADAR v1 (integrated_clean_db) | RADAR v2 (fq_resistance_db_v2) |
|--------|----------------|--------------------------------|--------------------------------|
| AMR genes detected | 96 | 295 | 295 |
| FQ allele positions | 28 | 44 | 83 |
| Wall-clock per sample | n/a | 84 s | 121 s |

(RADAR detects ~3× more AMR genes than the prior diamond_traditional run and is ~3× more sensitive on FQ. v2 vs v1 FQ DB triples the FQ allele coverage.)

---

## Batch Run Logistics

- **Launched**: 2026-05-05 in detached screen session `radar_batch` (`screen -S radar_batch` to attach).
- **Worker logs**: `nicu_resistome_analysis_v2/logs/batch/worker_gpu{0,1,2}.log`
- **Per-sample subprocess logs**: stdout/stderr captured in `worker_gpu*.log`
- **Estimated total time**: 214 samples × 121 s ≈ 7.2 hours per GPU running in parallel.
- **Resume support**: `run_radar_batch.sh` skips samples whose `_gene_stats.tsv` AND `_allele_frequencies.tsv` both already exist, so re-running is safe.

To re-attach:

```bash
screen -r radar_batch          # attach (Ctrl+A D to detach)
# or just tail the per-GPU logs:
tail -f /home/david/projects/ZCH_UCMC_Manuscript/nicu_resistome_analysis_v2/logs/batch/worker_gpu*.log
```
