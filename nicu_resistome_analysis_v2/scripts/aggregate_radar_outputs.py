"""Aggregate per-sample RADAR (unified_pipeline) outputs into long-format matrices.

The unified_pipeline writes per-sample TSVs where AMR metadata is packed into the
`protein_id` column as `id|gene_symbol|drug_class|accession`. The named columns
(gene_symbol, gene_family, drug_class, mechanism, gene_type) are empty because
`coverage_calculator` expects an 11-part DIAMOND header that the actual DB doesn't
use. We reconstruct the metadata via:

  1. Parsing the pipe-delimited protein_id (compact: id|gene|drug_class|accession)
  2. Joining `data/amr_protein_db/protein_details.json` for full metadata (id-keyed)
  3. Joining `data/amr_protein_db/amr_hierarchy_annotated.tsv` (accession-keyed)
     for ARO IDs, MEGARes Type/Class/Mechanism/Group, and NCBI hierarchy fields
     (subclass, type, scope, resistance_type).

Inputs (per sample):
  radar_runs/unified/<sid>/<sid>_gene_stats.tsv
  radar_runs/unified/<sid>/<sid>_amr_alignments.tsv  (raw DIAMOND hits, optional)
  radar_runs/unified/<sid>/<sid>_allele_frequencies.tsv
  radar_runs/unified/<sid>/<sid>_protein_matches.tsv
  radar_runs/unified/<sid>/<sid>_sample_stats.txt    (sample-level summary)

Outputs (matched to prior nicu_resistome_analysis schema):
  data/nicu_amr_only.tsv                  long-format gene-level matrix
  data/nicu_fq_allele_frequencies.csv     concatenated FQ allele freq
  data/nicu_fq_protein_matches.csv        concatenated FQ protein matches
  data/sample_summary.tsv                 per-sample summary stats
  data/amr_hierarchy_joined.tsv           debug: per-gene hierarchy lookup
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


def parse_protein_id(pid: str) -> dict:
    """Parse `protein_id` in either of the two formats RADAR emits.

    11-part NCBI form (newer DIAMOND DB, populates _gene_stats.tsv columns natively):
        gene_id|accession|x|y|gene_symbol|gene_family|mechanism|z|drug_class|type|description
    4-part compact form (current biogpu DB; named columns come back empty):
        id|gene_symbol|drug_class|accession
    """
    out = {"protein_id_raw": pid, "id": None, "gene_symbol": None, "drug_class": None, "accession": None}
    if not isinstance(pid, str):
        return out
    parts = pid.split("|")
    if len(parts) >= 11:
        try:
            out["id"] = int(parts[0])
        except ValueError:
            pass
        out["accession"] = parts[1] or None
        out["gene_symbol"] = parts[4] or None
        out["drug_class"] = parts[8] or None
    elif len(parts) >= 4:
        try:
            out["id"] = int(parts[0])
        except ValueError:
            pass
        out["gene_symbol"] = parts[1] or None
        out["drug_class"] = parts[2] or None
        out["accession"] = parts[3] or None
    return out


def load_metadata(metadata_json: Path) -> pd.DataFrame:
    with open(metadata_json) as f:
        items = json.load(f)
    md = pd.DataFrame(items)
    md.rename(columns={"id": "id", "gene_name": "gene_name_meta"}, inplace=True)
    return md


def load_hierarchy(hier_tsv: Path) -> pd.DataFrame:
    return pd.read_csv(hier_tsv, sep="\t")


def collect_amr(amr_root: Path, metadata: pd.DataFrame, hierarchy: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    abundance_frames = []
    summary_rows = []
    for sample_dir in sorted(amr_root.iterdir()):
        if not sample_dir.is_dir() or sample_dir.name.startswith("_"):
            continue
        sid = sample_dir.name
        gs_path = sample_dir / f"{sid}_gene_stats.tsv"
        if not gs_path.exists():
            continue
        df = pd.read_csv(gs_path, sep="\t")
        if df.empty:
            continue

        # Parse protein_id and fill in any blank annotation columns. Prefer the
        # raw column when populated (newer 11-part DB) and fall back to parsed
        # values when blank (current 4-part compact DB).
        parsed = df["protein_id"].apply(parse_protein_id).apply(pd.Series)
        for col in ("id", "gene_symbol", "drug_class", "accession"):
            if col not in df.columns:
                df[col] = parsed[col]
                continue
            raw = df[col]
            is_blank = raw.isna() | (raw.astype(str).str.len() == 0)
            df[col] = raw.where(~is_blank, parsed[col])

        # Filter to detected genes only (was 295 in smoke vs 440 total)
        if "is_detected" in df.columns:
            detected = df[df["is_detected"].astype(str).str.upper() == "TRUE"].copy()
        else:
            detected = df

        # Join metadata for gene_family
        if "id" in detected.columns and not metadata.empty:
            detected = detected.merge(
                metadata[["id", "gene_family", "gene_name_meta"]],
                on="id", how="left", suffixes=("", "_md"),
            )
            # Prefer metadata gene_family if blank in raw
            if "gene_family_md" in detected.columns:
                detected["gene_family"] = detected["gene_family"].where(
                    detected["gene_family"].astype(str).str.len() > 0,
                    detected["gene_family_md"],
                )

        # Join annotated hierarchy on accession for ARO + MEGARes + NCBI fields.
        # accession parsed from protein_id is the precise allele key — gives 1:1
        # match to the annotated TSV's per-accession rows.
        if not hierarchy.empty and "accession" in detected.columns:
            hcols = [c for c in (
                "accession",
                "subclass", "type", "scope", "resistance_type",
                "aro_id", "aro_match_method",
                "megares_type", "megares_class", "megares_mechanism",
                "megares_group", "megares_meg_id", "megares_class_source",
            ) if c in hierarchy.columns]
            hier = hierarchy[hcols].drop_duplicates(subset="accession")
            detected = detected.merge(hier, on="accession", how="left")

        detected.insert(0, "sample_name", sid)
        abundance_frames.append(detected)

        # Per-sample summary
        summary_rows.append(
            {
                "sample_name": sid,
                "total_amr_genes": len(detected),
                "total_amr_reads": int(detected.get("total_reads", pd.Series([0])).sum()),
                "total_amr_rpm": float(detected.get("rpm", pd.Series([0.0])).sum()),
                "high_conf_amr": int((detected.get("detection_confidence", pd.Series([])) == "HIGH").sum()),
                "num_drug_classes": detected["drug_class"].nunique() if "drug_class" in detected.columns else 0,
                "shannon_diversity": _shannon_from_rpm(detected),
            }
        )

    abundance = pd.concat(abundance_frames, ignore_index=True) if abundance_frames else pd.DataFrame()
    summary = pd.DataFrame(summary_rows)
    return abundance, summary


def _shannon_from_rpm(df: pd.DataFrame) -> float:
    import math
    if "rpm" not in df.columns or df.empty:
        return 0.0
    p = df["rpm"].astype(float)
    p = p[p > 0]
    total = p.sum()
    if total <= 0:
        return 0.0
    p = p / total
    return float(-sum(pi * math.log(pi) for pi in p))


def _read_allele_frequencies_csv(path: Path) -> pd.DataFrame:
    """Tolerant reader for allele_frequencies.csv.

    Multi-mutation rows pack unquoted commas into `mutation_summary`
    (e.g. `A460Q(12%),A460D(9%),A460E(78%)`), which trips the C parser. The
    column count is fixed by the header, so we merge any surplus fields back
    into `mutation_summary` (which is always the second-to-last column).
    """
    with open(path) as f:
        header = f.readline().rstrip("\n").split(",")
    n_cols = len(header)
    if "mutation_summary" not in header:
        return pd.read_csv(path)
    ms_idx = header.index("mutation_summary")
    rows = []
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split(",")
            if len(parts) > n_cols:
                n_extra = len(parts) - n_cols
                merged = ",".join(parts[ms_idx:ms_idx + n_extra + 1])
                parts = parts[:ms_idx] + [merged] + parts[ms_idx + n_extra + 1:]
            elif len(parts) < n_cols:
                parts = parts + [""] * (n_cols - len(parts))
            rows.append(parts)
    return pd.DataFrame(rows, columns=header)


def _read_fq_table(sample_dir: Path, sid: str, stem: str) -> pd.DataFrame | None:
    """Locate `<sid>_<stem>.{tsv,csv}` and read it with the matching separator.

    The unified pipeline emits .tsv (newer ZCH builds) or .csv (preterm_arg).
    Both layouts carry the same columns, just different delimiters. The
    allele_frequencies CSV needs the tolerant reader above.
    """
    tsv = sample_dir / f"{sid}_{stem}.tsv"
    csv = sample_dir / f"{sid}_{stem}.csv"
    if tsv.exists():
        return pd.read_csv(tsv, sep="\t")
    if csv.exists():
        if stem == "allele_frequencies":
            return _read_allele_frequencies_csv(csv)
        return pd.read_csv(csv, sep=",")
    return None


def collect_fq(fq_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    af_frames = []
    pm_frames = []
    for sample_dir in sorted(fq_root.iterdir()):
        if not sample_dir.is_dir() or sample_dir.name.startswith("_"):
            continue
        sid = sample_dir.name
        af = _read_fq_table(sample_dir, sid, "allele_frequencies")
        pm = _read_fq_table(sample_dir, sid, "protein_matches")
        if af is not None:
            af_frames.append(af)
        if pm is not None:
            pm_frames.append(pm)

    af = pd.concat(af_frames, ignore_index=True) if af_frames else pd.DataFrame()
    pm = pd.concat(pm_frames, ignore_index=True) if pm_frames else pd.DataFrame()
    return af, pm


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--v2-root", default="/home/david/projects/ZCH_UCMC_Manuscript/nicu_resistome_analysis_v2")
    ap.add_argument("--biogpu-root", default="/home/david/projects/biogpu")
    args = ap.parse_args()

    v2 = Path(args.v2_root)
    biogpu = Path(args.biogpu_root)
    data_dir = v2 / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    metadata = load_metadata(biogpu / "data" / "amr_protein_db" / "protein_details.json")
    hierarchy_path = biogpu / "data" / "amr_protein_db" / "amr_hierarchy_annotated.tsv"
    hierarchy = load_hierarchy(hierarchy_path) if hierarchy_path.exists() else pd.DataFrame()

    radar_root = v2 / "radar_runs" / "unified"
    amr, amr_summary = collect_amr(radar_root, metadata, hierarchy)
    af, pm = collect_fq(radar_root)

    print(f"AMR samples: {amr_summary.shape[0]}, gene rows: {amr.shape[0]}")
    if not af.empty:
        print(f"FQ samples: {af['sample_name'].nunique()}, allele rows: {af.shape[0]}")
    print(f"Protein matches: {pm.shape[0]}")

    if not amr.empty:
        amr.to_csv(data_dir / "nicu_amr_only.tsv", sep="\t", index=False)
    if not amr_summary.empty:
        amr_summary.to_csv(data_dir / "sample_summary.tsv", sep="\t", index=False)
    if not af.empty:
        af.to_csv(data_dir / "nicu_fq_allele_frequencies.csv", index=False)
    if not pm.empty:
        pm.to_csv(data_dir / "nicu_fq_protein_matches.csv", index=False)

    print(f"Wrote outputs to {data_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
