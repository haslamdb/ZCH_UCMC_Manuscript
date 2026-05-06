#!/bin/bash
# Rename radar_runs/unified/ZJH_N0X_* (single-digit subjects 1-9) to ZJH_NX_*,
# and rename the per-sample files inside each renamed directory so their
# basenames match. Idempotent — already-canonical dirs are skipped.

set -euo pipefail

UNIFIED="/home/david/projects/ZCH_UCMC_Manuscript/nicu_resistome_analysis_v2/radar_runs/unified"

if [[ ! -d "${UNIFIED}" ]]; then
  echo "missing: ${UNIFIED}"
  exit 1
fi

renamed=0
for d in "${UNIFIED}"/ZJH_N0[1-9]_*; do
  [[ -d "$d" ]] || continue
  old_base=$(basename "$d")
  # ZJH_N0X_a_b -> ZJH_NX_a_b
  new_base=$(echo "$old_base" | sed -E 's/^ZJH_N0([1-9])_/ZJH_N\1_/')
  if [[ "$old_base" == "$new_base" ]]; then
    continue
  fi
  new_dir="${UNIFIED}/${new_base}"
  if [[ -e "$new_dir" ]]; then
    echo "skip (target exists): $old_base -> $new_base"
    continue
  fi

  mv "$d" "$new_dir"
  # Rename files inside whose basename starts with the old sample id
  for f in "${new_dir}/${old_base}"*; do
    [[ -e "$f" ]] || continue
    fbase=$(basename "$f")
    new_fbase="${new_base}${fbase#${old_base}}"
    mv "$f" "${new_dir}/${new_fbase}"
  done
  renamed=$((renamed + 1))
done

echo "renamed: ${renamed} directories"
