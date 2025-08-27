#!/usr/bin/env bash
# Add a "file_size_restored" column by re-checking each path in the "bam_files" column.
# Parallelized. Does NOT modify BAMs. Missing files => "file_missing".
# Usage: script INPUT_TSV
# Tunables: MAX_JOBS (default: nproc), UNITS={bytes,KiB,MiB,GiB} (default: MiB), DECIMALS (default: 0)

INPUT="/Data1/git/meyer-nanopore/scripts/bam_input_metadata_8_18_2025_COND.txt"
OUTPUT="/Data1/git/meyer-nanopore/scripts/bam_input_metadata_8_18_2025_COND_estim.txt"

[[ -f "${INPUT}" ]] || { echo "error: ${INPUT} not found"; exit 1; }

MAX_JOBS="${MAX_JOBS:-$(command -v nproc >/dev/null 2>&1 && nproc || echo 4)}"
UNITS="${UNITS:-MiB}"      # bytes | KiB | MiB | GiB
DECIMALS="${DECIMALS:-0}"  # decimals in output

tmpdir="$(mktemp -d)"
cleanup(){ rm -rf "${tmpdir}"; }
trap cleanup EXIT

# ── Normalize line endings to LF into a temp copy ───────────────────
SAN="${tmpdir}/input_lf.tsv"
sed 's/\r$//' "${INPUT}" > "${SAN}"

# ── Find "bam_files" column index ───────────────────────────────────
bam_col_idx="$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="bam_files"){print i; exit}}' "${SAN}")"
[[ -n "${bam_col_idx}" ]] || { echo "error: header lacks bam_files column"; exit 1; }

# ── Write header with appended column ───────────────────────────────
awk -F'\t' 'NR==1{OFS="\t"; print $0,"file_size_restored"}' "${SAN}" > "${OUTPUT}"

# ── Build tasks: [lineno \t path] for data lines ────────────────────
awk -v FS='\t' -v c="${bam_col_idx}" 'NR>1{p=$c; print NR-1 "\t" p}' "${SAN}" > "${tmpdir}/tasks.tsv"
line_count="$(wc -l < "${tmpdir}/tasks.tsv" | tr -d ' ')"

# ── Unit conversion helper ──────────────────────────────────────────
convert_bytes() {
  local bytes="$1" units="$2" decimals="$3"
  awk -v b="${bytes}" -v u="${units}" -v d="${decimals}" 'BEGIN{
    scale=1
    if(u=="KiB") scale=1024
    else if(u=="MiB") scale=1024*1024
    else if(u=="GiB") scale=1024*1024*1024
    val=b/scale
    printf "%.*f", d, val
  }'
}

# ── Worker: compute size or file_missing ────────────────────────────
worker() {
  local lineno="$1" path="$2" outdir="$3" units="$4" decimals="$5"
  local result="file_missing"
  if [[ -f "$path" ]]; then
    if size_bytes=$(stat -c '%s' -- "$path" 2>/dev/null); then :; else
      size_bytes=$(stat -f '%z' -- "$path" 2>/dev/null || echo "")
    fi
    if [[ -n "${size_bytes}" ]]; then
      result="$(convert_bytes "${size_bytes}" "${units}" "${decimals}")"
    fi
  fi
  printf "%s\n" "${result}" > "${outdir}/${lineno}"
}

# ── Launch workers with throttle ────────────────────────────────────
while IFS=$'\t' read -r lineno fpath; do
  ( worker "${lineno}" "${fpath}" "${tmpdir}" "${UNITS}" "${DECIMALS}" ) &
  while (( $(jobs -r -p | wc -l) >= MAX_JOBS )); do wait -n || true; done
done < "${tmpdir}/tasks.tsv"
wait

# ── Collect results in-order ────────────────────────────────────────
sizes="${tmpdir}/sizes.txt"; : > "${sizes}"
for i in $(seq 1 "${line_count}"); do
  if [[ -f "${tmpdir}/${i}" ]]; then cat "${tmpdir}/${i}" >> "${sizes}"
  else echo "file_missing" >> "${sizes}"; fi
done

# ── Append column to sanitized data lines ───────────────────────────
tail -n +2 "${SAN}" | paste -d $'\t' - "${sizes}" >> "${OUTPUT}"

echo "Wrote: ${OUTPUT}"
echo "Units=${UNITS} Decimals=${DECIMALS} Jobs=${MAX_JOBS}"