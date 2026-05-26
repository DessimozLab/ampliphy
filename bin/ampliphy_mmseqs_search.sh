#!/usr/bin/env bash
set -euo pipefail

ID=""
INPUT=""
DB_TARGET=""
TMP_ROOT=""
MMSEQS_OPTIONS=""
MAX_DEPTH="5"
MAX_SEQS="0"
THREADS="4"
TAXONOMY_STATUS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --id)              ID="$2"; shift 2 ;;
    --input)           INPUT="$2"; shift 2 ;;
    --db-target)       DB_TARGET="$2"; shift 2 ;;
    --tmp-root)        TMP_ROOT="$2"; shift 2 ;;
    --mmseqs-options)  MMSEQS_OPTIONS="$2"; shift 2 ;;
    --max-depth)       MAX_DEPTH="$2"; shift 2 ;;
    --max-seqs)        MAX_SEQS="$2"; shift 2 ;;
    --threads)         THREADS="$2"; shift 2 ;;
    --taxonomy-status) TAXONOMY_STATUS="$2"; shift 2 ;;
    *)
      echo "[ampliphy-mmseqs-search] Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "${ID}" || -z "${INPUT}" || -z "${DB_TARGET}" || -z "${TMP_ROOT}" || -z "${TAXONOMY_STATUS}" ]]; then
  echo "[ampliphy-mmseqs-search] ERROR: --id, --input, --db-target, --tmp-root and --taxonomy-status are required" >&2
  exit 1
fi

TMP_DIR="${TMP_ROOT}/mmseqs_${ID}"
mkdir -p "${TMP_DIR}"

# Make a local uncompressed copy of the input for depth + identity checks
LOCAL_FA="${ID}.mmseqs.input.fa"
if [[ "${INPUT}" == *.gz ]]; then
  gzip -dc "${INPUT}" > "${LOCAL_FA}"
else
  cp "${INPUT}" "${LOCAL_FA}"
fi

# Initialize report with header
REPORT="${ID}.homolog_search_stats.tsv"

SEQ_DEPTH=$(grep -c '^>' "${LOCAL_FA}" || true)

write_report() {
  local homologs_found="$1"
  local homologs_added="$2"
  local amplified_sequences=$((SEQ_DEPTH + homologs_added))

  printf '%s\t%s\t%s\t%s\t%s\n' \
    "${ID}" \
    "${SEQ_DEPTH}" \
    "${homologs_found}" \
    "${homologs_added}" \
    "${amplified_sequences}" \
    > "${REPORT}"
}

RESULT_M8="${ID}.mmseqs.m8"

TAXONOMY_ENABLED="false"
if [[ "$(cut -f1 "${TAXONOMY_STATUS}")" == "enabled" ]]; then
  TAXONOMY_ENABLED="true"
fi

# Run MMseqs2 search against the prepared database
if [[ "${TAXONOMY_ENABLED}" == "true" ]]; then
  FORMAT_OUTPUT="target,taln,taxid"
else
  FORMAT_OUTPUT="target,taln"
fi
mmseqs easy-search "${INPUT}" "${DB_TARGET}" "${RESULT_M8}" "${TMP_DIR}" \
  --threads "${THREADS}" \
  --db-load-mode 2 \
  --format-output ${FORMAT_OUTPUT} \
  ${MMSEQS_OPTIONS}

# Unique hits by target, strip gaps from alignment sequence
RESULT_UNIQ="${ID}.mmseqs.uniq.m8"
if [[ "${TAXONOMY_ENABLED}" == "true" ]]; then
  awk '!f[$1]++{ gsub(/-/, "", $2); print $1 "\t" $2 "\t" $3 }' "${RESULT_M8}" > "${RESULT_UNIQ}"
else
  awk '!f[$1]++{ gsub(/-/, "", $2); print $1 "\t" $2 }' "${RESULT_M8}" > "${RESULT_UNIQ}"
fi

# Build canonical input sequences (one per record, concatenated, letters only)
INPUT_SEQ="${ID}.mmseqs.input.seq"
awk '
  /^>/ {
    if (seq != "") {
      gsub(/[^A-Za-z]/,"",seq); print seq; seq="";
    }
    next
  }
  { seq = seq $0 }
  END {
    if (seq != "") {
      gsub(/[^A-Za-z]/,"",seq); print seq;
    }
  }
' "${LOCAL_FA}" | sort -u > "${INPUT_SEQ}"

# Remove hits whose sequence is exactly identical to any input sequence
RESULT_FILTERED="${ID}.mmseqs.uniq.filtered.m8"
RESULT_FILTERED="${ID}.mmseqs.uniq.filtered.m8"
awk '
  FNR==NR { seen[$1]=1; next }
  {
    seq = $2;
    if (!(seq in seen))
      print $0;
  }
' "${INPUT_SEQ}" "${RESULT_UNIQ}" > "${RESULT_FILTERED}"

FILTERED_COUNT=$(wc -l < "${RESULT_FILTERED}" | tr -d ' ')

if [[ "${FILTERED_COUNT}" -eq 0 ]]; then
  echo "Warning: All MMseqs2 hits are identical to the input sequences; no non-identical homologs remain." >&2
  : > "${ID}.homologs.fa"
  if [[ "${TAXONOMY_ENABLED}" == "true" ]]; then
    : > "${ID}.homolog_taxids.tsv"
  fi
  write_report 0 0
  exit 0
fi

if [[ "${MAX_DEPTH}" == "inf" ]]; then
  MAX_FROM_DEPTH=2147483647
else
  MAX_FROM_DEPTH=$(awk -v d="${SEQ_DEPTH}" -v m="${MAX_DEPTH}" 'BEGIN { n=d*m; if (n<0) n=0; printf "%d\n", int(n) }')
fi

if [[ -n "${MAX_SEQS}" && "${MAX_SEQS}" != "0" ]]; then
  if [[ "${MAX_FROM_DEPTH}" -eq 0 || "${MAX_SEQS}" -lt "${MAX_FROM_DEPTH}" ]]; then
    FINAL_MAX="${MAX_SEQS}"
  else
    FINAL_MAX="${MAX_FROM_DEPTH}"
  fi
else
  FINAL_MAX="${MAX_FROM_DEPTH}"
fi

if [[ "${FINAL_MAX}" -le 0 ]]; then
  echo "Warning: Calculated maximum sequences to retrieve is ${FINAL_MAX}. No homologs will be retrieved." >&2
  : > "${ID}.homologs.fa"
  if [[ "${TAXONOMY_ENABLED}" == "true" ]]; then
    : > "${ID}.homolog_taxids.tsv"
  fi
  write_report "${FILTERED_COUNT}" 0
  exit 0
fi

EFFECTIVE_HITS="${FINAL_MAX}"
if [[ "${FINAL_MAX}" -gt "${FILTERED_COUNT}" ]]; then
  echo "Warning: Requested ${FINAL_MAX} sequences, but only ${FILTERED_COUNT} available. Retrieving all available sequences." >&2
  EFFECTIVE_HITS="${FILTERED_COUNT}"
fi

RESULT_SELECTED="${ID}.mmseqs.selected.m8"
head -n "${EFFECTIVE_HITS}" "${RESULT_FILTERED}" > "${RESULT_SELECTED}"
awk '{print ">" $1 "_enr\n" $2}' "${RESULT_SELECTED}" > "${ID}.homologs.fa"
if [[ "${TAXONOMY_ENABLED}" == "true" ]]; then
  awk 'BEGIN{OFS="\t"} {print $1, $3}' "${RESULT_SELECTED}" > "${ID}.homolog_taxids.tsv"
fi
write_report "${FILTERED_COUNT}" "${EFFECTIVE_HITS}"