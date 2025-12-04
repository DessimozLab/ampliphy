#!/usr/bin/env bash
set -euo pipefail

ID=""
INPUT=""
DATABASE=""
DB_PATH=""
CUSTOM_DB=""
TMP_ROOT=""
MMSEQS_OPTIONS=""
MAX_DEPTH="5"
MAX_SEQS="0"
THREADS="4"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --id)              ID="$2"; shift 2 ;;
    --input)           INPUT="$2"; shift 2 ;;
    --database)        DATABASE="$2"; shift 2 ;;
    --db-path)         DB_PATH="$2"; shift 2 ;;
    --custom-db)       CUSTOM_DB="$2"; shift 2 ;;
    --tmp-root)        TMP_ROOT="$2"; shift 2 ;;
    --mmseqs-options)  MMSEQS_OPTIONS="$2"; shift 2 ;;
    --max-depth)       MAX_DEPTH="$2"; shift 2 ;;
    --max-seqs)        MAX_SEQS="$2"; shift 2 ;;
    --threads)         THREADS="$2"; shift 2 ;;
    *)
      echo "[ampliphy-mmseqs] Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "${ID}" || -z "${INPUT}" ]]; then
  echo "[ampliphy-mmseqs] ERROR: --id and --input are required" >&2
  exit 1
fi

TMP_ROOT="${TMP_ROOT:-./tmp}"
mkdir -p "${TMP_ROOT}"

TMP_DIR="${TMP_ROOT}/mmseqs_${ID}"

echo "[ampliphy-mmseqs] DATABASE=${DATABASE}" >&2
echo "[ampliphy-mmseqs] DB_PATH=${DB_PATH}" >&2
echo "[ampliphy-mmseqs] CUSTOM_DB=${CUSTOM_DB}" >&2
echo "[ampliphy-mmseqs] TMP_DIR=${TMP_DIR}" >&2

mkdir -p "${TMP_DIR}"

# Choose or validate database
if [[ -n "${CUSTOM_DB}" ]]; then
  echo "Using custom MMseqs2 database at ${CUSTOM_DB}"
  if ! mmseqs view "${CUSTOM_DB}" --id-list 1 > /dev/null 2>&1; then
    echo "Error: Custom database ${CUSTOM_DB} is not a valid MMseqs2 database." >&2
    exit 1
  fi
  DB_TARGET="${CUSTOM_DB}"
else
  echo "Using MMseqs2 database: ${DATABASE}"
  mkdir -p "$(dirname "${DB_PATH}")"
  if [[ ! -e "${DB_PATH}" ]]; then
    echo "Downloading MMseqs2 database ${DATABASE} to ${DB_PATH}"
    mmseqs databases "${DATABASE}" "${DB_PATH}" "${TMP_DIR}"
  else
    echo "MMseqs2 database found at ${DB_PATH}"
  fi
  DB_TARGET="${DB_PATH}"
fi

# Load database into memory
mmseqs touchdb "${DB_TARGET}" || true

RESULT_M8="${ID}.mmseqs.m8"

# Run mmseqs easy-search
mmseqs easy-search "${INPUT}" "${DB_TARGET}" "${RESULT_M8}" "${TMP_DIR}" \
  --threads "${THREADS}" \
  --db-load-mode 2 \
  --format-output target,taln \
  ${MMSEQS_OPTIONS}

# Obtain unique hits and strip gaps
RESULT_UNIQ="${ID}.mmseqs.uniq.m8"
awk '!f[$1]++{ gsub(/-/, "", $2); print $1 "\t" $2 }' "${RESULT_M8}" > "${RESULT_UNIQ}"

# Calculate max sequences to retrieve
if [[ "${INPUT}" == *.gz ]]; then
  SEQ_DEPTH=$(gzip -dc "${INPUT}" | grep -c '^>' || true)
else
  SEQ_DEPTH=$(grep -c '^>' "${INPUT}" || true)
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
  exit 0
fi

# Retrieve homolog sequences
TOTAL_HITS=$(wc -l < "${RESULT_UNIQ}" | tr -d ' ')
EFFECTIVE_HITS="${FINAL_MAX}"
if [[ "${FINAL_MAX}" -gt "${TOTAL_HITS}" ]]; then
  echo "Warning: Requested ${FINAL_MAX} sequences, but only ${TOTAL_HITS} available. Retrieving all available sequences." >&2
  EFFECTIVE_HITS="${TOTAL_HITS}"
fi

head -n "${EFFECTIVE_HITS}" "${RESULT_UNIQ}" \
  | awk '{print ">" $1 "\n" $2}' \
  > "${ID}.homologs.fa"
