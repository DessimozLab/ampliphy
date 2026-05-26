#!/usr/bin/env bash
set -euo pipefail

DATABASE=""
DB_TARGET=""
CUSTOM_DB=""
TMP_ROOT=""
NO_TAXONOMY="false"
TAXONOMY_STATUS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --database)        DATABASE="$2"; shift 2 ;;
    --db-target)       DB_TARGET="$2"; shift 2 ;;
    --custom-db)       CUSTOM_DB="$2"; shift 2 ;;
    --tmp-root)        TMP_ROOT="$2"; shift 2 ;;
    --no-taxonomy)     NO_TAXONOMY="$2"; shift 2 ;;
    --taxonomy-status) TAXONOMY_STATUS="$2"; shift 2 ;;
    *)
      echo "[ampliphy-mmseqs-prepare] Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "${DATABASE}" || -z "${DB_TARGET}" || -z "${TMP_ROOT}" || -z "${TAXONOMY_STATUS}" ]]; then
  echo "[ampliphy-mmseqs-prepare] ERROR: --database, --db-target, --tmp-root and --taxonomy-status are required" >&2
  exit 1
fi

mkdir -p "${TMP_ROOT}"

if [[ -n "${CUSTOM_DB}" ]]; then
  echo "Using custom MMseqs2 database at ${CUSTOM_DB}"
  if ! mmseqs view "${CUSTOM_DB}" --id-list 1 > /dev/null 2>&1; then
    echo "Error: Custom database ${CUSTOM_DB} is not a valid MMseqs2 database." >&2
    exit 1
  fi
  mmseqs touchdb "${CUSTOM_DB}" || true
else
  echo "Using MMseqs2 database: ${DATABASE}"
  mkdir -p "$(dirname "${DB_TARGET}")"
  if [[ ! -e "${DB_TARGET}" ]]; then
    echo "Downloading MMseqs2 database ${DATABASE} to ${DB_TARGET}"
    mmseqs databases "${DATABASE}" "${DB_TARGET}" "${TMP_ROOT}"
  else
    echo "MMseqs2 database found at ${DB_TARGET}"
  fi
  mmseqs touchdb "${DB_TARGET}" || true
fi

if [[ "${NO_TAXONOMY}" == "true" ]]; then
  printf 'disabled\tuser_disabled\n' > "${TAXONOMY_STATUS}"
elif [[ -e "${DB_TARGET}_mapping" && -e "${DB_TARGET}_taxonomy" ]]; then
  printf 'enabled\t%s\n' "${DB_TARGET}" > "${TAXONOMY_STATUS}"
else
  echo "Warning: Taxonomy reporting disabled because taxonomy metadata was not found for ${DB_TARGET}." >&2
  printf 'disabled\ttaxonomy_metadata_unavailable\n' > "${TAXONOMY_STATUS}"
fi