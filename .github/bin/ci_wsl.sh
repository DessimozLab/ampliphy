#!/usr/bin/env bash
set -euo pipefail

# Run inside WSL (Ubuntu). Installs Nextflow + micromamba, creates env, runs .github/bin/ci_unix.sh

# micromamba root prefix
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/.micromamba}"
mkdir -p "${MAMBA_ROOT_PREFIX}"

apt-get update
apt-get install -y curl ca-certificates openjdk-17-jre-headless bzip2 tar

# Nextflow
if ! command -v nextflow >/dev/null 2>&1; then
  curl -s https://get.nextflow.io | bash
  mv nextflow /usr/local/bin/nextflow
fi

# micromamba
if ! command -v micromamba >/dev/null 2>&1; then
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
  mv bin/micromamba /usr/local/bin/micromamba
fi

# Create env (idempotent)
eval "$(micromamba shell hook -s bash -r "${MAMBA_ROOT_PREFIX}")"

if ! micromamba -r "${MAMBA_ROOT_PREFIX}" env list | awk '{print $1}' | grep -qx ampliphy; then
  micromamba -r "${MAMBA_ROOT_PREFIX}" create -y -n ampliphy -f envs/ampliphy.yml
fi

micromamba -r "${MAMBA_ROOT_PREFIX}" run -n ampliphy bash .github/bin/ci_unix.sh
