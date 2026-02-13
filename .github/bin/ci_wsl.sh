#!/usr/bin/env bash
# This script is vibe coded.
set -euo pipefail

# Run inside WSL (Ubuntu). Installs Nextflow + micromamba, creates env, runs .github/bin/ci_unix.sh

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
eval "$(micromamba shell hook -s bash)"
if ! micromamba env list | awk '{print $1}' | grep -qx ampliphy; then
  micromamba create -y -n ampliphy -f envs/ampliphy.yml
fi

micromamba run -n ampliphy bash .github/bin/ci_unix.sh
