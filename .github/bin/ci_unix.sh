#!/usr/bin/env bash
set -euo pipefail

# Run two modes:
#  - "standard": default algorithmic params (only CI toggles: --minimal/--no_rooting/--threads)
#  - "tuned": CI-tuned params to exercise overrides

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

rm -rf ci sample_input sample_output mmseqs_db tmp
mkdir -p sample_input mmseqs_db tmp ci/out ci/tmp

cat > sample_input/input1.fa <<'EOF'
>QRY1_A
MKTIIALSYIFCLVFADYKDDDDK
>QRY1_B
MKTVIALSYIFCLVFAEYKDDDDK
>QRY1_C
MKTVIALSYIFCLVFAEYKDDDDE
EOF

cat > sample_input/input2.fa <<'EOF'
>QRY2_A
GAVLIPFWYTSNQDEHKR
>QRY2_B
GAVLIPFWYTSNQDEHKQ
>QRY2_C
GAVLIPFWYTSNQDEHKK
EOF

# Small local MMseqs DB at the *default* location/prefix:
# params.database_dir=mmseqs_db and params.database=UniRef50
cat > tmp/ref.fa <<'EOF'
>REF1_A
MKTIIALSYIFCLVFADYKDDDDK
>REF1_A_var
MKTIIALSYIFCLVFADYKDDDDN
>REF1_B
MKTVIALSYIFCLVFAEYKDDDDK
>REF1_B_var
MKTVIALSYIFCLVFAEYKDDDDN
>REF1_C
MKTVIALSYIFCLVFAEYKDDDDE
>REF1_C_var
MKTVIALSYIFCLVFAEYKDDDDQ
>REF2_A
GAVLIPFWYTSNQDEHKR
>REF2_A_var
GAVLIPFWYTSNQDEHKA
>REF2_B
GAVLIPFWYTSNQDEHKQ
>REF2_B_var
GAVLIPFWYTSNQDEHKG
>REF2_C
GAVLIPFWYTSNQDEHKK
>REF2_C_var
GAVLIPFWYTSNQDEHKT
EOF

mmseqs createdb tmp/ref.fa mmseqs_db/UniRef50
mmseqs createindex mmseqs_db/UniRef50 ci/tmp --threads 2 || true

echo "Versions:"
nextflow -version || true
mafft --version | head -n 1 || true
mmseqs version || true
iqtree2 --version || true
gotree version || true

echo "=== Run: standard (defaults; CI toggles only) ==="
nextflow run ampliphy.nf -profile standard \
  --minimal true \
  --no_rooting true \
  --threads 2

for id in input1 input2; do
  test -f "sample_output/msa/${id}.msa.fa"
  test -f "sample_output/homologs/${id}.homologs.fa"
  test -f "sample_output/msa/${id}.amp.fa"
  test -f "sample_output/tree/${id}.amp.nwk"
  test -f "sample_output/tree/${id}.nwk"
done

echo "=== Run: tuned (CI-tuned overrides) ==="
nextflow run ampliphy.nf -profile standard \
  --minimal true \
  --no_rooting true \
  --threads 2 \
  --input_dir sample_input \
  --output_dir ci/out \
  --custom_database mmseqs_db/UniRef50 \
  --mmseqs_options "-e 10000 --min-seq-id 0.0" \
  --max_depth 1 \
  --max_seqs 10 \
  --iqtree_options "-m JTT+I+G4"

for id in input1 input2; do
  test -f "ci/out/msa/${id}.msa.fa"
  test -f "ci/out/homologs/${id}.homologs.fa"
  test -f "ci/out/msa/${id}.amp.fa"
  test -f "ci/out/tree/${id}.amp.nwk"
  test -f "ci/out/tree/${id}.nwk"
done
