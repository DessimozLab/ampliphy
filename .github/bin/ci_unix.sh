#!/usr/bin/env bash
set -euo pipefail

# Run two modes:
#  - "standard": default algorithmic params (only CI toggles: --minimal/--no_rooting/--threads)
#  - "tuned": CI-tuned params to exercise overrides

ROOT_DIR="$(git rev-parse --show-toplevel 2>/dev/null || true)"
if [[ -z "${ROOT_DIR}" ]]; then
  ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fi
cd "${ROOT_DIR}"

rm -rf ci sample_input sample_taxonomy sample_output ncbi_taxonomy mmseqs_db tmp
mkdir -p sample_input sample_taxonomy ncbi_taxonomy mmseqs_db tmp ci/out ci/tmp

cat > sample_input/input1.fa <<'EOF'
>QRY1_A
MKTIIALSYIFCLVFADYKDDDDK
>QRY1_B
MKTVIALSYIFCLVFAEYKDDDDK
>QRY1_C
MKTVIALSYIFCLVFAEYKDDDDE
EOF

cat > sample_taxonomy/input1.tax <<'EOF'
QRY1_A 9606
QRY1_B 9598
QRY1_C 9913
EOF

cat > sample_input/input2.fa <<'EOF'
>QRY2_A
GAVLIPFWYTSNQDEHKR
>QRY2_B
GAVLIPFWYTSNQDEHKQ
>QRY2_C
GAVLIPFWYTSNQDEHKK
EOF

cat > sample_taxonomy/input2.tax <<'EOF'
QRY2_A 13616
QRY2_B 9612
QRY2_C 9606
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

cat > tmp/taxonomy.tsv <<'EOF'
REF1_A 9606
REF1_A_var 10090
REF1_B 9598
REF1_B_var 10116
REF1_C 9913
REF1_C_var 13616
REF2_A 13616
REF2_A_var 9598
REF2_B 9544
REF2_B_var 9031
REF2_C 9606
REF2_C_var 10090
EOF

wget -q -O ncbi_taxonomy/taxdump.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf ncbi_taxonomy/taxdump.tar.gz -C ncbi_taxonomy

mmseqs createdb tmp/ref.fa mmseqs_db/UniRef50
mmseqs createtaxdb mmseqs_db/UniRef50 ci/tmp --ncbi-tax-dump ncbi_taxonomy --tax-mapping-file tmp/taxonomy.tsv --threads 2 || true
mmseqs createindex mmseqs_db/UniRef50 ci/tmp --threads 2 || true

echo "Versions:"
nextflow -version || true
mafft --version | head -n 1 || true
mmseqs version || true
iqtree2 --version || true
gotree version || true

echo "=== Run: standard (defaults; CI toggles only) ==="
nextflow run ./ampliphy.nf -profile standard \
  --minimal true \
  --no_rooting true \
  --threads 2

for id in input1 input2; do
  test -f "sample_output/msa/${id}.msa.fa"
  test -f "sample_output/homologs/${id}.homologs.fa"
  test -f "sample_output/msa/${id}.amp.fa"
  test -f "sample_output/tree/${id}.amp.nwk"
done

echo "=== Run: tuned (CI-tuned overrides) ==="
nextflow run ./ampliphy.nf -profile standard \
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
done

echo "=== Run: taxonomy (test NCBI taxonomy parsing) ==="
nextflow run ./ampliphy.nf -profile standard \
  --minimal true \
  --threads 2 \
  --input_dir sample_input \
  --input_taxonomy sample_taxonomy \
  --output_dir ci/out_taxonomy \
  --custom_database mmseqs_db/UniRef50

for id in input1 input2; do
  test -f "ci/out_taxonomy/msa/${id}.msa.fa"
  test -f "ci/out_taxonomy/homologs/${id}.homologs.fa"
  test -f "ci/out_taxonomy/msa/${id}.amp.fa"
  test -f "ci/out_taxonomy/tree/${id}.amp.nwk"
done
test -f "ci/out_taxonomy/report/homolog_search_report.tsv"
test -f "ci/out_taxonomy/report/homolog_taxonomy_report.tsv"
test -f "ci/out_taxonomy/report/amplified_taxonomy_report.tsv"
test -f "ci/out_taxonomy/report/amplified_tcs_report.tsv"
