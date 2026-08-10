# AmpliPhy

<p align="center"><img height="100" alt="AmpliPhy Logo" src="https://github.com/user-attachments/assets/6cd5e374-a50c-4374-b11e-3553c3eaf6bc"/></p>

**_AmpliPhy_** is a scalable, fully automated Nextflow pipeline for improving phylogenetic inference of gene families through database-driven homolog enrichment.

## Publication
Please cite:
> Kim, D., Gil, M., Katoh, K. & Dessimoz, C. (2026). AmpliPhy improves gene trees by adding homologous sequences without affecting alignments. _Bioinformatics Advances_, vbag222. [doi.org/10.1093/bioadv/vbag222](https://doi.org/10.1093/bioadv/vbag222)

## Quick start

Run with conda:

```bash
git clone https://github.com/DessimozLab/ampliphy.git
cd ampliphy

nextflow run ampliphy.nf -profile conda \
  --input_dir sample/input \
  --input_taxonomy sample/taxonomy \
  --output_dir sample/output
```

## Key features

- Automated workflow: MMseqs2 search → MAFFT amplification → IQ-TREE 2 inference → MAD rooting and pruning
- Curated or custom MMseqs2 databases
- Tunable homolog selection and sequence-addition limits
- Optional taxonomy and taxonomic congruence score (TCS) reports
- Portable Nextflow + Bioconda implementation for local or HPC execution

## Workflow overview

For each input gene family, AmpliPhy performs the following steps:

1. Aligns the input protein sequences with MAFFT to generate an initial MSA.
2. Searches for homologous sequences against an MMseqs2 database.
3. Removes exact sequence matches to the original input and selects homologs according to the configured limits.
4. Adds selected homologs to the original MSA with MAFFT `--addfragments --keeplength`.
5. Infers a tree from the amplified MSA with IQ-TREE 2.
6. Roots the amplified tree with MAD unless rooting is disabled.
7. Prunes added homologs from the rooted tree to produce a final tree containing the original input leaves.
8. Optionally generates taxonomy-based reports and TCS comparisons when input taxonomy is supplied.

## Requirements

- Nextflow (DSL2)
- Java
- Bash
- Either:
  - `-profile conda` (recommended), or
  - the following tools available in `PATH`:
    - `mafft`
    - `mmseqs`
    - `iqtree2`
    - `gotree`
    - Python 3 for reporting scripts

## Input sequences

`--input_dir` must contain one or more FASTA files with protein sequences. Each FASTA file represents one gene family.

Recognized extensions, optionally gzipped:

- `.fa`, `.fasta`, `.faa`, `.fna`, `.ffn`, `.frn`
- optional `.gz` suffix

The family identifier is derived from the file name by stripping the recognized FASTA extension and optional `.gz` suffix.

Example:

```bash
input/
├── fam1.fa
├── fam2.faa
└── fam3.fa.gz
```

These files are processed as families `fam1`, `fam2`, and `fam3`.

## Optional input taxonomy

When `--input_taxonomy` is supplied, the directory may contain one taxonomy file per input family:

```bash
taxonomy/
├── fam1.tax
├── fam2.tax
└── fam3.tax
```

Each `<family_id>.tax` file must contain two tab-separated columns:

```bash
sequence_id	species_taxid
```

For example, given `fam1.fa`:

```bash
>S1
AAA
>S2
WWW
```

the corresponding `fam1.tax` may be:

```bash
S1	9606
S2	10090
```

TaxIds must be valid NCBI Taxonomy identifiers at the species rank.

Behavior:

- Families with a valid corresponding `.tax` file are included in taxonomy- and TCS-based reports.
- Families without a corresponding `.tax` file are skipped with a warning.
- A provided `.tax` file that is malformed or contains invalid/non-species TaxIds causes the relevant reporting step to fail.
- If no family has usable taxonomy information, the requested taxonomy-dependent report cannot be generated.

## Output layout

All published outputs are written under `--output_dir`:

```bash
output_dir/
├── homologs/
│   └── {family_id}.homologs.fa
├── msa/
│   ├── {family_id}.msa.fa
│   └── {family_id}.amp.fa
├── tree/
│   ├── {family_id}.amp.nwk
│   └── {family_id}.amp.unpruned.nwk     # only with --keep_unpruned_tree true
└── report/
    ├── homolog_search_report.tsv
    ├── homolog_taxonomy_report.tsv      # when homolog taxonomy is available
    ├── amplified_taxonomy_report.tsv    # when --input_taxonomy is supplied
    └── amplified_tcs_report.tsv         # when --input_taxonomy is supplied
```

### Alignment outputs

- `{family_id}.msa.fa`  
  Initial MAFFT alignment of the input protein sequences.

- `{family_id}.amp.fa`  
  Amplified MSA after adding selected homologs with MAFFT.

### Homolog output

- `{family_id}.homologs.fa`  
  Selected homolog sequences added during MSA amplification.

Exact sequence matches to any input sequence are removed before selection.

### Tree outputs

- `{family_id}.amp.nwk`  
  Final amplified tree after optional MAD rooting and pruning of added homolog leaves. Its tips correspond to the original input sequences.

- `{family_id}.amp.unpruned.nwk`  
  IQ-TREE 2 tree inferred directly from the amplified MSA before pruning. Published only when `--keep_unpruned_tree true` is given.

### Report outputs

- `homolog_search_report.tsv`  
  Reports the number of original sequences, retrieved homologs, selected homologs, and amplified sequences per family.

- `homolog_taxonomy_report.tsv`  
  Summarizes the taxonomy of added homologs, including lineage diversity and least common ancestor information. Generated only when taxonomy metadata is available for the searched MMseqs2 database and taxonomy reporting is not disabled.

- `amplified_taxonomy_report.tsv`  
  Summarizes taxonomy for each original family and its amplified counterpart, using rows such as `fam1` and `fam1.amp`. Generated when `--input_taxonomy` is supplied and the required taxonomy information is available.

- `amplified_tcs_report.tsv`  
  Reports taxonomic congruence scores for original and amplified trees. Generated when `--input_taxonomy` is supplied.

Nextflow intermediate task files remain in the `work/` directory.

## Main options

### General

- `--input_dir`  
  Directory containing input protein FASTA files.  
  Default: `sample_input`

- `--output_dir`  
  Output directory.  
  Default: `sample_output`

- `--threads`  
  Number of threads assigned to supported processes.  
  Default: `4`

- `--max_memory`  
  Memory allocation for MMseqs2-labelled processes.  
  Default: `16 GB`

### MAFFT alignment and amplification

- `--mafft_preset`  
  MAFFT preset used for the initial alignment and amplified MSA construction.  
  Default: `auto`  
  Available values: `auto`, `fast`, `linsi`, `ginsi`, `einsi`

- `--mafft_options`  
  Additional MAFFT options appended to the selected preset.  
  Default: empty

Example:

```bash
nextflow run ampliphy.nf -profile conda \
  --input_dir sample_input \
  --output_dir results \
  --mafft_preset linsi
```

### MMseqs2 homolog search

- `--database`  
  Curated MMseqs2 database to use for homolog searching. Database names are accepted case-insensitively.  
  Default: `UniRef50`

  Supported database choices include:

  - `UniRef100`
  - `UniRef90`
  - `UniRef50`
  - `UniProtKB`
  - `UniProtKB/Swiss-Prot`
  - `NR`

- `--custom_database`  
  Path prefix to an existing local MMseqs2 database. If supplied, this overrides `--database`.  
  Default: empty

- `--database_dir`  
  Directory in which downloaded MMseqs2 databases are cached.  
  Default: `mmseqs_db`

- `--tmp_dir`  
  Temporary directory used by MMseqs2 and auxiliary resources.  
  Default: `./tmp`

- `--mmseqs_options`  
  Additional options passed to `mmseqs easy-search`, for example E-value, sequence identity, or alignment coverage thresholds.  
  Default: empty

- `--max_depth`  
  Relative cap on the number of selected homologs. For an input family containing `N` sequences, at most `max_depth × N` homologs are selected.  
  Default: `5`  
  Use `inf` to disable the depth-based cap.

- `--max_seqs`  
  Absolute cap on the number of selected homologs. If both `--max_depth` and `--max_seqs` impose limits, the smaller limit is applied.  
  Default: `0` (disabled)

Example using a custom database and selection limits:

```bash
nextflow run ampliphy.nf -profile conda \
  --input_dir sample_input \
  --output_dir results \
  --custom_database /path/to/mmseqs_db_prefix \
  --max_depth 3 \
  --max_seqs 500 \
  --mmseqs_options "-e 1e-5 --min-seq-id 0.3"
```

### Taxonomy reporting

- `--input_taxonomy`  
  Directory containing optional `<family_id>.tax` files for the original input sequences. When supplied, enables amplified-taxonomy and TCS reporting for families with valid taxonomy files.  
  Default: empty

- `--ncbi_dir`  
  Directory in which NCBI Taxonomy dump files are stored or downloaded for taxonomy-based reporting.  
  Default: `./ncbi_taxonomy`

- `--no_taxonomy`  
  Disable taxonomy reporting based on MMseqs2 database annotations.  
  Default: `false`

Notes:

- `homolog_taxonomy_report.tsv` requires taxonomy metadata associated with the searched MMseqs2 database.
- `amplified_taxonomy_report.tsv` requires valid input taxonomy files and homolog taxonomy information.
- `amplified_tcs_report.tsv` uses input taxonomy to generate complete lineage files for the original tree tips.

Example:

```bash
nextflow run ampliphy.nf -profile conda \
  --input_dir sample_input \
  --input_taxonomy input_taxonomy \
  --ncbi_dir ncbi_taxonomy \
  --output_dir results
```

### IQ-TREE 2 inference

- `--iqtree_options`  
  Options passed to IQ-TREE 2.  
  Default: `-m JTT+I+G4 -B 1000`

Example:

```bash
nextflow run ampliphy.nf -profile conda \
  --iqtree_options "-m LG+G4 -B 1000"
```

If an MSA contains fewer than four sequences, bootstrap options are disabled for that family because IQ-TREE cannot perform bootstrap analysis on fewer than four sequences.

### Rooting and pruning

- `--no_rooting`  
  Skip MAD rooting and prune directly from the inferred amplified tree.  
  Default: `false`

- `--mad_script`  
  Path to a custom MAD executable.

- `--keep_unpruned_tree`  
  Publish `{family_id}.amp.unpruned.nwk`, the tree inferred from the amplified MSA before rooting/pruning.  
  Default: `false`

## Profiles

- `standard`  
  Uses tools available in the local environment.

- `conda`  
  Enables the Bioconda-based environment defined in `envs/ampliphy.yml`.

## Citation

When using AmpliPhy, please cite the AmpliPhy preprint:

- AmpliPhy preprint: [doi.org/10.64898/2026.01.26.701724](https://www.biorxiv.org/content/10.64898/2026.01.26.701724)

AmpliPhy also relies on the following tools:

- Nextflow
- MMseqs2
- MAFFT
- IQ-TREE 2
- MAD
- gotree
