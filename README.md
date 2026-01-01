# AmpliPhy

**_AmpliPhy_** is a scalable, fully automated Nextflow pipeline for improving phylogenetic inference through database-driven homolog enrichment. 

## Quick start

Run with conda:

```bash
git clone https://github.com/DessimozLab/ampliphy.git
nextflow run ampliphy.nf -profile conda \
  --input_dir sample_input \
  --output_dir results
```

## Key features

- Fully automated pipeline: MMseqs2 search → MAFFT MSA + amplification → IQ-TREE 2 inference → optional root/prune  
- Curated or custom MMseqs2 databases to search for homologs
- Tunable homolog selection: MMseqs2 thresholds + relative/absolute caps on added sequences  
- Scalable Nextflow + Bioconda implementation; runs locally or on HPC

## Requirements

- Nextflow (DSL2)
- Either:
  - `-profile conda` (recommended), or
  - Following tools available in `PATH`:
    - `mafft`, `mmseqs2`, `iqtree`, `gotree`
    - `mad` ([link](https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen/mad2-2.zip); configured via `--mad_script`)

## Input

`--input_dir` must contain one or more FASTA files with protein sequences.

Recognized extensions (optionally gzipped):

- `.fa`, `.fasta`, `.faa`, `.fna`, `.ffn`, `.frn`
- with optional `.gz` suffix

The pipeline derives `id` from the file name by stripping the extension (and `.gz`).


## Output layout

All outputs are published under `--output_dir`:

- `output_dir/msa/`
  - `{id}.msa.fa`   (initial MAFFT MSA)
  - `{id}.amp.fa`   (amplified MSA after adding homologs)

- `output_dir/homologs/`
  - `{id}.homologs.fa` (homolog sequences selected from MMseqs2 hits)

- `output_dir/tree/`
  - `{id}.amp.nwk` (IQ-TREE 2 tree from amplified MSA)
  - `{id}.nwk`     (final rooted+pruned tree)

Nextflow’s `work/` directory contains intermediate task folders.

## Main options

### General

- `--input_dir` (default: `sample_input`)
- `--output_dir` (default: `sample_output`)
- `--threads` (default: `4`)
- `--max_memory` (default: `16 GB`)  
  Used by MMseqs2-labeled processes (`withLabel: 'mmseqs'`).

### MAFFT (alignment + amplification)

- `--mafft_preset` (default: `auto`)  
  One of: `auto`, `fast`, `linsi`, `ginsi`, `einsi`
- `--mafft_options` (default: empty)  
  Additional MAFFT flags appended to the preset.

### MMseqs2

- `--database` (default: `UniRef50`)  
  Database name (case-insensitive aliases supported in `mmseqs_prepare_db`).
- `--custom_database` (default: empty)  
  Path prefix to a local MMseqs2 database. If set, `--database` is ignored.
- `--database_dir` (default: `mmseqs_db`)  
  Where downloaded MMseqs2 databases are cached.
- `--tmp_dir` (default: `./tmp`)  
  Temp directory used for MMseqs2 work dirs.
- `--mmseqs_options` (default: empty)  
  Extra arguments appended to `mmseqs easy-search`.
- `--max_depth` (default: `5`)  
  Relative cap on the number of homologs: at most `max_depth * N`, where `N` is the number of sequences in the input file.  
  Use `inf` for no depth-based cap.
- `--max_seqs` (default: `0`)  
  Absolute cap on homolog count. If both caps are set, the smaller wins.  
  Use `0` to disable the absolute cap.

Notes:
- Homolog hits are deduplicated by target ID.
- Hits whose sequence is **exactly identical** to any input sequence are removed before slicing. If none remain, an empty `{id}.homologs.fa` is produced and the pipeline continues.

### IQ-TREE 2

- `--iqtree_options` (default: `-m JTT+I+G4 -B 1000`)  
  Passed directly to `iqtree2`.

### Rooting + pruning

- `--no_rooting` (default: `false`)  
  If true, rooting is skipped and the amplified tree is used as-is for pruning.
- `--mad_script` (default: `bin/mad.py`)  
  Path to MAD executable/script used for rooting.

Example (use platform-specific MAD executable):

```bash
nextflow run ampliphy.nf -profile conda \
  --mad_script bin/mad.osx
```

---

## Profiles

- `standard` (default): local executor
- `conda`: enables conda and uses `envs/ampliphy.yml`

Example:

```bash
nextflow run ampliphy.nf -profile conda
```

---

## Example runs

Use a custom MMseqs2 database:

```bash
nextflow run ampliphy.nf -profile conda \
  --custom_database /path/to/mmseqs_db_prefix \
  --input_dir sample_input \
  --output_dir results
```

Tune MMseqs2 hit limits:

```bash
nextflow run ampliphy.nf -profile conda \
  --max_depth 3 \
  --max_seqs 500
```

Change MAFFT preset:

```bash
nextflow run ampliphy.nf -profile conda \
  --mafft_preset linsi
```

Change IQ-TREE model/options:

```bash
nextflow run ampliphy.nf -profile conda \
  --iqtree_options "-m LG+G4 -B 1000"
```

Disable rooting:

```bash
nextflow run ampliphy.nf -profile conda \
  --no_rooting true
```
