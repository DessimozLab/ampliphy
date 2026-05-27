#!/usr/bin/env python3

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import Dict, List

from ampliphy_homolog_taxonomy import (
    ensure_taxdump,
    lineage,
    load_taxonomy,
    resolve_taxid,
)

FASTA_SUFFIX = re.compile(
    r"\.(fa|fasta|faa|fna|ffn|frn)(\.gz)?$",
    flags=re.IGNORECASE,
)


def fail(message: str) -> None:
    raise SystemExit(f"Error: {message}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate full NCBI lineage files for AmpliPhy TCS calculation."
    )
    parser.add_argument("--ncbi-dir", required=True, type=Path)
    parser.add_argument("--fasta-dir", required=True, type=Path)
    parser.add_argument("--taxonomy-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def family_id_from_fasta(path: Path) -> str:
    family_id = FASTA_SUFFIX.sub("", path.name)
    if family_id == path.name:
        fail(f"Unrecognized FASTA extension: {path.name}")
    return family_id


def open_fasta(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def read_fasta_ids(path: Path) -> List[str]:
    seq_ids: List[str] = []
    seen = set()

    with open_fasta(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.startswith(">"):
                continue

            header = line[1:].strip()
            seq_id = header.split()[0] if header else ""

            if not seq_id:
                fail(f"{path.name}:{line_number}: empty FASTA identifier.")
            if seq_id in seen:
                fail(f"{path.name}: duplicate FASTA identifier '{seq_id}'.")

            seq_ids.append(seq_id)
            seen.add(seq_id)

    if not seq_ids:
        fail(f"{path.name}: no FASTA records found.")

    return seq_ids


def read_taxonomy_file(
    path: Path,
    expected_ids: List[str],
    parent: Dict[str, str],
    rank: Dict[str, str],
    merged: Dict[str, str],
) -> Dict[str, str]:
    taxids: Dict[str, str] = {}

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                fail(f"{path.name}:{line_number}: blank lines are not permitted.")

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2 or not fields[0] or not fields[1]:
                fail(
                    f"{path.name}:{line_number}: expected exactly two tab-separated "
                    "fields: sequence_id and species_taxid."
                )

            seq_id, raw_taxid = fields

            if seq_id in taxids:
                fail(f"{path.name}: duplicate sequence identifier '{seq_id}'.")
            if not raw_taxid.isdigit():
                fail(f"{path.name}:{line_number}: invalid NCBI TaxId '{raw_taxid}'.")

            resolved_taxid = resolve_taxid(raw_taxid, parent, merged)
            if resolved_taxid is None:
                fail(f"{path.name}:{line_number}: unknown NCBI TaxId '{raw_taxid}'.")
            if rank.get(resolved_taxid) != "species":
                fail(
                    f"{path.name}:{line_number}: TaxId '{raw_taxid}' resolves to "
                    f"rank '{rank.get(resolved_taxid, 'unknown')}', not species."
                )

            taxids[seq_id] = resolved_taxid

    expected = set(expected_ids)
    provided = set(taxids)

    if expected != provided:
        missing = sorted(expected - provided)
        extra = sorted(provided - expected)

        details = []
        if missing:
            details.append(f"missing sequence IDs: {', '.join(missing)}")
        if extra:
            details.append(f"unexpected sequence IDs: {', '.join(extra)}")

        fail(f"{path.name}: taxonomy entries do not match its FASTA file ({'; '.join(details)}).")

    return taxids


def main() -> None:
    args = parse_args()

    nodes_file, names_file, merged_file = ensure_taxdump(args.ncbi_dir)
    parent, rank, _, merged = load_taxonomy(nodes_file, names_file, merged_file)

    fasta_by_family = {
        family_id_from_fasta(path): path
        for path in args.fasta_dir.iterdir()
        if FASTA_SUFFIX.search(path.name)
    }
    taxonomy_by_family = {
        path.name.removesuffix(".tax"): path
        for path in args.taxonomy_dir.glob("*.tax")
    }

    if not fasta_by_family:
        fail("No input FASTA files were staged for lineage generation.")
    if not taxonomy_by_family:
        fail("No <family_id>.tax files were found in --input_taxonomy.")

    missing = sorted(set(fasta_by_family) - set(taxonomy_by_family))
    extra = sorted(set(taxonomy_by_family) - set(fasta_by_family))

    if missing:
        print(
            "Warning: Input taxonomy is unavailable for the following families; "
            "TCS reporting will skip them: " + ", ".join(missing),
            file=sys.stderr,
        )

    if extra:
        print(
            "Warning: Ignoring taxonomy files without corresponding input FASTA files: "
            + ", ".join(extra),
            file=sys.stderr,
        )

    reportable = sorted(set(fasta_by_family) & set(taxonomy_by_family))

    if not reportable:
        fail("No input family has a corresponding <family_id>.tax file; no TCS report can be generated.")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    for family_id in reportable:
        fasta_ids = read_fasta_ids(fasta_by_family[family_id])
        taxids = read_taxonomy_file(
            taxonomy_by_family[family_id],
            fasta_ids,
            parent,
            rank,
            merged,
        )

        output = args.output_dir / f"{family_id}.lineage.tsv"
        with output.open("w") as handle:
            handle.write("query\tlineage\n")

            for seq_id in fasta_ids:
                taxid = taxids[seq_id]
                full_lineage = lineage(taxid, parent)
                formatted = ", ".join(
                    f"{node} ({rank.get(node, 'no rank')})"
                    for node in full_lineage
                )
                handle.write(f"{seq_id}\t{formatted}\n")


if __name__ == "__main__":
    main()