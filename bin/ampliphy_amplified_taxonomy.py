#!/usr/bin/env python3

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from ampliphy_homolog_taxonomy import (
    ensure_taxdump,
    find_lca,
    lineage,
    load_taxonomy,
    resolve_taxid,
    taxon_at_rank,
)

FASTA_SUFFIX = re.compile(
    r"\.(fa|fasta|faa|fna|ffn|frn)(\.gz)?$",
    flags=re.IGNORECASE,
)

RANK_COLUMNS = [
    ("species", "distinct_species"),
    ("genus", "distinct_genera"),
    ("family", "distinct_families"),
    ("order", "distinct_orders"),
    ("class", "distinct_classes"),
    ("phylum", "distinct_phyla"),
    ("kingdom", "distinct_kingdoms"),
    ("superkingdom", "distinct_domains"),
]


def fail(message: str) -> None:
    raise SystemExit(f"Error: {message}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Report taxonomy of original and amplified AmpliPhy families."
    )
    parser.add_argument("--ncbi-dir", required=True, type=Path)
    parser.add_argument("--fasta-dir", required=True, type=Path)
    parser.add_argument("--input-taxonomy-dir", required=True, type=Path)
    parser.add_argument("--homolog-taxids-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
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

            seq_id = line[1:].strip().split()[0] if line[1:].strip() else ""
            if not seq_id:
                fail(f"{path.name}:{line_number}: empty FASTA identifier.")
            if seq_id in seen:
                fail(f"{path.name}: duplicate FASTA identifier '{seq_id}'.")

            seq_ids.append(seq_id)
            seen.add(seq_id)

    if not seq_ids:
        fail(f"{path.name}: no FASTA records found.")

    return seq_ids


def read_input_taxonomy(
    path: Path,
    expected_seq_ids: Sequence[str],
    parent: Dict[str, str],
    rank: Dict[str, str],
    merged: Dict[str, str],
) -> List[str]:
    seq_ids: List[str] = []
    species_taxids: List[str] = []
    seen = set()

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

            if seq_id in seen:
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

            seq_ids.append(seq_id)
            species_taxids.append(resolved_taxid)
            seen.add(seq_id)

    if seq_ids != list(expected_seq_ids):
        missing = sorted(set(expected_seq_ids) - set(seq_ids))
        extra = sorted(set(seq_ids) - set(expected_seq_ids))

        details = []
        if missing:
            details.append(f"missing IDs: {', '.join(missing)}")
        if extra:
            details.append(f"unexpected IDs: {', '.join(extra)}")
        if not missing and not extra:
            details.append("sequence order does not match the FASTA file")

        fail(f"{path.name}: taxonomy entries do not match the FASTA file ({'; '.join(details)}).")

    return species_taxids


def read_homolog_taxids(
    path: Path,
    parent: Dict[str, str],
    merged: Dict[str, str],
) -> Tuple[int, List[str]]:
    sequence_count = 0
    resolved_taxids: List[str] = []

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                fail(
                    f"{path.name}:{line_number}: expected two tab-separated "
                    "fields: homolog_id and taxid."
                )

            sequence_count += 1
            raw_taxid = fields[1]

            if not raw_taxid or raw_taxid in {"0", "NA", "N/A"}:
                continue

            resolved_taxid = resolve_taxid(raw_taxid, parent, merged)
            if resolved_taxid is None:
                print(
                    f"Warning: {path.name}:{line_number}: ignoring unresolved "
                    f"homolog TaxId '{raw_taxid}'.",
                    file=sys.stderr,
                )
                continue

            resolved_taxids.append(resolved_taxid)

    return sequence_count, resolved_taxids


def summarize(
    family_id: str,
    sequence_count: int,
    taxids: List[str],
    parent: Dict[str, str],
    rank: Dict[str, str],
    name: Dict[str, str],
) -> List[object]:
    lineages = [lineage(taxid, parent) for taxid in taxids]

    rank_counts = []
    for target_rank, _ in RANK_COLUMNS:
        observed = {
            node
            for one_lineage in lineages
            if (node := taxon_at_rank(one_lineage, rank, target_rank)) is not None
        }
        rank_counts.append(len(observed))

    lca_taxid = find_lca(taxids, parent)
    lca_name = name.get(lca_taxid, "") if lca_taxid else ""
    lca_rank = rank.get(lca_taxid, "") if lca_taxid else ""

    return [
        family_id,
        sequence_count,
        len(taxids),
        *rank_counts,
        lca_taxid or "",
        lca_name,
        lca_rank,
    ]


def main() -> None:
    args = parse_args()

    nodes_file, names_file, merged_file = ensure_taxdump(args.ncbi_dir)
    parent, rank, name, merged = load_taxonomy(nodes_file, names_file, merged_file)

    fasta_files = sorted(
        path for path in args.fasta_dir.iterdir() if FASTA_SUFFIX.search(path.name)
    )
    if not fasta_files:
        fail("No input FASTA files were staged for amplified taxonomy reporting.")

    fasta_by_family = {
        family_id_from_fasta(path): path
        for path in fasta_files
    }
    taxonomy_by_family = {
        path.name[:-4]: path
        for path in args.input_taxonomy_dir.glob("*.tax")
    }
    homolog_by_family = {
        path.name.removesuffix(".homolog_taxids.tsv"): path
        for path in args.homolog_taxids_dir.glob("*.homolog_taxids.tsv")
    }

    if not taxonomy_by_family:
        fail("No <family_id>.tax files were found in --input_taxonomy.")

    expected = set(fasta_by_family)
    provided = set(taxonomy_by_family)
    homolog_available = set(homolog_by_family)

    missing_input_taxonomy = sorted(expected - provided)
    extra_input_taxonomy = sorted(provided - expected)

    if missing_input_taxonomy:
        print(
            "Warning: Input taxonomy is unavailable for the following families; "
            "they will be skipped: "
            + ", ".join(missing_input_taxonomy),
            file=sys.stderr,
        )

    if extra_input_taxonomy:
        print(
            "Warning: Ignoring taxonomy files without corresponding input FASTA files: "
            + ", ".join(extra_input_taxonomy),
            file=sys.stderr,
        )

    missing_homolog_taxonomy = sorted((expected & provided) - homolog_available)
    if missing_homolog_taxonomy:
        print(
            "Warning: Homolog taxonomy is unavailable for the following annotated "
            "families; they will be skipped: "
            + ", ".join(missing_homolog_taxonomy),
            file=sys.stderr,
        )

    reportable_families = sorted(expected & provided & homolog_available)

    if not reportable_families:
        fail(
            "No input family has both input taxonomy annotations and homolog "
            "taxonomy data. No amplified taxonomy report can be generated."
        )

    header = [
        "family_id",
        "sequence_count",
        "sequences_with_taxid",
        *[column for _, column in RANK_COLUMNS],
        "lca_taxid",
        "lca_name",
        "lca_rank",
    ]

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)

        for family_id in reportable_families:
            fasta_ids = read_fasta_ids(fasta_by_family[family_id])

            input_taxids = read_input_taxonomy(
                taxonomy_by_family[family_id],
                fasta_ids,
                parent,
                rank,
                merged,
            )

            homolog_count, homolog_taxids = read_homolog_taxids(
                homolog_by_family[family_id],
                parent,
                merged,
            )

            writer.writerow(
                summarize(
                    f"{family_id}.ori",
                    len(fasta_ids),
                    input_taxids,
                    parent,
                    rank,
                    name,
                )
            )

            writer.writerow(
                summarize(
                    f"{family_id}.amp",
                    len(fasta_ids) + homolog_count,
                    input_taxids + homolog_taxids,
                    parent,
                    rank,
                    name,
                )
            )


if __name__ == "__main__":
    main()