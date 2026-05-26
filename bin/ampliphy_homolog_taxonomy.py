#!/usr/bin/env python3

import argparse
import csv
import shutil
import sys
import tarfile
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

RANK_COLUMNS = [
    ("species", "distinct_species"),
    ("genus", "distinct_genera"),
    ("family", "distinct_families"),
    ("order", "distinct_orders"),
    ("class", "distinct_classes"),
    ("phylum", "distinct_phyla"),
    ("kingdom", "distinct_kingdoms"),
    ("domain", "distinct_domains"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize NCBI taxonomy of homologs added by AmpliPhy."
    )
    parser.add_argument("--ncbi-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("inputs", nargs="+", type=Path)
    return parser.parse_args()


def ensure_taxdump(ncbi_dir: Path) -> Tuple[Path, Path, Optional[Path]]:
    ncbi_dir.mkdir(parents=True, exist_ok=True)

    nodes = ncbi_dir / "nodes.dmp"
    names = ncbi_dir / "names.dmp"
    merged = ncbi_dir / "merged.dmp"

    if nodes.exists() and names.exists():
        return nodes, names, merged if merged.exists() else None

    archive = ncbi_dir / "taxdump.tar.gz"
    partial = ncbi_dir / "taxdump.tar.gz.part"
    wanted = {"nodes.dmp", "names.dmp", "merged.dmp"}

    def validate_archive(path: Path) -> bool:
        try:
            with tarfile.open(path, "r:gz") as handle:
                members = {Path(member.name).name for member in handle.getmembers()}
            return {"nodes.dmp", "names.dmp"}.issubset(members)
        except (tarfile.TarError, EOFError, OSError):
            return False

    def download_archive() -> None:
        partial.unlink(missing_ok=True)
        print(f"Downloading NCBI Taxonomy dump to {archive}", file=sys.stderr)

        try:
            urllib.request.urlretrieve(TAXDUMP_URL, partial)

            if not validate_archive(partial):
                raise RuntimeError(
                    "Downloaded NCBI Taxonomy archive is incomplete or invalid."
                )

            partial.replace(archive)

        except Exception:
            partial.unlink(missing_ok=True)
            raise

    if archive.exists() and not validate_archive(archive):
        print(
            f"Warning: Removing incomplete or corrupted taxonomy archive: {archive}",
            file=sys.stderr,
        )
        archive.unlink()

    if not archive.exists():
        download_archive()

    if not validate_archive(archive):
        archive.unlink(missing_ok=True)
        download_archive()

    with tarfile.open(archive, "r:gz") as handle:
        members = {
            Path(member.name).name: member
            for member in handle.getmembers()
            if Path(member.name).name in wanted
        }

        for name in wanted:
            member = members.get(name)

            if member is None:
                if name == "merged.dmp":
                    continue
                raise RuntimeError(f"{name} was not found in {archive}")

            source = handle.extractfile(member)
            if source is None:
                raise RuntimeError(f"Unable to extract {name} from {archive}")

            with source, (ncbi_dir / name).open("wb") as target:
                shutil.copyfileobj(source, target)

    return nodes, names, merged if merged.exists() else None


def fields_from_dmp(line: str) -> List[str]:
    return [field.strip() for field in line.rstrip("\n").split("|")]


def load_taxonomy(
    nodes_file: Path,
    names_file: Path,
    merged_file: Optional[Path],
) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str], Dict[str, str]]:
    parent: Dict[str, str] = {}
    rank: Dict[str, str] = {}
    name: Dict[str, str] = {}
    merged: Dict[str, str] = {}

    with nodes_file.open() as handle:
        for line in handle:
            fields = fields_from_dmp(line)
            if len(fields) >= 3:
                taxid, parent_taxid, node_rank = fields[:3]
                parent[taxid] = parent_taxid
                rank[taxid] = node_rank

    with names_file.open() as handle:
        for line in handle:
            fields = fields_from_dmp(line)
            if len(fields) >= 4 and fields[3] == "scientific name":
                name[fields[0]] = fields[1]

    if merged_file is not None:
        with merged_file.open() as handle:
            for line in handle:
                fields = fields_from_dmp(line)
                if len(fields) >= 2:
                    merged[fields[0]] = fields[1]

    return parent, rank, name, merged


def resolve_taxid(
    taxid: str,
    parent: Dict[str, str],
    merged: Dict[str, str],
) -> Optional[str]:
    current = taxid
    visited = set()

    while current not in parent and current in merged and current not in visited:
        visited.add(current)
        current = merged[current]

    return current if current in parent else None


def lineage(taxid: str, parent: Dict[str, str]) -> List[str]:
    result: List[str] = []
    current = taxid
    visited = set()

    while current in parent and current not in visited:
        visited.add(current)
        result.append(current)

        next_taxid = parent[current]
        if next_taxid == current:
            break
        current = next_taxid

    return result


def taxon_at_rank(
    lineage_taxids: Iterable[str],
    rank: Dict[str, str],
    target_rank: str,
) -> Optional[str]:
    for taxid in lineage_taxids:
        if rank.get(taxid) == target_rank:
            return taxid
    return None


def find_lca(taxids: List[str], parent: Dict[str, str]) -> Optional[str]:
    if not taxids:
        return None

    lineages = [lineage(taxid, parent) for taxid in taxids]
    common = set(lineages[0])

    for one_lineage in lineages[1:]:
        common.intersection_update(one_lineage)

    for taxid in lineages[0]:
        if taxid in common:
            return taxid

    return None


def family_id_from_path(path: Path) -> str:
    suffix = ".homolog_taxids.tsv"
    return path.name[:-len(suffix)] if path.name.endswith(suffix) else path.stem


def read_selected_taxids(path: Path) -> Tuple[int, List[str]]:
    added = 0
    taxids: List[str] = []

    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue

            added += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[1] not in {"", "0", "NA", "N/A"}:
                taxids.append(fields[1])

    return added, taxids


def main() -> None:
    args = parse_args()
    nodes_file, names_file, merged_file = ensure_taxdump(args.ncbi_dir)
    parent, rank, name, merged = load_taxonomy(nodes_file, names_file, merged_file)

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

    with args.output.open("w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)

        for input_file in sorted(args.inputs, key=lambda path: family_id_from_path(path)):
            family_id = family_id_from_path(input_file)
            homologs_added, raw_taxids = read_selected_taxids(input_file)

            resolved_taxids = [
                resolved
                for raw_taxid in raw_taxids
                if (resolved := resolve_taxid(raw_taxid, parent, merged)) is not None
            ]

            lineages = [lineage(taxid, parent) for taxid in resolved_taxids]

            counts = []
            for target_rank, _ in RANK_COLUMNS:
                observed = {
                    node
                    for one_lineage in lineages
                    if (node := taxon_at_rank(one_lineage, rank, target_rank)) is not None
                }
                counts.append(len(observed))

            lca_taxid = find_lca(resolved_taxids, parent)
            lca_name = name.get(lca_taxid, "") if lca_taxid else ""
            lca_rank = rank.get(lca_taxid, "") if lca_taxid else ""

            writer.writerow([
                family_id,
                homologs_added,
                len(resolved_taxids),
                *counts,
                lca_taxid or "",
                lca_name,
                lca_rank,
            ])


if __name__ == "__main__":
    main()