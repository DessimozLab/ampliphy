#!/usr/bin/env -S python3 -u
"""
Usage: tcs_overlap.py <tree.nwk> <lineage.tsv>

Calculates the taxonomy overlap score over the NCBI taxonomy tree of the phylogenetic tree.
  - tree.nwk: A Newick formatted tree file, where each leaf node is labeled with a Uniprot ID.
    Example: ((P12345, P67890), (Q12345, Q67890));

  - lineage.tsv: A tab-separated file containing the lineage information for each Uniprot ID.
    Example:
    query    lineage
    P03007   3379134 (kingdom), 1224 (phylum), 1236 (class), 91347 (order), 543 (family), 561 (genus), 562 (species)

Original source code written by David Moi (@cactuskid)
Modified and written by Dongwook Kim (@endixk)

Please cite the following paper if you use this code:
Moi, D. et al. (2023). bioRxiv, https://doi.org/10.1101/2023.09.19.558401
"""

import csv
import sys
from pathlib import Path

VERBOSE = False


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class TreeNode:
    """Minimal tree node providing the interface required by the scoring logic."""

    def __init__(self, name="", children=None):
        self.name = name
        self.children = children or []

    def add_feature(self, name, value):
        setattr(self, name, value)

    def is_leaf(self):
        return not self.children

    def get_children(self):
        return self.children

    def iter_leaves(self):
        if self.is_leaf():
            yield self
        else:
            for child in self.children:
                yield from child.iter_leaves()


class Tree:
    """Minimal wrapper retaining the previous tree.treenode access pattern."""

    def __init__(self, treenode):
        self.treenode = treenode


class NewickParser:
    """Parse the Newick features used by this script without external packages."""

    def __init__(self, text):
        self.text = text
        self.pos = 0
        self.length = len(text)

    def parse(self):
        self._skip_whitespace_and_comments()
        if self.pos >= self.length:
            raise ValueError("Newick input is empty.")

        root = self._parse_subtree()
        self._skip_whitespace_and_comments()

        if self._peek() == ";":
            self.pos += 1
            self._skip_whitespace_and_comments()

        if self.pos != self.length:
            raise ValueError(
                f"Unexpected content in Newick tree at position {self.pos}: "
                f"{self.text[self.pos:self.pos + 20]!r}"
            )
        return Tree(root)

    def _parse_subtree(self):
        self._skip_whitespace()
        if self._peek() == "(":
            self.pos += 1
            children = [self._parse_subtree()]
            while True:
                self._skip_whitespace()
                char = self._peek()
                if char == ",":
                    self.pos += 1
                    children.append(self._parse_subtree())
                elif char == ")":
                    self.pos += 1
                    break
                else:
                    raise ValueError(
                        f"Expected ',' or ')' in Newick tree at position {self.pos}."
                    )
            node = TreeNode(children=children)
            node.name = self._parse_label()
        else:
            label = self._parse_label()
            node = TreeNode(name=label)

        self._parse_suffix()
        return node

    def _parse_label(self):
        self._skip_whitespace()
        if self._peek() == "'":
            self.pos += 1
            chars = []
            while self.pos < self.length:
                char = self.text[self.pos]
                self.pos += 1
                if char == "'":
                    if self._peek() == "'":
                        chars.append("'")
                        self.pos += 1
                    else:
                        return "".join(chars)
                else:
                    chars.append(char)
            raise ValueError("Unterminated quoted Newick label.")

        start = self.pos
        while self.pos < self.length:
            char = self.text[self.pos]
            if char in ":,();[]" or char.isspace():
                break
            self.pos += 1
        return self.text[start:self.pos]

    def _parse_suffix(self):
        while True:
            self._skip_whitespace()
            char = self._peek()
            if char == "[":
                self._skip_comment()
            elif char == ":":
                self.pos += 1
                self._skip_branch_length()
            else:
                return

    def _skip_branch_length(self):
        self._skip_whitespace()
        start = self.pos
        while self.pos < self.length:
            char = self.text[self.pos]
            if char in ",();[" or char.isspace():
                break
            self.pos += 1
        if start == self.pos:
            raise ValueError(f"Missing branch length at position {self.pos}.")

    def _skip_whitespace_and_comments(self):
        while True:
            self._skip_whitespace()
            if self._peek() == "[":
                self._skip_comment()
            else:
                return

    def _skip_whitespace(self):
        while self.pos < self.length and self.text[self.pos].isspace():
            self.pos += 1

    def _skip_comment(self):
        if self._peek() != "[":
            return
        depth = 0
        while self.pos < self.length:
            char = self.text[self.pos]
            self.pos += 1
            if char == "[":
                depth += 1
            elif char == "]":
                depth -= 1
                if depth == 0:
                    return
        raise ValueError("Unterminated Newick comment.")

    def _peek(self):
        if self.pos >= self.length:
            return ""
        return self.text[self.pos]


def parse_newick(text):
    return NewickParser(text).parse()


def read_lineage_tsv(lineage_file):
    try:
        with open(lineage_file, "r", newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                eprint("Error: lineage TSV is empty.")
                sys.exit(1)

            required_columns = {"query", "lineage"}
            if not required_columns.issubset(reader.fieldnames):
                eprint("Error: lineage TSV must include 'query' and 'lineage' columns.")
                sys.exit(1)

            rows = list(reader)
    except (OSError, csv.Error) as exc:
        eprint(f"Error reading lineage file: {exc}")
        sys.exit(1)

    return rows


def make_lineages(uniprot_rows):
    return {
        row["query"]: set(row["lineage"].split(","))
        for row in uniprot_rows
    }


def get_species(uniprot_rows):
    return {
        row["query"]: row["lineage"].split(",")[-1].split("(")[0].strip()
        for row in uniprot_rows
    }


def label_leaves( tree , leaf_lineages , species_map):
    """
    Adds lineage information to the leaves of a tree.

    Parameters:
    tree (Tree): A tree object containing a ``treenode`` root.
    leaf_lineages (dict): A dictionary mapping leaf names to lineage information.

    Returns:
    Tree: The input tree object with the added lineage information.
    """
    species_count = {}
    #takes lineage records with lineage info from uniprot
    for n in tree.treenode.iter_leaves():
        if n.name in leaf_lineages:
            n.add_feature( 'lineage' ,   leaf_lineages[n.name] )
            if species_map[n.name] in species_count:
                species_count[species_map[n.name]] += 1
            else:
                species_count[species_map[n.name]] = 1
            spcount = species_count[species_map[n.name]]
            #pad so that it is a 3 digit number
            spcount = str(spcount).zfill(3)
            n.add_feature( 'sp_num' ,   species_map[n.name]+ '_' + spcount)
        else:
            n.add_feature( 'lineage' ,   None )
            n.add_feature( 'sp_num' ,   None )
    return tree

#taxonomy overlap score
def getTaxOverlap(node):

    """
    Calculate the taxonomy overlap score for the given node in a phylogenetic tree.

    The taxonomy overlap score is defined as the number of taxonomic labels shared by all the leaf nodes
    descended from the given node, plus the sum of the scores of all its children. If a leaf node has no
    taxonomic label, it is not counted towards the score. The function also calculates the size of the
    largest loss in lineage length, defined as the difference between the length of the set of taxonomic
    labels shared by all the leaf nodes and the length of the longest set of taxonomic labels among the
    children of the node.

    The function adds the following features to the node object:
    - 'score': the taxonomy overlap score.
    - 'size': the largest loss in lineage length.
    - 'lineage': the set of taxonomic labels shared by all the leaf nodes descended from the node.

    Parameters:
    node (TreeNode): The node in a phylogenetic tree.

    Returns:
    set: The set of taxonomic labels shared by all the leaf nodes descended from the node, or `None` if
    the node has no children with taxonomic labels.
    """

    if node.is_leaf() == True:
        node.add_feature( 'score' ,  0 )
        node.add_feature( 'size' ,  0 )
        node.add_feature( 'score_x_frac' , 0)
        return node.lineage
    else:
        lengths = []
        total = 0
        redtotal = 0
        fractotal = 0
        sets = []
        scores = []
        for i,c in enumerate(node.get_children()):
            sets.append( getTaxOverlap(c))
            total += c.score
        sets = [s for s in sets if s]
        if len(sets)> 0:
            for i,cset in enumerate(sets):
                if i == 0:
                    nset = cset
                else:
                    nset = cset.intersection(nset)
                lengths.append(len(cset))
            #add the number of unique lineages
            score = len(nset) + total
            #add the number of unique lineages weighted by the fraction of the tree

            node.add_feature( 'score' ,  score )

            #show the biggest loss in lineage length
            node.add_feature( 'size' ,  abs( len(nset) - max(lengths) ) )

        else:
            nset = None
            node.add_feature( 'size' ,  0 )
            node.add_feature( 'score' ,  0 )

        node.add_feature( 'lineage' ,  nset )
        #only in the case of a leaf with no label
    return nset

def calc_score(t, lineage):
    try:
        tree = parse_newick(t)
    except ValueError as exc:
        eprint(f"Error parsing Newick tree: {exc}")
        sys.exit(1)

    uniprot_rows = read_lineage_tsv(lineage)

    if VERBOSE:
        eprint(uniprot_rows)

    lineages = make_lineages(uniprot_rows)
    species = get_species(uniprot_rows)
    tree = label_leaves(tree, lineages, species)
    getTaxOverlap(tree.treenode)
    return tree.treenode.score


if __name__ == "__main__":
    if "-h" in sys.argv or "--help" in sys.argv:
        print(__doc__)
        sys.exit(0)

    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    tree_file = sys.argv[1]
    try:
        tree = Path(tree_file).read_text(encoding="utf-8").strip()
    except FileNotFoundError:
        print(f"Error: The file '{tree_file}' does not exist.")
        sys.exit(1)
    except OSError as exc:
        print(f"Error reading tree file '{tree_file}': {exc}")
        sys.exit(1)

    if not tree:
        print("Error: The tree file is empty.")
        sys.exit(1)

    score = calc_score(tree, lineage=sys.argv[2])
    print(f"{tree_file}\t{score}")
