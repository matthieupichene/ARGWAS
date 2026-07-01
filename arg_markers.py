"""Incremental generation of ARG-derived genotype vectors.

This module identifies candidate ARG branches whose descendant sets may need to
be retested when moving from one local tree to the next. Each retained branch is
encoded as a diploid genotype vector, with values 0, 1, or 2 indicating how many
haplotypes of each individual descend from that branch.

The code assumes that diploid individuals are represented by two consecutive
sample nodes, so haplotype/sample node ``0`` and ``1`` belong to individual 0,
``2`` and ``3`` belong to individual 1, and so on.
"""

from collections.abc import Iterator

import numpy as np


MIN_DERIVED_ALLELE_FREQUENCY = 0.05


def _parent_array(tree, nodes: set[int], max_node: int) -> list[int]:
    """Return parent IDs for ``nodes`` in ``tree`` using -1 for missing/root nodes."""

    parents = [-1] * (max_node + 1)
    for node in nodes:
        parents[node] = tree.parent(node)
    return parents


def _candidate_nodes_changed_between_trees(tree1, tree2) -> set[int]:
    """Find nodes in ``tree2`` that should be re-encoded after a tree transition.

    Only part of the ARG changes between two adjacent local trees. This function
    compares parent relationships in the previous and current trees, then walks
    upward from affected parents to collect branches whose descendant sets may
    have changed.

    Parameters
    ----------
    tree1, tree2
        Consecutive ``tskit.Tree`` objects, where ``tree1`` is the previous tree
        and ``tree2`` is the current tree.

    Returns
    -------
    set[int]
        Node IDs in ``tree2`` to convert into genotype vectors.
    """

    nodes1 = set(tree1.nodes())
    nodes2 = set(tree2.nodes())
    all_nodes = nodes1 | nodes2

    max_node = max(max(nodes1, default=0), max(nodes2, default=0))
    parent1 = _parent_array(tree1, nodes1, max_node)
    parent2 = _parent_array(tree2, nodes2, max_node)

    # Parents involved in edges that differ between the two local trees.
    parents_from_tree1 = {
        parent1[node]
        for node in all_nodes
        if parent1[node] != -1 and parent1[node] != parent2[node]
    }
    parents_from_tree2 = {
        parent2[node]
        for node in all_nodes
        if parent2[node] != -1 and parent2[node] != parent1[node]
    }

    candidate_nodes: set[int] = set()

    # Add ancestors of parents introduced in tree2.
    for node in parents_from_tree2:
        while node != -1:
            candidate_nodes.add(node)
            node = parent2[node]

    # Add tree2 ancestors above the portion affected in tree1. This preserves
    # the original algorithm's behavior while making the traversal explicit.
    for node in parents_from_tree1:
        ancestor = node
        while ancestor != -1 and ancestor in parents_from_tree1:
            ancestor = parent1[ancestor]

        while ancestor != -1:
            candidate_nodes.add(ancestor)
            ancestor = parent2[ancestor]

    return candidate_nodes


def _node_to_diploid_genotype(tree, node: int, n_individuals: int) -> np.ndarray:
    """Encode descendants of ``node`` as a diploid 0/1/2 genotype vector."""

    leaves = np.fromiter(tree.leaves(node), dtype=np.int32)
    individuals = leaves // 2

    genotype = np.zeros(n_individuals, dtype=np.int8)
    np.add.at(genotype, individuals, 1)
    return genotype


def node_recalculate_genotypes(tree1, tree2, nIndiv: int) -> Iterator[np.ndarray]:
    """Yield ARG-derived genotype vectors affected by a local tree transition.

    Parameters
    ----------
    tree1, tree2
        Consecutive ``tskit.Tree`` objects. Candidate branches are encoded using
        the descendant sets in ``tree2``.
    nIndiv
        Number of diploid individuals. The current implementation assumes two
        haploid sample nodes per individual and converts sample node IDs to
        individual IDs using ``sample_node_id // 2``.

    Yields
    ------
    numpy.ndarray
        One int8 genotype vector of length ``nIndiv`` per candidate branch. Each
        value is 0, 1, or 2, corresponding to the number of descendant haplotypes
        carried by that individual.

    Notes
    -----
    Candidate branches with derived allele frequency below 5% or above 95% are
    skipped. This removes nearly fixed or nearly absent putative mutations before
    downstream association testing.
    """

    n_haplotypes = nIndiv * 2
    min_derived_count = int(n_haplotypes * MIN_DERIVED_ALLELE_FREQUENCY)
    max_derived_count = n_haplotypes - min_derived_count

    for node in _candidate_nodes_changed_between_trees(tree1, tree2):
        n_descendant_haplotypes = tree2.num_samples(node)

        if n_descendant_haplotypes < min_derived_count:
            continue
        if n_descendant_haplotypes > max_derived_count:
            continue

        yield _node_to_diploid_genotype(tree2, node, nIndiv)
