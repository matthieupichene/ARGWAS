"""Ancestral-allele inference utilities for ARGWAS.

ARGWAS can either use a user-provided ancestral sequence or infer ancestral
states from an ARG. ARG-based inference is intentionally conservative: an allele
is treated as derived only when its carrier haplotypes exactly form a clade in
that local tree. Ambiguous cases are reported rather than silently guessed.
"""

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd
import tskit


VALID_AMBIGUOUS_POLICIES = {"skip", "major", "reference"}
VALID_COORDINATE_SYSTEMS = {"one-based", "zero-based"}


@dataclass(frozen=True)
class InferredAncestralAlleles:
    """Result of ARG-based ancestral allele inference."""

    alleles: list[str | None]
    status: list[str]
    summary: Counter


def map_positions_to_indices(
    map_df: pd.DataFrame,
    sequence_length: int,
    coordinate_system: str = "one-based",
) -> np.ndarray:
    """Convert PLINK MAP BP coordinates to zero-based sequence indices.

    Parameters
    ----------
    map_df
        PLINK MAP table with physical positions in column 4.
    sequence_length
        Length of the target sequence/array being indexed.
    coordinate_system
        ``"one-based"`` means BP=1 maps to index 0. ``"zero-based"`` means BP=0
        maps to index 0.
    """
    if coordinate_system not in VALID_COORDINATE_SYSTEMS:
        raise ValueError("coordinate_system must be 'one-based' or 'zero-based'.")
    if map_df.shape[1] < 4:
        raise ValueError("MAP file must contain at least four columns: CHR SNP CM BP.")

    positions = map_df.iloc[:, 3].to_numpy(dtype=np.int64, copy=True)
    if coordinate_system == "one-based":
        if np.any(positions < 1):
            bad = int(positions[positions < 1][0])
            raise ValueError(
                "MAP positions contain a value < 1 while using one-based coordinates "
                f"(first offending position: {bad}). Use --ancestral-coordinate-system "
                "zero-based if your MAP positions are zero-based."
            )
        indices = positions - 1
    else:
        if np.any(positions < 0):
            bad = int(positions[positions < 0][0])
            raise ValueError(f"MAP positions contain a negative coordinate ({bad}).")
        indices = positions

    if np.any(indices >= sequence_length):
        bad_position = int(positions[indices >= sequence_length][0])
        bad_index = int(indices[indices >= sequence_length][0])
        raise ValueError(
            "MAP position is outside the indexed sequence: "
            f"position {bad_position} maps to Python index {bad_index}, but length is "
            f"{sequence_length}. Check chromosome/coordinate consistency."
        )
    return indices


def _allele_carriers_from_ped_row(
    ped_df: pd.DataFrame,
    snp_index: int,
) -> dict[str, frozenset[int]]:
    """Return haploid sample-node carriers for each allele at one PED SNP.

    This preserves the original ARGWAS assumption that diploid individuals are
    encoded as consecutive haploid sample nodes: individual i has sample nodes
    2*i and 2*i+1.
    """
    genotype_columns = ped_df.iloc[:, 6:]
    allele_1 = genotype_columns.iloc[:, 2 * snp_index].to_numpy()
    allele_2 = genotype_columns.iloc[:, 2 * snp_index + 1].to_numpy()

    carriers: dict[str, set[int]] = {}
    for individual_index, (a1, a2) in enumerate(zip(allele_1, allele_2)):
        if a1 != "0":
            carriers.setdefault(str(a1), set()).add(2 * individual_index)
        if a2 != "0":
            carriers.setdefault(str(a2), set()).add(2 * individual_index + 1)
    return {allele: frozenset(nodes) for allele, nodes in carriers.items()}


def _is_exact_clade(tree: tskit.Tree, carrier_nodes: frozenset[int], total_samples: int) -> bool:
    """Return True if carrier sample nodes exactly match a non-root clade."""
    n_carriers = len(carrier_nodes)
    if n_carriers == 0 or n_carriers == total_samples:
        return False

    # The MRCA of the carriers is the only possible node whose descendant sample
    # set could equal the carrier set. This avoids scanning every node.
    carrier_list = list(carrier_nodes)
    mrca = carrier_list[0]
    for node in carrier_list[1:]:
        mrca = tree.mrca(mrca, node)
        if mrca == -1:
            return False

    descendant_samples = frozenset(tree.samples(mrca))
    return descendant_samples == carrier_nodes


def infer_ancestral_allele_for_snp(
    tree: tskit.Tree,
    allele_carriers: dict[str, frozenset[int]],
    total_samples: int,
) -> tuple[str | None, str]:
    """Infer one SNP ancestral allele from local tree topology.

    Returns ``(allele, status)``. ``status`` is ``"inferred"`` when exactly one
    allele forms a clean derived clade, ``"ambiguous"`` when orientation is not
    unique, and ``"missing"`` when the SNP is not biallelic or has missing calls.
    """
    if len(allele_carriers) != 2:
        return None, "missing"

    alleles = list(allele_carriers)
    derived_candidates = [
        allele
        for allele in alleles
        if _is_exact_clade(tree, allele_carriers[allele], total_samples=total_samples)
    ]

    if len(derived_candidates) != 1:
        return None, "ambiguous"

    derived = derived_candidates[0]
    ancestral = next(allele for allele in alleles if allele != derived)
    return ancestral, "inferred"


def _fallback_allele(
    allele_carriers: dict[str, frozenset[int]],
    policy: str,
) -> tuple[str | None, str]:
    """Choose a fallback ancestral allele for ambiguous ARG inference."""
    if policy == "skip":
        return None, "skipped"
    if len(allele_carriers) == 0:
        return None, "missing"

    # In the absence of a reliable ancestral state, the major allele is the least
    # surprising approximation for derived-dosage coding. ``reference`` uses the
    # first allele in PED order, which is sometimes useful for debugging only.
    if policy == "major":
        return max(allele_carriers, key=lambda allele: len(allele_carriers[allele])), "fallback_major"
    if policy == "reference":
        return next(iter(allele_carriers)), "fallback_reference"
    raise ValueError(f"Unknown ambiguous ancestral policy: {policy}")


def infer_ancestral_alleles_from_args(
    ped_df: pd.DataFrame,
    map_df: pd.DataFrame,
    trees_dir,
    tree_files,
    coordinate_system: str = "one-based",
    ambiguous_policy: str = "skip",
) -> InferredAncestralAlleles:
    """Infer ancestral alleles for observed PED SNPs using chromosome ARGs.

    MAP column 1 is interpreted as the chromosome number matching the natural
    sort order of the supplied ``.trees`` files: chromosome 1 -> first tree file,
    chromosome 2 -> second tree file, and so on. MAP BP coordinates are converted
    to zero-based tree-sequence coordinates using ``coordinate_system``.
    """
    if ambiguous_policy not in VALID_AMBIGUOUS_POLICIES:
        raise ValueError(
            "ambiguous_policy must be one of: " + ", ".join(sorted(VALID_AMBIGUOUS_POLICIES))
        )

    genotype_columns = ped_df.iloc[:, 6:]
    n_snps = genotype_columns.shape[1] // 2
    if genotype_columns.shape[1] % 2 != 0:
        raise ValueError("PED genotype columns should contain exactly two allele columns per SNP.")
    if len(map_df) != n_snps:
        raise ValueError(f"MAP has {len(map_df)} rows but PED contains {n_snps} SNPs.")

    sorted_tree_files = [Path(p) for p in tree_files]
    chromosome_to_tree = {chromosome: path for chromosome, path in enumerate(sorted_tree_files, start=1)}

    alleles: list[str | None] = [None] * n_snps
    status: list[str] = ["missing"] * n_snps
    summary: Counter = Counter()
    loaded_trees: dict[int, tskit.TreeSequence] = {}

    for chromosome_value, row_indices in map_df.groupby(map_df.iloc[:, 0]).groups.items():
        try:
            chromosome = int(chromosome_value)
        except ValueError as exc:
            raise ValueError(
                f"ARG-based ancestral inference currently expects numeric MAP chromosomes; "
                f"got {chromosome_value!r}."
            ) from exc

        tree_path = chromosome_to_tree.get(chromosome)
        if tree_path is None:
            for snp_index in row_indices:
                allele_carriers = _allele_carriers_from_ped_row(ped_df, int(snp_index))
                allele, fallback_status = _fallback_allele(allele_carriers, ambiguous_policy)
                alleles[int(snp_index)] = allele
                status[int(snp_index)] = "missing_tree" if allele is None else fallback_status
                summary[status[int(snp_index)]] += 1
            continue

        ts = loaded_trees.get(chromosome)
        if ts is None:
            ts = tskit.load(str(tree_path))
            loaded_trees[chromosome] = ts

        sub_map = map_df.loc[list(row_indices)]
        tree_indices = map_positions_to_indices(
            sub_map,
            sequence_length=int(ts.sequence_length),
            coordinate_system=coordinate_system,
        )

        for snp_index, tree_index in zip(row_indices, tree_indices):
            snp_index = int(snp_index)
            allele_carriers = _allele_carriers_from_ped_row(ped_df, snp_index)
            try:
                tree = ts.at(float(tree_index))
                allele, call_status = infer_ancestral_allele_for_snp(
                    tree,
                    allele_carriers,
                    total_samples=2 * ped_df.shape[0],
                )
            except ValueError:
                # Usually means the position is outside the tree sequence.
                # Other exceptions should surface, because they may indicate a
                # real bug rather than an uninferable ancestral state.
                allele, call_status = None, "missing"

            if allele is None:
                allele, fallback_status = _fallback_allele(allele_carriers, ambiguous_policy)
                call_status = call_status if allele is None else fallback_status

            alleles[snp_index] = allele
            status[snp_index] = call_status
            summary[call_status] += 1

    return InferredAncestralAlleles(alleles=alleles, status=status, summary=summary)


def write_ancestral_inference_report(
    output_path: str | Path,
    map_df: pd.DataFrame,
    result: InferredAncestralAlleles,
) -> Path:
    """Write a per-SNP report for ARG-inferred ancestral alleles."""
    output_path = Path(output_path)
    report = pd.DataFrame(
        {
            "chromosome": map_df.iloc[:, 0].to_numpy(),
            "marker_id": map_df.iloc[:, 1].to_numpy(),
            "position": map_df.iloc[:, 3].to_numpy(),
            "ancestral_allele": result.alleles,
            "ancestral_status": result.status,
            "used_for_whitening": [allele is not None for allele in result.alleles],
        }
    )
    report.to_csv(output_path, index=False)
    return output_path
