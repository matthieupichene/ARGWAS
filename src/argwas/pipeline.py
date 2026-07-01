"""Command-line entry point for ARGWAS.

ARGWAS scans ancestral recombination graph (ARG) tree sequences for putative
QTLs. For each local tree, it generates genotype-like vectors for candidate
mutations on ARG branches, whitens those vectors using a variance-component
model fitted from observed SNPs, and keeps the best-scoring ARG-derived marker
per local tree.

The public command-line interface lives in :mod:`argwas.cli`; this module
contains the reusable pipeline functions.
"""

import gc
import glob
import os
import re
from multiprocessing import Pool
from pathlib import Path
from typing import Iterable, NamedTuple

import h5py
import numpy as np
import pandas as pd
import tskit

from .hdf5_io import hdf5_block_reader, hdf5_close, hdf5_stream_writer, hdf5_write_snp
from .ancestry import infer_ancestral_alleles_from_args, write_ancestral_inference_report
from .arg_markers import node_recalculate_genotypes
from .whitening import compute_grm, ped_to_genotype_matrix, reml_vc


BLOCK_SIZE = 1_000
MIN_VARIANCE = 1e-8


class ChromosomeResult(NamedTuple):
    """Result files and association summaries for one chromosome/tree sequence."""

    chromosome: int
    marker_ids: list[str]
    positions: list[int]
    r_squared_values: list[float]
    betas: list[float]
    map_path: Path
    ped_path: Path


def natural_sort_key(text: str) -> list[int | str]:
    """Return a key that sorts strings containing numbers in human order.

    Example: ``chr2.trees`` is sorted before ``chr10.trees``.
    """

    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", text)]


def generate_arg_marker_hdf5(ts: tskit.TreeSequence, n_individuals: int, output_path: Path) -> None:
    """Generate ARG-derived genotype vectors and write them to HDF5.

    The first local tree receives a reference vector coded as homozygous derived
    for all individuals. For each following tree, only candidate nodes affected
    by the transition from the previous tree are recalculated.
    """

    # hdf5_stream_writer opens files in append mode, so remove stale output from
    # previous runs before writing a fresh chromosome-specific marker file.
    if output_path.exists():
        output_path.unlink()

    writer = hdf5_stream_writer(output_path)
    try:
        for tree_index in range(ts.num_trees):
            tree = ts.at_index(tree_index)

            if tree_index == 0:
                genotype = np.full(n_individuals, 2, dtype=np.int8)
                hdf5_write_snp(writer, genotype, tree_index)
                continue

            previous_tree = ts.at_index(tree_index - 1)
            for genotype in node_recalculate_genotypes(previous_tree, tree, n_individuals):
                hdf5_write_snp(writer, genotype, tree_index)
    finally:
        hdf5_close(writer)


def scan_arg_markers(
    marker_hdf5: Path,
    n_trees: int,
    u_file: Path,
    inv_sqrt_file: Path,
    phenotype_file: Path,
    block_size: int = BLOCK_SIZE,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Whiten ARG markers and keep the best R-squared value per tree.

    Returns
    -------
    best_r_squared
        Best coefficient of determination observed for each tree, computed
        in the whitened phenotype space.
    best_beta
        Regression effect estimate for the best-scoring marker in each tree.
    best_marker_index
        Row index in the HDF5 marker matrix for the best marker in each tree,
        or -1 if no valid marker was found.
    """

    U = np.load(u_file, mmap_mode="r")
    inv_sqrt_v = np.load(inv_sqrt_file, mmap_mode="r")
    phenotypes_w = np.load(phenotype_file)

    best_r_squared = np.zeros(n_trees)
    best_beta = np.zeros(n_trees)
    best_marker_index = np.full(n_trees, -1, dtype=int)

    phenotype_sum_squares = float(phenotypes_w @ phenotypes_w)
    if phenotype_sum_squares <= MIN_VARIANCE:
        raise ValueError("Whitened phenotype has near-zero variance; cannot compute R-squared.")

    global_start = 0

    for genotype_block, tree_ids in hdf5_block_reader(marker_hdf5, block_size=block_size):
        # hdf5_block_reader returns markers as columns: n_individuals x n_markers.
        genotype_block = genotype_block.astype(np.float32)
        genotype_block -= genotype_block.mean(axis=0)

        rotated = U.T @ genotype_block
        whitened = inv_sqrt_v[:, None] * rotated

        numerator = phenotypes_w @ whitened
        variance = np.sum(whitened * whitened, axis=0)
        valid = variance > MIN_VARIANCE

        if not np.any(valid):
            global_start += genotype_block.shape[1]
            continue

        # Keep the original HDF5 row index after filtering. This avoids selecting
        # the wrong SNP vector when some columns are removed by the variance mask.
        marker_indices = global_start + np.arange(genotype_block.shape[1])

        numerator = numerator[valid]
        variance = variance[valid]
        tree_ids = tree_ids[valid]
        marker_indices = marker_indices[valid]

        beta = numerator / variance
        # The raw score is the regression sum of squares: (x'y)^2 / (x'x).
        # Dividing by y'y gives the usual coefficient of determination for a
        # one-marker linear regression in the whitened space. Do not divide by
        # phenotype variance squared; R^2 is a fraction of sums of squares.
        r_squared = (numerator * numerator / variance) / phenotype_sum_squares

        for tree_id, marker_index, marker_beta, marker_r_squared in zip(tree_ids, marker_indices, beta, r_squared):
            if marker_r_squared > best_r_squared[tree_id]:
                best_r_squared[tree_id] = marker_r_squared
                best_beta[tree_id] = marker_beta
                best_marker_index[tree_id] = marker_index

        global_start += genotype_block.shape[1]

    return best_r_squared, best_beta, best_marker_index


def write_chromosome_outputs(
    ts: tskit.TreeSequence,
    chromosome: int,
    marker_hdf5: Path,
    best_r_squared: np.ndarray,
    best_beta: np.ndarray,
    best_marker_index: np.ndarray,
    map_path: Path,
    ped_path: Path,
) -> tuple[list[str], list[int], list[float], list[float]]:
    """Write temporary MAP/PED-like files for the best ARG marker per tree."""

    marker_ids: list[str] = []
    positions: list[int] = []
    r_squared_values: list[float] = []
    betas: list[float] = []

    with open(map_path, "w") as map_file, open(ped_path, "w") as ped_file, h5py.File(marker_hdf5, "r") as h5_file:
        marker_dataset = h5_file["SNPs"]

        for tree_index in range(ts.num_trees):
            marker_index = best_marker_index[tree_index]
            if marker_index == -1:
                continue

            tree = ts.at_index(tree_index)
            position = int(tree.interval.mid)
            marker_name = f"ARG{tree_index}_{chromosome}"

            marker_ids.append(marker_name)
            positions.append(position)
            r_squared_values.append(float(best_r_squared[tree_index]))
            betas.append(float(best_beta[tree_index]))

            map_file.write(f"{chromosome} {marker_name} 0 {position}\n")

            marker_vector = marker_dataset[marker_index, :]
            ped_file.write(f"{marker_name} {' '.join(map(str, marker_vector))}\n")

    return marker_ids, positions, r_squared_values, betas


def process_chromosome(args: tuple[str, int, int, str, str, str, str, str, str]) -> ChromosomeResult:
    """Process one chromosome/tree sequence file.

    This function is top-level so it can be used by ``multiprocessing.Pool``.
    """

    (
        tree_sequence_path,
        chromosome,
        n_individuals,
        u_file,
        inv_sqrt_file,
        phenotype_file,
        map_path,
        ped_path,
        marker_hdf5,
    ) = args

    ts = tskit.load(tree_sequence_path)
    marker_hdf5 = Path(marker_hdf5)

    generate_arg_marker_hdf5(ts, n_individuals, marker_hdf5)
    best_r_squared, best_beta, best_marker_index = scan_arg_markers(
        marker_hdf5=marker_hdf5,
        n_trees=ts.num_trees,
        u_file=Path(u_file),
        inv_sqrt_file=Path(inv_sqrt_file),
        phenotype_file=Path(phenotype_file),
    )
    marker_ids, positions, r_squared_values, betas = write_chromosome_outputs(
        ts=ts,
        chromosome=chromosome,
        marker_hdf5=marker_hdf5,
        best_r_squared=best_r_squared,
        best_beta=best_beta,
        best_marker_index=best_marker_index,
        map_path=Path(map_path),
        ped_path=Path(ped_path),
    )

    return ChromosomeResult(chromosome, marker_ids, positions, r_squared_values, betas, Path(map_path), Path(ped_path))


def fit_whitening_from_ped(
    ped_table: pd.DataFrame,
    map_table: pd.DataFrame,
    ancestral_alleles,
    output_dir: Path,
    ancestral_coordinate_system: str = "one-based",
    ancestral_is_per_snp: bool = False,
) -> tuple[Path, Path, Path]:
    """Fit the GRM/REML whitening transform from observed SNP genotypes."""

    # ``to_numpy`` can return a read-only view depending on the pandas
    # backend/input source. Make an explicit writable copy before centering.
    phenotypes = ped_table[5].to_numpy(dtype=float, copy=True)
    phenotypes -= np.mean(phenotypes)

    print("Whitening observed SNP genotypes...")
    genotype_matrix = ped_to_genotype_matrix(
        ped_df=ped_table,
        map_df=map_table,
        ancestral_alleles=ancestral_alleles,
        coordinate_system=ancestral_coordinate_system,
        ancestral_is_per_snp=ancestral_is_per_snp,
    )
    grm = compute_grm(genotype_matrix)
    (sigma_g2, sigma_e2), eigvals, U = reml_vc(phenotypes, grm)

    inv_sqrt_v = 1.0 / np.sqrt(sigma_g2 * eigvals + sigma_e2 + 1e-10)
    phenotypes_w = inv_sqrt_v * (U.T @ phenotypes)

    u_file = output_dir / "U.npy"
    inv_sqrt_file = output_dir / "inv_sqrt_v.npy"
    phenotype_file = output_dir / "phenotypes_w.npy"

    np.save(u_file, U)
    np.save(inv_sqrt_file, inv_sqrt_v)
    np.save(phenotype_file, phenotypes_w)

    del genotype_matrix, grm, U, inv_sqrt_v
    gc.collect()

    return u_file, inv_sqrt_file, phenotype_file


def merge_temporary_outputs(results: Iterable[ChromosomeResult], output_dir: Path) -> tuple[Path, Path]:
    """Merge chromosome-level temporary MAP/PED files into project outputs."""

    map_path = output_dir / "results.map"
    ped_path = output_dir / "results.ped"

    with open(map_path, "w") as output_map, open(ped_path, "w") as output_ped:
        for result in results:
            output_map.write(result.map_path.read_text())
            output_ped.write(result.ped_path.read_text())

    return map_path, ped_path


def write_final_plink_files(ped_table: pd.DataFrame, map_path: Path, ped_path: Path) -> None:
    """Convert temporary ARG marker files into PLINK-style MAP/PED files.

    ARG-derived markers are represented as biallelic putative mutations with
    ancestral allele ``A`` and derived allele ``T``.
    """

    map_df = pd.read_csv(map_path, sep=r"\s+", header=None)
    ped_df = pd.read_csv(ped_path, sep=r"\s+", header=None)

    # Sort ARG-derived markers by chromosome and position, keeping the matching
    # rows in the temporary marker-dosage table in the same order.
    sort_order = map_df.sort_values([0, 3]).index
    map_df = map_df.loc[sort_order].reset_index(drop=True)
    ped_df = ped_df.loc[sort_order].reset_index(drop=True)

    output_ped = ped_table.iloc[:, :6].copy()
    genotypes = ped_df.iloc[:, 1:].to_numpy(copy=True).T

    n_individuals, n_markers = genotypes.shape
    alleles = np.empty((n_individuals, n_markers * 2), dtype=object)
    alleles[:, 0::2] = np.where(genotypes >= 1, "T", "A")
    alleles[:, 1::2] = np.where(genotypes == 2, "T", "A")

    genotype_df = pd.DataFrame(alleles)
    output_ped = pd.concat([output_ped.reset_index(drop=True), genotype_df], axis=1)

    map_df.to_csv(map_path, sep="\t", header=False, index=False)
    output_ped.to_csv(ped_path, sep="\t", header=False, index=False)



def run_argwas(
    output: str | Path,
    ped: str | Path,
    trees: str | Path,
    ances: str | Path | None = None,
    workers: int = 4,
    map_file: str | Path | None = None,
    ancestral_coordinate_system: str = "one-based",
    ancestral_mode: str = "file",
    ambiguous_ancestral: str = "skip",
) -> None:
    """Run the complete ARGWAS workflow.

    Parameters
    ----------
    output
        Directory where output files will be written.
    ped
        Input PLINK PED file containing metadata, phenotype, and observed SNPs.
    trees
        Directory containing chromosome-level ``.trees`` files.
    ances
        File containing the genome-wide ancestral allele sequence. Required when
        ``ancestral_mode="file"`` and ignored when ``ancestral_mode="arg"``.
    workers
        Number of worker processes used to scan tree-sequence files.
    map_file
        PLINK MAP file matching the input PED file. If omitted, ARGWAS looks for
        a file with the same prefix as the PED and suffix ``.map``.
    ancestral_coordinate_system
        Coordinate convention for MAP physical positions. ``"one-based"`` is the
        default and corresponds to standard VCF/PLINK coordinates, where BP=1
        indexes the first ancestral allele / tree-sequence coordinate 0. Use
        ``"zero-based"`` only for MAP files generated with zero-based positions.
    ancestral_mode
        ``"file"`` uses the provided ancestral sequence. ``"arg"`` infers the
        ancestral allele for each observed SNP from the local ARG topology.
    ambiguous_ancestral
        Policy for ARG-inferred ambiguous sites: ``"skip"`` omits the observed
        SNP from the whitening GRM, ``"major"`` uses the major observed allele,
        and ``"reference"`` uses the first observed PED allele. This affects
        only the real SNPs used for whitening, not the ARG-derived markers tested
        by ARGWAS.
    """
    output_dir = Path(output)
    ped_path = Path(ped)
    map_path = Path(map_file) if map_file is not None else ped_path.with_suffix(".map")
    trees_dir = Path(trees)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not map_path.exists():
        raise FileNotFoundError(
            f"Could not find MAP file {map_path}. Provide it explicitly with --map."
        )

    ped_table = pd.read_csv(ped_path, header=None, sep=r"\s+")
    map_table = pd.read_csv(map_path, header=None, sep=r"\s+")

    tree_files = glob.glob(str(trees_dir / "*.trees"))
    tree_files.sort(key=natural_sort_key)
    if not tree_files:
        raise FileNotFoundError(f"No .trees files found in {trees_dir}")

    ancestral_is_per_snp = False
    if ancestral_mode == "file":
        if ances is None:
            raise ValueError("--ances is required when --ancestral-mode file is used.")
        ancestral_alleles = list(Path(ances).read_text().strip())
        ancestral_is_per_snp = False
    elif ancestral_mode == "arg":
        print("Inferring ancestral alleles from ARG topology...")
        inference = infer_ancestral_alleles_from_args(
            ped_df=ped_table,
            map_df=map_table,
            trees_dir=trees_dir,
            tree_files=tree_files,
            coordinate_system=ancestral_coordinate_system,
            ambiguous_policy=ambiguous_ancestral,
        )
        ancestral_alleles = inference.alleles
        ancestral_is_per_snp = True
        report_path = write_ancestral_inference_report(
            output_dir / "ancestral_inference.csv", map_table, inference
        )
        print("Ancestral inference summary:")
        for key in sorted(inference.summary):
            print(f"  {key}: {inference.summary[key]}")
        print(f"Wrote ancestral inference report: {report_path}")
    else:
        raise ValueError("ancestral_mode must be 'file' or 'arg'.")

    n_individuals = len(ped_table)
    u_file, inv_sqrt_file, phenotype_file = fit_whitening_from_ped(
        ped_table=ped_table,
        map_table=map_table,
        ancestral_alleles=ancestral_alleles,
        output_dir=output_dir,
        ancestral_coordinate_system=ancestral_coordinate_system,
        ancestral_is_per_snp=ancestral_is_per_snp,
    )

    pool_args = []
    for chromosome_index, tree_file in enumerate(tree_files, start=1):
        pool_args.append(
            (
                tree_file,
                chromosome_index,
                n_individuals,
                str(u_file),
                str(inv_sqrt_file),
                str(phenotype_file),
                str(output_dir / f"temp_chr{chromosome_index}.map"),
                str(output_dir / f"temp_chr{chromosome_index}.ped"),
                str(output_dir / f"temp_snp{chromosome_index}.h5"),
            )
        )

    if workers == 1:
        chromosome_results = [process_chromosome(args) for args in pool_args]
    else:
        with Pool(processes=workers) as pool:
            chromosome_results = pool.map(process_chromosome, pool_args)

    summary_rows = []
    for result in chromosome_results:
        for marker_id, position, r_squared, beta in zip(
            result.marker_ids, result.positions, result.r_squared_values, result.betas
        ):
            summary_rows.append((result.chromosome, marker_id, position, r_squared, beta))

    summary = pd.DataFrame(
        summary_rows, columns=["chromosome", "marker_id", "position", "r_squared", "beta"]
    )
    summary.to_csv(output_dir / "results.csv", index=False)

    map_path, ped_path = merge_temporary_outputs(chromosome_results, output_dir)
    write_final_plink_files(ped_table, map_path, ped_path)
