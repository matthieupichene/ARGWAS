"""Utilities for storing and reading ARG-derived marker genotypes in HDF5.

The ARGWAS pipeline can generate many putative markers from local trees in an
ancestral recombination graph. Keeping all markers in memory can be expensive,
so this module provides a simple streaming writer and block reader.

Storage layout
--------------
The HDF5 file contains two datasets by default:

``SNPs``
    Two-dimensional integer dataset with shape ``(n_markers, n_individuals)``.
    Each row is one ARG-derived marker / putative mutation and each column is an
    individual genotype dosage.

``tree_id``
    One-dimensional integer dataset with shape ``(n_markers,)``. Entry ``i``
    gives the local tree index from which marker ``i`` was generated.
"""

from pathlib import Path

import h5py
import numpy as np


DEFAULT_SNP_DATASET = "SNPs"
DEFAULT_TREE_DATASET = "tree_id"
DEFAULT_BLOCK_SIZE = 1_000
HDF5_COMPRESSION = "gzip"

def hdf5_block_reader(
    hdf5_file: str | Path,
    snp_dataset: str = DEFAULT_SNP_DATASET,
    tree_dataset: str = DEFAULT_TREE_DATASET,
    block_size: int = DEFAULT_BLOCK_SIZE,
):
    """Yield genotype markers from an HDF5 file in column-oriented blocks.

    Parameters
    ----------
    hdf5_file
        Path to the HDF5 file produced by :func:`hdf5_stream_writer` and
        :func:`hdf5_write_snp`.
    snp_dataset
        Name of the genotype dataset. The expected shape is
        ``(n_markers, n_individuals)``.
    tree_dataset
        Name of the dataset containing the tree index for each marker. The
        expected shape is ``(n_markers,)``.
    block_size
        Number of markers to read per block.

    Yields
    ------
    tuple[np.ndarray, np.ndarray]
        ``(G_block, tree_block)`` where ``G_block`` has shape
        ``(n_individuals, n_markers_in_block)`` and ``tree_block`` has shape
        ``(n_markers_in_block,)``.

    Notes
    -----
    The transpose is intentional: markers are stored as rows on disk, but many
    matrix operations in the pipeline expect individuals as rows and markers as
    columns.
    """
    if block_size <= 0:
        raise ValueError("block_size must be a positive integer.")

    with h5py.File(hdf5_file, "r") as hdf5:
        if snp_dataset not in hdf5:
            raise KeyError(f"Dataset {snp_dataset!r} was not found in {hdf5_file!s}.")
        if tree_dataset not in hdf5:
            raise KeyError(f"Dataset {tree_dataset!r} was not found in {hdf5_file!s}.")

        genotypes = hdf5[snp_dataset]
        tree_ids = hdf5[tree_dataset]

        if genotypes.ndim != 2:
            raise ValueError(
                f"Dataset {snp_dataset!r} must be two-dimensional; "
                f"found shape {genotypes.shape}."
            )
        if tree_ids.ndim != 1:
            raise ValueError(
                f"Dataset {tree_dataset!r} must be one-dimensional; "
                f"found shape {tree_ids.shape}."
            )
        if genotypes.shape[0] != tree_ids.shape[0]:
            raise ValueError(
                "The genotype and tree-id datasets contain different numbers "
                f"of markers: {genotypes.shape[0]} vs {tree_ids.shape[0]}."
            )

        n_markers = genotypes.shape[0]
        for start in range(0, n_markers, block_size):
            end = min(start + block_size, n_markers)

            # Stored on disk as (marker, individual), returned as
            # (individual, marker) for downstream matrix operations.
            genotype_block = genotypes[start:end, :].T
            tree_block = tree_ids[start:end]

            yield genotype_block, tree_block


def hdf5_stream_writer(
    hdf5_file: str | Path,
    dataset: str = DEFAULT_SNP_DATASET,
    tree_dataset: str = DEFAULT_TREE_DATASET,
):
    """Open an HDF5 file for streaming marker-by-marker genotype writing.

    Parameters
    ----------
    hdf5_file
        Path to the HDF5 file to create or append to.
    dataset
        Name of the genotype dataset to create/use.
    tree_dataset
        Name of the tree-id dataset to create/use.

    Returns
    -------
    dict
        A lightweight writer object used by :func:`hdf5_write_snp` and
        :func:`hdf5_close`.

    Notes
    -----
    The file is opened in append mode. If the datasets already exist, new
    markers are appended. Delete the output file before calling this function if
    you want a fresh run.
    """
    hdf5 = h5py.File(hdf5_file, "a")

    if dataset in hdf5:
        genotype_dataset = hdf5[dataset]
        if genotype_dataset.ndim != 2:
            hdf5.close()
            raise ValueError(
                f"Existing dataset {dataset!r} must be two-dimensional; "
                f"found shape {genotype_dataset.shape}."
            )
        initialized = True
    else:
        genotype_dataset = hdf5.create_dataset(
            dataset,
            shape=(0, 0),
            maxshape=(None, None),
            dtype=np.int8,
            chunks=True,
            compression=HDF5_COMPRESSION,
        )
        initialized = False

    if tree_dataset in hdf5:
        tree_id_dataset = hdf5[tree_dataset]
        if tree_id_dataset.ndim != 1:
            hdf5.close()
            raise ValueError(
                f"Existing dataset {tree_dataset!r} must be one-dimensional; "
                f"found shape {tree_id_dataset.shape}."
            )
    else:
        tree_id_dataset = hdf5.create_dataset(
            tree_dataset,
            shape=(0,),
            maxshape=(None,),
            dtype=np.int32,
            chunks=True,
            compression=HDF5_COMPRESSION,
        )

    if genotype_dataset.shape[0] != tree_id_dataset.shape[0]:
        hdf5.close()
        raise ValueError(
            "Existing genotype and tree-id datasets have inconsistent lengths: "
            f"{genotype_dataset.shape[0]} vs {tree_id_dataset.shape[0]}."
        )

    return {
        "file": hdf5,
        "dset": genotype_dataset,
        "tree": tree_id_dataset,
        "initialized": initialized,
    }


def hdf5_write_snp(writer, genotype, tree_id):
    """Append one marker genotype vector to an open HDF5 writer.

    Parameters
    ----------
    writer
        Writer returned by :func:`hdf5_stream_writer`.
    genotype
        One-dimensional genotype dosage vector with one value per individual.
    tree_id
        Index of the ARG local tree from which this marker was generated.
    """
    genotype = np.asarray(genotype, dtype=np.int8)
    if genotype.ndim != 1:
        raise ValueError(
            "genotype must be a one-dimensional vector with one value per "
            f"individual; found shape {genotype.shape}."
        )

    genotype_dataset = writer["dset"]
    tree_id_dataset = writer["tree"]
    n_individuals = genotype.shape[0]

    if not writer["initialized"]:
        genotype_dataset.resize((0, n_individuals))
        writer["initialized"] = True
    elif genotype_dataset.shape[1] != n_individuals:
        raise ValueError(
            "All genotype vectors written to the same HDF5 dataset must have "
            f"the same length. Existing length: {genotype_dataset.shape[1]}; "
            f"new length: {n_individuals}."
        )

    old_n_markers = genotype_dataset.shape[0]
    new_n_markers = old_n_markers + 1

    genotype_dataset.resize((new_n_markers, n_individuals))
    tree_id_dataset.resize((new_n_markers,))

    genotype_dataset[old_n_markers, :] = genotype
    tree_id_dataset[old_n_markers] = int(tree_id)


def hdf5_close(writer):
    """Close an HDF5 writer returned by :func:`hdf5_stream_writer`."""
    writer["file"].close()
