"""Utilities for genotype whitening and variance-component estimation.

This module contains helper functions used by ARGWAS to account for relatedness
between individuals before testing ARG-derived markers. The workflow is:

1. Convert a PLINK PED table into a genotype dosage matrix.
2. Build a genomic relationship matrix (GRM) from observed SNPs.
3. Estimate genetic and residual variance components with a simple REML model.
4. Reuse the eigendecomposition of the GRM for phenotype/marker whitening.

The functions here are intentionally small and independent from the ARG-specific
code, so they can be tested separately.
"""

import numpy as np
from scipy.optimize import minimize


EPSILON = 1e-8
MIN_VARIANCE = 1e-12
THETA_BOUNDS = (-20.0, 20.0)


from .ancestry import map_positions_to_indices

def ped_to_genotype_matrix(
    ped_df,
    map_df,
    ancestral_alleles,
    coordinate_system: str = "one-based",
    ancestral_is_per_snp: bool = False,
) -> np.ndarray:
    """Convert PLINK PED alleles into derived-allele dosages.

    Ancestral alleles are used only to polarize the *observed SNPs* used for
    GRM construction and spectral whitening. ARG-derived putative markers are
    generated and tested later, and do not depend on this function.

    Parameters
    ----------
    ped_df
        PLINK PED table. The first six columns are metadata; remaining columns
        are allele calls, two columns per SNP.
    map_df
        Matching PLINK MAP table. Row order must match the SNP order in the PED.
    ancestral_alleles
        Either a genome-wide ancestral sequence indexed by MAP BP positions
        (file mode), or a per-SNP list returned by ARG-based inference
        (ARG mode). Unknown per-SNP entries should be ``None``.
    coordinate_system
        Coordinate convention for MAP BP positions when ``ancestral_is_per_snp``
        is False. With the default ``"one-based"``, MAP BP=1 maps to Python
        index 0 in the ancestral sequence.
    ancestral_is_per_snp
        Set to True when ``ancestral_alleles`` already has one entry per SNP in
        PED/MAP order, as with ARG-based ancestral inference. Set to False for
        a genome-wide ancestral sequence from a file.

    Returns
    -------
    numpy.ndarray
        Genotype dosage matrix with one column per usable observed SNP. Sites
        whose ancestral allele is unknown are skipped entirely, so they do not
        enter the GRM as all-zero columns.
    """
    genotype_columns = ped_df.iloc[:, 6:]
    n_individuals = ped_df.shape[0]
    n_snps = genotype_columns.shape[1] // 2

    if genotype_columns.shape[1] % 2 != 0:
        raise ValueError("PED genotype columns should contain exactly two allele columns per SNP.")

    if len(map_df) != n_snps:
        raise ValueError(
            "The MAP file must contain exactly one row per SNP in the PED file: "
            f"{len(map_df)} MAP rows for {n_snps} PED SNPs."
        )

    if ancestral_is_per_snp:
        if len(ancestral_alleles) != n_snps:
            raise ValueError(
                "Per-SNP ancestral alleles must contain exactly one entry per PED SNP: "
                f"{len(ancestral_alleles)} ancestral entries for {n_snps} SNPs."
            )
        ancestral_lookup = list(ancestral_alleles)
    else:
        ancestral_indices = map_positions_to_indices(
            map_df=map_df,
            sequence_length=len(ancestral_alleles),
            coordinate_system=coordinate_system,
        )
        ancestral_lookup = [ancestral_alleles[int(index)] for index in ancestral_indices]

    genotype_columns_to_keep = []

    for snp_index, ancestral_allele in enumerate(ancestral_lookup):
        if ancestral_allele is None:
            # Unknown ancestral state: omit this observed SNP from whitening.
            # This does not remove any ARG-derived marker from the association
            # scan; it only reduces the real-SNP panel used to build the GRM.
            continue

        allele_1 = genotype_columns.iloc[:, 2 * snp_index].to_numpy()
        allele_2 = genotype_columns.iloc[:, 2 * snp_index + 1].to_numpy()
        dosage = (
            (allele_1 != ancestral_allele).astype(np.int8)
            + (allele_2 != ancestral_allele).astype(np.int8)
        )
        genotype_columns_to_keep.append(dosage)

    if not genotype_columns_to_keep:
        raise ValueError(
            "No observed SNPs have a usable ancestral allele for whitening. "
            "Provide an ancestral file, use a different --ambiguous-ancestral policy, "
            "or check ARG ancestral inference."
        )

    return np.column_stack(genotype_columns_to_keep).astype(np.int8, copy=False)

def compute_grm(genotype_matrix: np.ndarray) -> np.ndarray:
    """Compute a genomic relationship matrix from genotype dosages.

    Parameters
    ----------
    genotype_matrix
        Array of shape ``(n_individuals, n_snps)`` containing diploid genotype
        dosages encoded as 0, 1, or 2.

    Returns
    -------
    numpy.ndarray
        Genomic relationship matrix of shape ``(n_individuals, n_individuals)``.

    Notes
    -----
    SNPs are standardized using allele frequency ``p`` and variance
    ``2p(1-p)``. A small epsilon is added to avoid division by zero for nearly
    fixed markers.
    """
    genotypes = genotype_matrix.astype(np.float32, copy=False)

    if genotypes.ndim != 2:
        raise ValueError("genotype_matrix must be a 2D array.")
    if genotypes.shape[1] == 0:
        raise ValueError("Cannot compute a GRM from zero SNPs.")

    allele_frequencies = genotypes.mean(axis=0) / 2.0
    expected_variance = 2.0 * allele_frequencies * (1.0 - allele_frequencies)
    standardized_genotypes = (genotypes - 2.0 * allele_frequencies) / np.sqrt(
        expected_variance + EPSILON
    )

    return (standardized_genotypes @ standardized_genotypes.T) / standardized_genotypes.shape[1]


def reml_vc(y: np.ndarray, K: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Estimate variance components using a simple REML likelihood.

    Parameters
    ----------
    y
        Phenotype vector of length ``n_individuals``.
    K
        Genomic relationship matrix of shape ``(n_individuals, n_individuals)``.

    Returns
    -------
    tuple
        ``(variance_components, eigvals, U)`` where:

        - ``variance_components`` is ``[sigma_g2, sigma_e2]``;
        - ``eigvals`` are the eigenvalues of ``K``;
        - ``U`` contains the eigenvectors of ``K``.

    Notes
    -----
    The optimization is performed on log-variance parameters so that the returned
    variance components are positive. The eigendecomposition is returned because
    downstream whitening can reuse it without recomputing it.
    """
    phenotype = np.asarray(y, dtype=np.float64)
    grm = np.asarray(K, dtype=np.float64)

    if phenotype.ndim != 1:
        raise ValueError("y must be a 1D phenotype vector.")
    if grm.shape != (phenotype.size, phenotype.size):
        raise ValueError(
            "K must be a square matrix with one row/column per phenotype value."
        )

    eigvals, eigenvectors = np.linalg.eigh(grm)
    rotated_phenotype = eigenvectors.T @ phenotype

    def negative_reml_log_likelihood(theta: np.ndarray) -> float:
        """REML objective in the eigenspace of the GRM."""
        theta = np.clip(theta, *THETA_BOUNDS)
        sigma_g2, sigma_e2 = np.exp(theta)

        variance = sigma_g2 * eigvals + sigma_e2
        variance = np.maximum(variance, MIN_VARIANCE)

        return 0.5 * (
            np.sum(np.log(variance))
            + np.sum(rotated_phenotype**2 / variance)
            + np.log(np.sum(1.0 / variance))
        )

    result = minimize(
        negative_reml_log_likelihood,
        x0=np.log([1.0, 1.0]),
        method="L-BFGS-B",
    )

    if not result.success:
        raise RuntimeError(f"REML optimization failed: {result.message}")

    variance_components = np.exp(result.x)
    return variance_components, eigvals, eigenvectors
