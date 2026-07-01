"""Utilities for preparing ARGWAS results for LD clumping with PLINK."""

from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_RESULT_COLUMNS = {"chromosome", "position", "r_squared"}


def prepare_plink_clump_file(
    results_csv: str | Path,
    output: str | Path,
    *,
    marker_column: str = "marker_id",
    chromosome_column: str = "chromosome",
    position_column: str = "position",
    r_squared_column: str = "r_squared",
    beta_column: str = "beta",
) -> Path:
    """Write a PLINK-compatible association file for clumping ARGWAS results.

    PLINK's ``--clump`` command expects an association file containing a marker
    identifier and a field where smaller values represent stronger signals. ARGWAS
    reports ``r_squared`` values, where larger values represent stronger signals,
    so this helper writes a ranking-only pseudo p-value:

    ``P = 1 - r_squared``

    This value is only intended to rank markers during clumping. It is not a
    calibrated p-value and should not be interpreted as statistical significance.

    Parameters
    ----------
    results_csv
        Path to the ARGWAS ``results.csv`` file.
    output
        Path of the clumping association file to write.
    marker_column, chromosome_column, position_column, r_squared_column, beta_column
        Column names in ``results_csv``. The defaults match ARGWAS output.

    Returns
    -------
    pathlib.Path
        Path to the written clumping file.
    """

    results_csv = Path(results_csv)
    output = Path(output)

    results = pd.read_csv(results_csv)
    required = {chromosome_column, position_column, r_squared_column}
    missing = required.difference(results.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"Missing required column(s) in {results_csv}: {missing_text}")

    if marker_column not in results.columns:
        raise ValueError(
            f"Missing marker column '{marker_column}' in {results_csv}. "
            "Recent ARGWAS output includes a marker_id column matching results.map. "
            "For older results, add marker IDs before clumping or rerun ARGWAS."
        )

    r_squared = results[r_squared_column].to_numpy(dtype=float)
    if np.any(~np.isfinite(r_squared)):
        raise ValueError(f"Column '{r_squared_column}' contains NaN or infinite values")

    clipped_r_squared = np.clip(r_squared, 0.0, 1.0)
    pseudo_p = 1.0 - clipped_r_squared

    clump_table = pd.DataFrame(
        {
            "CHR": results[chromosome_column].astype(int),
            "SNP": results[marker_column].astype(str),
            "BP": results[position_column].astype(int),
            "P": pseudo_p,
            "R2_ARGWAS": results[r_squared_column].astype(float),
        }
    )

    if beta_column in results.columns:
        clump_table["BETA_ARGWAS"] = results[beta_column].astype(float)

    output.parent.mkdir(parents=True, exist_ok=True)
    clump_table.to_csv(output, sep="\t", index=False)
    return output
