# LD clumping ARGWAS results

ARGWAS tests ARG-derived markers one at a time. This is useful for detecting candidate signals, but it means that several highly correlated markers can share the apparent effect of the same underlying QTL. For interpretation, it is often useful to keep the full `results.csv` file and create a second, LD-clumped set of approximately independent signals.

ARGWAS does not clump results by default because clumping is a post-processing and interpretation step. The raw output answers: "which ARG-derived markers explain phenotype variation?" A clumped output answers: "which approximately independent regions/signals remain after LD pruning?"

## Prepare a PLINK clumping file

ARGWAS reports `r_squared`, where larger values mean stronger association. PLINK `--clump` expects a field where smaller values are better, usually a p-value. ARGWAS therefore provides a helper that writes a ranking-only pseudo p-value:

```text
P = 1 - r_squared
```

This value is only used to rank markers for PLINK clumping. It is not a calibrated p-value and should not be interpreted as statistical significance.

Run:

```bash
argwas prepare-clump results/results.csv \
  --output results/results_for_clump.txt
```

The output contains columns compatible with PLINK clumping:

```text
CHR SNP BP P R2_ARGWAS BETA_ARGWAS
```

`SNP` matches the marker IDs in `results.map` and in PLINK files derived from `results.map`/`results.ped`.

## Run PLINK clumping

If you have converted the ARGWAS `results.map` and `results.ped` to a PLINK binary dataset called `results`, you can run for example:

```bash
plink \
  --bfile results \
  --clump results/results_for_clump.txt \
  --clump-field P \
  --clump-p1 1 \
  --clump-r2 0.1 \
  --clump-kb 500 \
  --out clumped
```

The exact `--clump-r2` and `--clump-kb` values are analysis choices. A stricter `--clump-r2` keeps fewer, more independent signals; a larger window merges signals over longer genomic distances.

## Recommended reporting

For reproducibility, report both:

1. the complete ARGWAS `results.csv`, and
2. the clumped signal list with the chosen PLINK parameters.

The clumped file should be treated as a summarized candidate-signal list, not as a replacement for the full ARGWAS scan.
