# ARGWAS

ARGWAS is a Python tool for association testing using ancestral recombination graphs (ARGs). Instead of testing only observed SNPs, ARGWAS generates putative mutation markers from ARG branches and tests those markers for phenotype association.

The motivating idea is that causal QTLs may not be directly observed as SNPs. They may instead be in linkage disequilibrium with observed variation. By reconstructing the evolutionary history of sampled individuals, an ARG can suggest branch-level putative mutations that can be encoded as genotype-like vectors and tested in a GWAS-like framework.

## Current status

This package is an early research release. The code has been reorganized from a script into an installable package while preserving the original workflow and function names where possible.

## Installation

### Recommended on HPC clusters with conda

From a local checkout:

```bash
module load conda
source $(conda info --base)/etc/profile.d/conda.sh
conda env create -f environment.yml
conda activate argwas_env
argwas --help
```

For later runs, activate the same environment again:

```bash
module load conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate argwas_env
```

### Alternative pip installation

```bash
pip install -e .
```

## Command-line usage

```bash
argwas run \
  --output results \
  --ped input.ped \
  --map input.map \
  --trees trees_directory \
  --ances ancestral_alleles.txt \
  --workers 4
```

For backward compatibility, the old direct form is also supported:

```bash
argwas \
  --output results \
  --ped input.ped \
  --map input.map \
  --trees trees_directory \
  --ances ancestral_alleles.txt \
  --workers 4
```

The old source-checkout workflow is also kept:

```bash
python main.py -o results -p input.ped --map input.map -t trees_directory -a ancestral_alleles.txt -w 4
```

## Inputs

ARGWAS currently expects:

- a PLINK-style `.ped` file containing individual metadata, phenotype, and observed SNP genotypes;
- the matching PLINK `.map` file, either supplied with `--map` or inferred from the PED prefix;
- a directory containing `.trees` files readable by `tskit`;
- either a text file containing the genome-wide ancestral allele sequence for the MAP coordinate system, or `--ancestral-mode arg` to infer ancestral alleles from the ARG topology.


### Ancestral allele modes

By default, ARGWAS uses an explicit ancestral sequence:

```bash
argwas run \
  --output results \
  --ped input.ped \
  --map input.map \
  --trees trees_directory \
  --ances ancestral_alleles.txt \
  --ancestral-mode file
```

For datasets where the true ancestral sequence is unavailable, ARGWAS can also
infer ancestral alleles from the ARG:

```bash
argwas run \
  --output results \
  --ped input.ped \
  --map input.map \
  --trees trees_directory \
  --ancestral-mode arg \
  --ambiguous-ancestral skip
```

ARG-based inference is conservative. For each observed SNP, ARGWAS looks at the
local tree covering the SNP position. An allele is treated as derived only if the
haplotypes carrying that allele exactly form a clade in the local tree. The other
allele is then used as the ancestral allele. If neither allele, or both alleles,
fit cleanly, the site is marked as ambiguous.

Ambiguous sites are controlled by:

```text
--ambiguous-ancestral skip       omit ambiguous SNPs from the real-SNP GRM step
--ambiguous-ancestral major      use the major observed allele as a fallback
--ambiguous-ancestral reference  use the first observed PED allele as a fallback
```

The default is `skip`, because it avoids silently mixing guessed ancestral states
with inferred ones. When `--ancestral-mode arg` is used, ARGWAS writes:

```text
ancestral_inference.csv
```

with one row per observed SNP and columns describing the inferred allele,
status (`inferred`, `ambiguous`, `missing`, `fallback_major`, etc.), and whether
the SNP was used for whitening. SNPs with unknown ancestral state are skipped
only from the real-SNP whitening GRM; ARG-derived putative QTL markers are still
generated and tested independently.

### Ancestral allele coordinates

ARGWAS uses the physical position column of the PLINK MAP file to select the
ancestral allele for each observed SNP used to build the whitening GRM. By
default, MAP positions are treated as standard one-based genomic coordinates:

```text
MAP BP = 1  -> ancestral sequence index 0
MAP BP = 36 -> ancestral sequence index 35
```

This is the expected convention for PLINK/VCF-derived files. If your simulated
MAP file was generated with zero-based positions, run with:

```bash
argwas run ... --ancestral-coordinate-system zero-based
```

If `--map` is omitted, ARGWAS looks for a MAP file with the same prefix as the
PED file, for example `input.ped` -> `input.map`. The same coordinate convention
is also used when locating SNPs in tree sequences for `--ancestral-mode arg`: with
`one-based`, MAP BP=1 is converted to tree-sequence coordinate 0.

## Outputs

The main outputs are written to the output directory:

- `results.csv`: best ARG-derived marker per local tree, with `chromosome`, `marker_id`, `position`, `r_squared`, and `beta`;
- `results.map`: MAP file for selected ARG-derived markers;
- `results.ped`: PED file containing selected ARG-derived markers encoded as biallelic markers.

Intermediate files such as whitening matrices and temporary HDF5 marker stores are also written to the output directory.

### Association summary columns

For each local tree, ARGWAS keeps the candidate ARG-derived marker with the largest `r_squared`. The reported `beta` is the single-marker regression effect estimate after whitening. The reported `r_squared` is computed as the marker regression sum of squares divided by the whitened phenotype sum of squares, so it is a coefficient of determination rather than the previous raw association score.

## Optional LD clumping with PLINK

ARGWAS tests ARG-derived markers individually. Several correlated markers can therefore tag the same underlying QTL and share the apparent effect. The full `results.csv` file should be kept as the primary scan output, but it can be useful to create a second, LD-clumped candidate-signal list.

ARGWAS provides a helper to prepare a PLINK-compatible clumping association file:

```bash
argwas prepare-clump results/results.csv \
  --output results/results_for_clump.txt
```

This writes a ranking-only pseudo p-value:

```text
P = 1 - r_squared
```

This `P` column is only intended to let PLINK rank markers during clumping. It is not a calibrated p-value and should not be interpreted as statistical significance.

Example PLINK command:

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

See `docs/clumping.md` for more detail.

## Package layout

```text
src/argwas/
├── cli.py          # command-line interface
├── clump.py        # helper for PLINK clumping input files
├── ancestry.py     # ancestral allele lookup / ARG-based inference
├── pipeline.py     # complete ARGWAS workflow
├── arg_markers.py  # ARG branch/node to genotype-vector logic
├── hdf5_io.py      # streaming marker storage
└── whitening.py    # GRM, REML, and whitening utilities
```

## Scientific assumptions to review

The current ARG marker generation assumes that diploid individuals are represented by consecutive haploid sample nodes, so sample nodes `0` and `1` are individual 0, sample nodes `2` and `3` are individual 1, and so on. This matches the original script, but should be documented clearly for users and may need to be generalized later using `tskit` individual metadata.

The current whitening model is intentionally minimal and should be reviewed before publication, especially regarding fixed effects such as an intercept and any additional covariates.

## Development

Run basic syntax checks with:

```bash
python -m compileall src
```

Suggested next development steps:

1. add small example input files;
2. add tests for marker generation, HDF5 writing/reading, and PED conversion;
3. document the statistical model in `docs/method.md`;
4. add a citation section once a manuscript or preprint is available.
