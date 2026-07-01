# ARGWAS method notes

ARGWAS tests genotype-like markers induced by branches of local trees in an ancestral recombination graph.

At a high level:

1. observed SNP genotypes are used to estimate a genomic relationship matrix;
2. a variance-component model is fit from the GRM and phenotype;
3. phenotypes and ARG-derived marker vectors are whitened using the fitted covariance structure;
4. candidate branch markers are scanned tree by tree;
5. the marker with the highest whitened-space coefficient of determination, `r_squared`, is retained for each local tree and exported as PLINK-style output.

For a whitened phenotype vector `y` and whitened marker vector `x`, the effect estimate is `beta = (x'y) / (x'x)`. The raw association score used internally is the regression sum of squares, `(x'y)^2 / (x'x)`. ARGWAS reports `r_squared = ((x'y)^2 / (x'x)) / (y'y)`, i.e. the fraction of whitened phenotype sum of squares explained by that one marker.

This document should be expanded with the formal model, assumptions, and notation before release.

## Association output and clumping

ARGWAS reports the best ARG-derived marker per local tree using `r_squared` and `beta`. The `r_squared` column is a coefficient of determination computed after whitening; larger values indicate stronger single-marker association. The `beta` column is the corresponding single-marker effect estimate after whitening.

Because ARGWAS tests markers individually, multiple ARG-derived markers in linkage disequilibrium may capture the same underlying QTL. For interpretation, it can be useful to perform LD clumping as a post-processing step. Clumping is not applied automatically because it changes the question being answered: the raw scan reports all tested marker-level signals, while the clumped result is a summarized list of approximately independent candidate signals.

The `argwas prepare-clump` helper writes a PLINK association file with `P = 1 - r_squared`. This `P` value is only a ranking variable for PLINK clumping and is not a calibrated p-value.


## Ancestral allele lookup

For observed SNPs in the input PED/MAP files, ARGWAS defines the derived-allele
dosage relative to the ancestral allele at the SNP physical position. The MAP
file is therefore required for correct whitening: SNP column `j` in the PED is
matched to row `j` in the MAP, and the MAP BP coordinate is used to index the
ancestral sequence.

By default, BP coordinates are interpreted as one-based, so BP=1 corresponds to
Python index 0 in the ancestral sequence. This matches standard PLINK/VCF-style
coordinates. For zero-based simulated coordinates, use
`--ancestral-coordinate-system zero-based`.

This avoids the incorrect shortcut of using `ancestral_alleles[j]`, which only
works if the PED contains every genomic site in order from the start of the
ancestral sequence.


## ARG-based ancestral allele inference

For simulated data, the true ancestral allele at each genomic coordinate may be
known exactly. For real data, ancestral states are often uncertain and may depend
on outgroups or ancestral-state reconstruction. ARGWAS therefore supports an
optional ARG-based ancestral inference mode.

For each observed SNP in the input PED/MAP files, ARGWAS identifies the local
tree covering the SNP position. It then compares the carrier haplotypes of each
observed allele with the clades in that tree. If exactly one allele is carried by
haplotypes that form an exact clade, that allele is interpreted as the derived
state and the other allele is used as ancestral. If both alleles or neither allele
can be oriented in this way, the SNP is marked as ambiguous.

This is deliberately conservative. The method should not be interpreted as a
fully general ancestral-state reconstruction algorithm. It is a practical option
for polarizing observed SNPs used in the whitening GRM when a trusted ancestral
sequence is unavailable. The per-SNP calls are written to
`ancestral_inference.csv` so users can inspect how many alleles were inferred,
ambiguous, missing, skipped from whitening, or filled by a fallback policy.
Unknown ancestral states affect only the observed SNP panel used to build the
whitening GRM. They do not remove ARG-derived putative QTL markers from the
association scan, because those markers are already oriented by the ARG branch
on which the putative mutation occurs.

Recommended first-pass usage is:

```bash
argwas run ... --ancestral-mode arg --ambiguous-ancestral skip
```

Fallback modes such as `--ambiguous-ancestral major` are available for exploratory
work, but they should be reported explicitly because they mix ARG-inferred states
with approximate allele-frequency-based states.
