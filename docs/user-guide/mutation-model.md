# Mutation Model

antigen-prime couples genetic sequence evolution with antigenic phenotype changes. Each virus carries both a nucleotide sequence and coordinates in antigenic space. Mutation events simultaneously change the sequence and may move the virus in antigenic space.

## K80 Nucleotide Mutation

Mutations follow the K80 (Kimura 1980) model parameterized by a transition/transversion ratio $\kappa$.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `transitionTransversionRatio` | 5.0 | Bias toward transitions vs transversions |

Default $\kappa = 5.0$ matches empirical observations in influenza (Rabadan et al. 2006, Bloom & Glassman 2009).

During mutation:

1. A nucleotide site is randomly selected
2. A new nucleotide is drawn from the K80 probability distribution
3. If mutation creates a stop codon, reject and try another site

## Epitope vs Non-Epitope Sites

The initial sequence encodes a protein. Users define a subset of amino-acid sites as "epitope sites" (remaining sites are "non-epitope").

| Parameter | Default | Description |
|-----------|---------|-------------|
| `epitopeSites` | "epitopeSites.txt" | File listing epitope site indices (1-indexed) |

For influenza HA, we typically use the 49 epitope sites from Luksza & Lassig (2014).

## Antigenic Space Movement

The effect of a mutation depends on whether and how it changes the protein:

| Mutation Type | Antigenic Effect |
|---------------|------------------|
| Synonymous | No movement |
| Non-synonymous at epitope site | Large movement |
| Non-synonymous at non-epitope site | Very small movement |

### Step Size Distribution

Movement step sizes are drawn from gamma distributions:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `meanStepEpitope` | 0.6 | Mean step size for epitope mutations (antigenic units) |
| `sdStepEpitope` | 0.3 | Standard deviation for epitope mutations |
| `meanStep` | 1e-5 | Mean step size for non-epitope mutations |
| `sdStep` | 0.3 | Standard deviation for non-epitope mutations |

Default epitope step size of 0.6 antigenic units produces ~1.6 AU/year antigenic drift matching empirical HI assay data (Smith et al. 2004, Koel et al. 2013).

### Step Direction

The direction $\theta$ of movement in antigenic space is uniformly random.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `mut2D` | false | If true, allow full 360° arc; if false, 1D movement only |

## Acceptance Rates

Users can apply acceptance/rejection filtering to model selection:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `epitopeAcceptance` | 1.0 | Probability of accepting epitope mutations |
| `nonEpitopeAcceptance` | 1.0 | Probability of accepting non-epitope mutations |

Setting different rates allows modeling differential selection between site types. Note: only applied to non-synonymous mutations.

## High/Low Epitope Sites

For finer control, epitope sites can be subdivided into "high" and "low" categories with different step size distributions:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `proportionHighSites` | 0.2 | Fraction of epitope sites designated as "high" |
| `meanStepEpitopeLow` | 0.3 | Mean step for "low" epitope sites |
| `meanStepEpitopeHigh` | 0.3 | Mean step for "high" epitope sites |

## Predefined Vectors

By default, step sizes and directions are drawn randomly. Alternatively, users can predefine mutation vectors for each site/amino-acid pair:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `predefinedVectors` | false | Use precomputed site-specific mutation effects |

## Implementation

The mutation logic is implemented in `GeometricSeqPhenotype.mutate()`:

```
src/main/java/org/antigen/phenotype/GeometricSeqPhenotype.java
```

Key steps:

1. Select random nucleotide site
2. Draw mutant nucleotide from K80 distribution
3. Reject if creates stop codon
4. Determine if synonymous; if so, return with sequence change only
5. Apply acceptance filter based on site type
6. Update mutation counts (epitope/non-epitope)
7. Calculate antigenic step (predefined or random)
8. Return new phenotype with updated sequence and coordinates

## References

- Kimura M. (1980). A simple method for estimating evolutionary rates of base substitutions. J Mol Evol.
- Smith DJ et al. (2004). Mapping the antigenic and genetic evolution of influenza virus. Science.
- Luksza M, Lassig M. (2014). A predictive fitness model for influenza. Nature.
- Koel BF et al. (2013). Substitutions near the receptor binding site determine major antigenic change during influenza virus evolution. Science.
