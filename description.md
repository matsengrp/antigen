## GeometricSeqPhenotype.java
`GeometricSeqPhenotype` is a subclass of `GeometricPhenotype`.
Along with the parameters that determine the position of the virus in antigenic space,
it stores the virus's nucleotide sequence and the number of epitope and non-epitope mutations that have accumulated over the virus's evolutionary history.
A key question in developing this class was how to map changes in a virus's nucleotide sequence to changes in its location in antigenic space.
Our basic strategy for doing so is as follows.
First, before the simulation starts, we precompute antigenic effects for all amino-acid mutations at all sites in the protein.
As with the original version of `antigen`, mutation effects are encoded by vectors that move the virus in antigenic space, 
where vector angles are randomly drawn from a uniform distribution and vector magnitudes are randomly drawn from a gamma distribution.
In this updated version, we specify two gamma distributions: one for epitope sites and one for non-epitope sites, 
allowing average mutational effects to differ between the two groups of sites.
The number of epitope sites, which range from 0 to (the starting sequence length / 3), is parameterized in that file as well.
Sites are chosen to be epitopes randomly and will remain epitopes until the program exits.
`testGammaDistribution.py` can be used to plot the distribution of mutation effects for individual sites, 
helping to validate and visualize differences between epitope and non-epitope sites. 
This validation is consistent and possible due to the precomputing of the vectors.
Even for over 200 sites, precomputing all mutational effects is quite fast.

### mutate()
In the original version `antigen`, there is a rate at which individual viruses mutate in antigenic phenotype.
In our updated version, mutations also result in a change to the virus's nucleotide sequence.
Specifically, when there is a mutation event, the `mutate()` method randomly selects an index in the nucleotide sequence to mutate.
Based on the transition/transversion ratio, 
a cumulative sum distribution array is used to determine what the nucleotide will mutate to.
If the mutant amino acid is a stop codon, the mutation is discarded, and the process repeats starting with randomly selecting an index to mutate.
Once a valid mutation occurs, where the object moves in space is determined by precomputed mutational effect vectors described above.
`this` does not update. Instead, the code creates a new `GeometricSeqPhenotype` object with
updated entries for the virus's nucleotide sequence, location in antigenic space, and number of epitope and non-epitope mutations.
There must be a mutation if this method is called, so the nucleotide sequence must be different
and the number of epitope mutations or non-epitope mutations must be updated by 1.
However, the new position parameters might not change since a nucleotide mutation does not necessarily 
cause a mutation in the protein sequence. 

### Requirements
- DMS data has to have the same number of amino-acid sites as the input sequence
- sites must occur in sequential numbering (the site column in the DMS data is ignored)
- the file has a single header line (giving column names)
- amino acids are assumed to occur in a particular order (right? if so what is this order?)

epitope sites should be indexed at 1.0