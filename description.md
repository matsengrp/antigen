## GeometricSeqPhenotype.java
`GeometricSeqPhenotype` is a subclass of `GeometricPhenotype`. 
Along with the parameters that determines the position of the object in Euclidean space, it stores nucleotide sequence information and the number of epitope and non-epitope mutations in fields. 
Only the nucleotide sequence is being stored to improve space and time complexity. 
The translation of a nucleotide sequence to a protein sequence can be done after the simulation. 
However, this does not mean amino acids don't play an important role in our model. 
For each site in the sequence, a matrix of vectors is precomputed before the simulation runs. 
The vectors are drawn from a gamma distribution, whose parameters can be changed in `parameters.yml`.
The number of epitope sites, which range from 0 to (the starting sequence length / 3), is parameterized in that file as well.
Sites are chosen to be epitopes randomly and will remain epitopes until the program exits. 
Epitope sites have corresponding matrices with vectors drawn from a gamma distribution with different parameters than non-epiope sites. 
`testGammaDistribution.py` validates the graph shapes of each site. 
This validation is consistent and possible due to the precomputing of the vectors.
Even for over 200 sites, precomputing all the matrices are quite fast.

### mutate()
This method randomly selects an index in the nucleotide sequence to mutate.
Based on the transition and transversion ratio, a cumulative sum distribution array is used to determine what the nucleotide will mutate to. 
If the mutant amino acid is a stop codon, the process repeats starting with randomly selecting an index to mutate.
Once a valid mutation occurs, where the object moves in space is determined by the vectors in the site's matrix at entry $m$, $n$ which represents the index of the wild type and mutant amino acids in `Parameters.NUCLEOTIDES`, respectively. 
`this` does not update. 
Instead, a new `GeometricSeqPhenotype`with the updated nucleotide sequence, new position parameters, and number of epitope and non-epitope mutations.
There must be a mutation if this method is called, so the nucleotide sequence must be different and the number of epitope mutations xor non-epitope mutations must be updated by 1.
However, the new position parameters might not change since a nucleotide mutation does not necessarily cause a mutation in the protein sequence. 
