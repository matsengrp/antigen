/**
 * <b>GeometricSeqPhenotype</b> represents a phenotype.
 * GeometricSeqPhenotypes are identified by their data contents. Each instance has a GeometricPhenotype.
 * GeometricSeqPhenotypes permit all data, except for <b>""</b> and <b>null</b>.
 *
 * @author Thien Tran
 */
public class GeometricSeqPhenotype implements Phenotype {
    /**
     * The valid letters that make up the sequence of this SequencePhenotype.
     *
     * Valid letters depend on the representation of the sequence (nucleotides or amino acids)
     */
    private final char[] ALPHABET = Parameters.alphabet.toCharArray();

    /**
     * The corresponding GeometricPhenotype of this GeometricSeqPhenotype
     */
    private GeometricPhenotype correspondingGeometricPhenotype;

    /**
     * The sequence of this SequencePhenotype
     */
    private String nucleotideSequence;

    private String aminoAcidSequence;

    /**
     * Constructor that creates a new GeometricSeqPhenotype.
     *
     * @spec.requires sequence != null && sequence.length() > 0
     * @spec.effects Constructs a new GeometricSeqPhenotype that has a random starting sequence of length Parameters.sequence.length()
     * and GeometricPhenotype
     */
    public GeometricSeqPhenotype() {
        this.correspondingGeometricPhenotype = new GeometricPhenotype(1, 1);
        this.nucleotideSequence = startingSequenceGenerator();

        translateSequenceToAminoAcids();
    }

    /**
     * Constructor that creates a new GeometricSeqPhenotype.
     *
     * @param tA
     * @param tB
     * @spec.requires sequence != null && sequence.length() > 0
     * @spec.effects Constructs a new GeometricSeqPhenotype that has a random starting sequence of length Parameters.sequence.length()
     * and GeometricPhenotype with the data content of the given parameters.
     */
    public GeometricSeqPhenotype(double tA, double tB) {
        this.correspondingGeometricPhenotype = new GeometricPhenotype(tA, tB);

        this.nucleotideSequence = startingSequenceGenerator();
        translateSequenceToAminoAcids();
    }

    /**
     * Constructor that creates a new GeometricSeqPhenotype.
     *
     * @param tA
     * @param tB
     * @param startingSequence
     * @spec.requires startingSequence != null && sequence.length() > 0
     * @spec.effects Constructs a new GeometricSeqPhenotype that has a starting sequence
     * and GeometricPhenotype with the data content of the given parameters.
     */
    public GeometricSeqPhenotype(double tA, double tB, String startingSequence) {
        this.correspondingGeometricPhenotype = new GeometricPhenotype(tA, tB);

        startingSequence = startingSequence.toUpperCase();
        this.nucleotideSequence = startingSequence;

        translateSequenceToAminoAcids();
    }

    private void translateSequenceToAminoAcids() {
        StringBuilder translatedSequence = new StringBuilder();
        for (int i = 0; i < nucleotideSequence.length(); i += 3) {
            String triplet = nucleotideSequence.substring(i, i+3);
            String translatedAminoAcid = Simulation.codonMap.get(triplet);

            translatedSequence.append(translatedAminoAcid);
        }

        this.aminoAcidSequence = translatedSequence.toString();
    }

    // Generates a random sequence of nucleotides or amino acids of length Parameters.startingSequence.length()
    private String startingSequenceGenerator() {
        String startingSequence = "";
        for (int i = 0; i < Parameters.startingSequence.length(); i++) {
            int indexAlphabet = random(this.ALPHABET.length);
            startingSequence += this.ALPHABET[indexAlphabet];
        }
        return startingSequence;
    }

    public double getTraitA() {
        return correspondingGeometricPhenotype.getTraitA();
    }

    public double getTraitB() {
        return correspondingGeometricPhenotype.getTraitB();
    }

    public void setTraitA(double tA) {
        correspondingGeometricPhenotype.setTraitA(tA);
    }

    public void setTraitB(double tB) {
        correspondingGeometricPhenotype.setTraitB(tB);
    }

    public String getSequence() {
        return this.nucleotideSequence;
    }

    public GeometricPhenotype getCorrespondingGeometricPhenotype() {
        return this.correspondingGeometricPhenotype;
    }

    public double distance(Phenotype p) {
        GeometricPhenotype geometricP = ((GeometricSeqPhenotype) p).getCorrespondingGeometricPhenotype();
        return correspondingGeometricPhenotype.distance(geometricP);
    }

    public double riskOfInfection( Phenotype[] history) {
        GeometricPhenotype[] geometricHistory = new GeometricPhenotype[history.length];
        for (int i = 0; i < history.length; i++) {
            geometricHistory[i] = ((GeometricSeqPhenotype) history[i]).getCorrespondingGeometricPhenotype();
        }
        return correspondingGeometricPhenotype.riskOfInfection(geometricHistory);
    }

    public Phenotype mutate() {
        int mutationIndexSite = random(this.nucleotideSequence.length());
        int indexAlphabet = random(this.ALPHABET.length);

        int tripletStartIndex = 3 * (mutationIndexSite / 3);

        String originalAminoAcid = Simulation.codonMap.get(this.nucleotideSequence.substring(tripletStartIndex, tripletStartIndex + 3));

        // substitute a random index of sequence with a random character in the ALPHABET
        StringBuilder mutated = new StringBuilder(this.nucleotideSequence);
        mutated.setCharAt(mutationIndexSite, this.ALPHABET[indexAlphabet]);

        String mutatedNucleotideSequence = mutated.toString();

        this.nucleotideSequence = mutatedNucleotideSequence;
        String mutatedAminoAcid = translateTripleSequenceToAminoAcids(mutationIndexSite);

        int mSiteMutationVectors = Parameters.AlphabetType.AMINO_ACIDS.getValidCharacters().indexOf(originalAminoAcid);
        int nSiteMutationVectors = Parameters.AlphabetType.AMINO_ACIDS.getValidCharacters().indexOf(mutatedAminoAcid);

        double[] mutations = Simulation.siteMutationVectors.get(mutationIndexSite)[mSiteMutationVectors][nSiteMutationVectors];

        this.correspondingGeometricPhenotype = new GeometricPhenotype(this.correspondingGeometricPhenotype.getTraitA() + mutations[0], this.correspondingGeometricPhenotype.getTraitB() + mutations[1]);

        return new GeometricSeqPhenotype(this.correspondingGeometricPhenotype.getTraitA(), this.correspondingGeometricPhenotype.getTraitB(), this.nucleotideSequence);
    }

    public String toString() {
        String fullString = String.format("%s, %.4f,%.4f", this.nucleotideSequence, this.correspondingGeometricPhenotype.getTraitA(),
                                          this.correspondingGeometricPhenotype.getTraitB());
        return fullString;
    }

    /**
     * Translates a codon, a sequence of three DNA or RNA nucleotides, to its corresponding amino acid
     * until it's not a stop signal.
     *
     * @param mutationIndexSite the overall index of the nucleotide to be mutated
     * @return the specified protein after mutating the nucleotide at mutationIndexSite in nucleotideSequence
     */
    private String translateTripleSequenceToAminoAcids(int mutationIndexSite) {
        int mutationAminoAcidSiteNumber = mutationIndexSite / 3;
        int mutationTripleStartIndex = mutationIndexSite % 3;

        String codon = this.nucleotideSequence.substring(mutationAminoAcidSiteNumber * 3, (mutationAminoAcidSiteNumber * 3) + 3);
        String translatedAminoAcidFromTriplet = Simulation.codonMap.get(codon);

        boolean nucleotideSequenceChanged = false;

        // Mutate the nucleotide at the original mutationIndexSite until the amino acid is not "STOP"
        while (translatedAminoAcidFromTriplet.equals("STOP")) {
            int indexAlphabet = random(this.ALPHABET.length);

            StringBuilder mutated = new StringBuilder(codon);
            mutated.setCharAt(mutationTripleStartIndex, this.ALPHABET[indexAlphabet]);

            codon = mutated.toString();
            translatedAminoAcidFromTriplet = Simulation.codonMap.get(codon);
            nucleotideSequenceChanged = true;
        }

        // If translatedAminoAcidFromTriplet was ever "STOP" then update this.nucleotideSequence
        if (nucleotideSequenceChanged) {
            if (mutationAminoAcidSiteNumber == 0) {
                // first site
                this.nucleotideSequence = codon + this.nucleotideSequence.substring(3);
            } else if (mutationAminoAcidSiteNumber == (this.nucleotideSequence.length() / 3) - 1) {
                // last site
                this.nucleotideSequence = this.nucleotideSequence.substring(0, mutationAminoAcidSiteNumber * 3) + codon;
            } else {
                this.nucleotideSequence = this.nucleotideSequence.substring(0, 3) + codon + this.nucleotideSequence.substring((mutationAminoAcidSiteNumber + 1) * 3);
            }
        }

        StringBuilder updatedAminoAcids = new StringBuilder(this.aminoAcidSequence);
        updatedAminoAcids.setCharAt(mutationAminoAcidSiteNumber, translatedAminoAcidFromTriplet.charAt(0));

        this.aminoAcidSequence = updatedAminoAcids.toString();

        return translatedAminoAcidFromTriplet;
    }

    private int random(int maxRange) {
        Random random = new Random();
        return random.nextInt(0, maxRange - 1);
    }
}
