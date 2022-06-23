/**
 * <b>GeometricSeqPhenotype</b> represents a phenotype.
 * Distance and cross-immunity can be calculated between two GeometricSeqPhenotypes.
 * Multiple Viruses can reference a single GeometricSeqPhenotype which is identified by its data contents.
 *
 * @author Thien Tran
 */
public class GeometricSeqPhenotype implements Phenotype {
    /**
     * The valid letters that make up the nucleotide sequence of this GeometricSeqPhenotype.
     *
     * Valid letters are {A, G, T , C}
     */
    private final char[] ALPHABET = Parameters.AlphabetType.NUCLEOTIDES.getValidCharacters().toCharArray();

    /**
     * The corresponding GeometricPhenotype of this GeometricSeqPhenotype
     */
    private GeometricPhenotype correspondingGeometricPhenotype;

    /**
     * The nucleotide sequence of this GeometricSeqPhenotype
     */
    private String nucleotideSequence;

    /**
     * The amino acid sequence of this GeometricSeqPhenotype
     */
    private String aminoAcidSequence;

    /**
     * Constructor that creates a new GeometricSeqPhenotype.
     *
     * @spec.effects Constructs a new GeometricSeqPhenotype that has a random starting sequence of length Parameters.sequence.length()
     * and GeometricPhenotype
     */
    public GeometricSeqPhenotype() {
        this.correspondingGeometricPhenotype = new GeometricPhenotype(1, 1);

        this.nucleotideSequence = startingSequenceGenerator();
        translateSequenceToAminoAcids();

        while (this.aminoAcidSequence.contains("STOP")) {
            this.nucleotideSequence = startingSequenceGenerator();
            translateSequenceToAminoAcids();
        }

        if (this.aminoAcidSequence.contains("STOP")) {
            System.out.println("hdjfaf" + this.aminoAcidSequence);
        }
    }

    /**
     * Constructor that creates a new GeometricSeqPhenotype.
     *
     * @param tA
     * @param tB
     * @spec.requires nucleotideSequence != null && nucleotideSequence.length() > 0 && nucleotideSequence.length() % 3 == 0
     * @spec.effects Constructs a new GeometricSeqPhenotype that has a random starting sequence of length Parameters.sequence.length()
     * and GeometricPhenotype with the data content of the given parameters.
     */
    public GeometricSeqPhenotype(double tA, double tB) {
        this.correspondingGeometricPhenotype = new GeometricPhenotype(tA, tB);

        this.nucleotideSequence = startingSequenceGenerator();
        translateSequenceToAminoAcids();

        while (this.aminoAcidSequence.contains("STOP")) {
            this.nucleotideSequence = startingSequenceGenerator();
            translateSequenceToAminoAcids();
        }

        if (this.aminoAcidSequence.contains("STOP")) {
            System.out.println("hdjfaf" + this.aminoAcidSequence);
        }
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
        char originalNucleotide = this.nucleotideSequence.charAt(mutationIndexSite);
        char randomNucleotide = getRandomMutatedNucleotide(originalNucleotide, mutationIndexSite);


        int tripletStartIndex = 3 * (mutationIndexSite / 3);
        String originalAminoAcid = Simulation.codonMap.get(this.nucleotideSequence.substring(tripletStartIndex, tripletStartIndex + 3));

        // substitute a random index of sequence with a random character in the ALPHABET
        StringBuilder mutatedNucleotideSequence = new StringBuilder(this.nucleotideSequence);
        mutatedNucleotideSequence.setCharAt(mutationIndexSite, randomNucleotide);
        this.nucleotideSequence = mutatedNucleotideSequence.toString();


        String mutatedAminoAcid = translateTripleSequenceToAminoAcids(mutationIndexSite);

        int mSiteMutationVectors = Parameters.AlphabetType.AMINO_ACIDS.getValidCharacters().indexOf(originalAminoAcid);
        int nSiteMutationVectors = Parameters.AlphabetType.AMINO_ACIDS.getValidCharacters().indexOf(mutatedAminoAcid);

        if (mSiteMutationVectors == -1 || nSiteMutationVectors == -1) {
            System.out.println(originalAminoAcid);
            System.out.println(mutatedAminoAcid);
        }

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

        String codon = this.nucleotideSequence.substring(mutationAminoAcidSiteNumber * 3, (mutationAminoAcidSiteNumber * 3) + 3);
        String translatedAminoAcidFromTriplet = Simulation.codonMap.get(codon);

        StringBuilder updatedAminoAcids = new StringBuilder(this.aminoAcidSequence);
        updatedAminoAcids.setCharAt(mutationAminoAcidSiteNumber, translatedAminoAcidFromTriplet.charAt(0));

        this.aminoAcidSequence = updatedAminoAcids.toString();

        return translatedAminoAcidFromTriplet;
    }

    private char getRandomMutatedNucleotide(char originalNucleotide, int mutationIndexSite) {
        int indexAlphabet = random(this.ALPHABET.length);
        char randomMutatedNucleotide = this.ALPHABET[indexAlphabet];

        int mutationAminoAcidSiteNumber = mutationIndexSite / 3;

        String codon = this.nucleotideSequence.substring(mutationAminoAcidSiteNumber * 3, (mutationAminoAcidSiteNumber * 3) + 3);

        StringBuilder newCodon = new StringBuilder(codon);
        newCodon.setCharAt(mutationIndexSite % 3, randomMutatedNucleotide);

        String translatedAminoAcidFromTriplet = Simulation.codonMap.get(newCodon.toString());

        while (originalNucleotide == randomMutatedNucleotide || translatedAminoAcidFromTriplet.equals("STOP")) {
            indexAlphabet = random(this.ALPHABET.length);
            randomMutatedNucleotide = this.ALPHABET[indexAlphabet];

            newCodon.setCharAt(mutationIndexSite % 3, randomMutatedNucleotide);

            translatedAminoAcidFromTriplet = Simulation.codonMap.get(newCodon.toString());
        }

        return randomMutatedNucleotide;
    }

    private int random(int maxRange) {
        Random random = new Random();
        return random.nextInt(0, maxRange - 1);
    }
}
