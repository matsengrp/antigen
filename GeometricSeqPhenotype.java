import org.biojava.nbio.core.sequence.transcription.*;

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
        this.correspondingGeometricPhenotype = new GeometricPhenotype();
        this.nucleotideSequence = startingSequenceGenerator();

        this.aminoAcidSequence = translateSequenceToAminoAcids(this.nucleotideSequence);
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

        this.aminoAcidSequence = translateSequenceToAminoAcids(this.nucleotideSequence);
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

        this.aminoAcidSequence = translateSequenceToAminoAcids(this.nucleotideSequence);
    }

    /*
    public GeometricSeqPhenotype(String startingSequence, Phenotype correspondingGeometricPhenotype) {
        startingSequence = startingSequence.toUpperCase();
        this.sequence = startingSequence;
        this.correspondingGeometricPhenotype = (GeometricPhenotype) correspondingGeometricPhenotype;
    }
    */

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
        int indexSite = random(this.nucleotideSequence.length());
        int indexAlphabet = random(this.ALPHABET.length);

        // substitute a random index of sequence with a random character in the ALPHABET
        StringBuilder mutated = new StringBuilder(this.nucleotideSequence);
        mutated.setCharAt(indexSite, this.ALPHABET[indexAlphabet]);

        this.nucleotideSequence = mutated.toString();
        this.aminoAcidSequence = translateSequenceToAminoAcids(this.nucleotideSequence);

        GeometricPhenotype mutatedCorrespondingGeometricPhenotype = (GeometricPhenotype) correspondingGeometricPhenotype.mutate();

        System.out.println(this.aminoAcidSequence + "       " + this.nucleotideSequence);

        // TODO: GeometricSeqPhenotype(String startingSequence, Phenotype correspondingGeometricPhenotype) instead?
        return new GeometricSeqPhenotype(mutatedCorrespondingGeometricPhenotype.getTraitA(),
                                         mutatedCorrespondingGeometricPhenotype.getTraitB(), mutated.toString());
    }

    private int random(int maxRange) {
        Random random = new Random();
        return random.nextInt(0, maxRange - 1);
    }

    public String toString() {
        return this.nucleotideSequence;
    }

    private String translateSequenceToAminoAcids(String nucleotideSequence) {
        StringBuilder translatedSequence = new StringBuilder();
        for (int i = 0; i < nucleotideSequence.length(); i += 3) {
            String triplet = nucleotideSequence.substring(i, i+3);
            String translatedAminoAcid = Simulation.codonMap.get(triplet);

            translatedSequence.append(translatedAminoAcid);
        }

        return translatedSequence.toString();
    }

}
