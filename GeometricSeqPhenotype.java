import java.util.Arrays;

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
    private char[] nucleotideSequence;

    /**
     * Run expensive tests iff DEBUG == true.
     */
    public static final boolean DEBUG = false;

    // Abstraction Function:
    // A GeometricSeqPhenotype, gsp, is null if gsp.sequence = null,
    // otherwise gsp.sequence = sequence for sequence.length() > 0 and sequence.length() % 3 == 0
    // gsp1.sequence and gsp2.sequence are allowed to be equal for two given GeometricSeqPhenotype

    // Representation invariant for every GeometricSeqPhenotype:
    // gsp.sequence != null && gsp.sequence.length() > 0 && gsp.sequence.length() % 3 == 0 && gsp.sequence.length() == startingSequence.length()
    // for all indexSite such that gsp.sequence.charAt(i):
    //     indexSite is a character in ALPHABET
    // for all codon such that gsp.sequence.charAt(i) + gsp.sequence.charAt(i + 1) + gsp.sequence.charAt(i + 2):
    //     codon is not "STOP"

    /**
     * Constructor that creates a new GeometricSeqPhenotype.
     *
     * @effects Constructs a new GeometricSeqPhenotype that has a random starting sequence of length Parameters.sequence.length()
     * and an empty GeometricPhenotype
     */
    public GeometricSeqPhenotype() {
        this.correspondingGeometricPhenotype = new GeometricPhenotype();

        this.nucleotideSequence = new char[Parameters.startingSequence.length()];
        startingSequenceGenerator();
        checkRep();
    }

    /**
     * Constructor that creates a new GeometricSeqPhenotype.
     *
     * @param tA
     * @param tB
     * @effects Constructs a new GeometricSeqPhenotype that has a random starting sequence of length Parameters.sequence.length()
     * and GeometricPhenotype with the data content of the given parameters.
     */
    public GeometricSeqPhenotype(double tA, double tB) {
        this.correspondingGeometricPhenotype = new GeometricPhenotype(tA, tB);

        this.nucleotideSequence = new char[Parameters.startingSequence.length()];
        startingSequenceGenerator();
        checkRep();
    }

    // Generates a random sequence of nucleotides of length Parameters.startingSequence.length()
    private void startingSequenceGenerator() {
        for (int siteIndex = 0; siteIndex < Parameters.startingSequence.length(); siteIndex += 3) {
            String aminoAcid = "";

            while (aminoAcid.equals("") || aminoAcid.equals("STOP")) {
                String codon = "";
                for (int codonIndex = 0; codonIndex < 3; codonIndex++) {
                    int indexAlphabet = random(this.ALPHABET.length);
                    char nucleotide = this.ALPHABET[indexAlphabet];
                    codon += nucleotide;
                    this.nucleotideSequence[siteIndex + codonIndex] = nucleotide;
                }
                aminoAcid = Simulation.codonMap.get(codon);
            }
        }
        System.out.println(String.valueOf(this.nucleotideSequence));
    }

    /**
     * Constructor that creates a new GeometricSeqPhenotype.
     *
     * @param tA
     * @param tB
     * @param startingSequence
     * @requires nucleotideSequence != null && nucleotideSequence.length() > 0 && nucleotideSequence.length() % 3 == 0
     * @effects Constructs a new GeometricSeqPhenotype that has a starting sequence
     * and GeometricPhenotype with the data content of the given parameters.
     */
    public GeometricSeqPhenotype(double tA, double tB, char[] startingSequence) {
        this.correspondingGeometricPhenotype = new GeometricPhenotype(tA, tB);
        this.nucleotideSequence = startingSequence;
        checkRep();
    }

    public double getTraitA() {
        return correspondingGeometricPhenotype.getTraitA();
    }

    public double getTraitB() {
        return correspondingGeometricPhenotype.getTraitB();
    }

    public void setTraitA(double tA) {
        correspondingGeometricPhenotype.setTraitA(tA);
        checkRep();
    }

    public void setTraitB(double tB) {
        correspondingGeometricPhenotype.setTraitB(tB);
        checkRep();
    }

    public String getSequence() {
        return Arrays.toString(this.nucleotideSequence);
    }

    public GeometricPhenotype getCorrespondingGeometricPhenotype() {
        return this.correspondingGeometricPhenotype;
    }

    public double distance(Phenotype p) {
        GeometricPhenotype geometricP = ((GeometricSeqPhenotype) p).getCorrespondingGeometricPhenotype();
        checkRep();
        return correspondingGeometricPhenotype.distance(geometricP);
    }

    public double riskOfInfection( Phenotype[] history) {
        GeometricPhenotype[] geometricHistory = new GeometricPhenotype[history.length];
        for (int i = 0; i < history.length; i++) {
            geometricHistory[i] = ((GeometricSeqPhenotype) history[i]).getCorrespondingGeometricPhenotype();
        }
        checkRep();
        return correspondingGeometricPhenotype.riskOfInfection(geometricHistory);
    }

    /**
     * Returns a mutated copy of this GeometricSeqPhenotype (point substitution), original GeometricSeqPhenotype is unharmed
     *
     * @return a mutated copy of this GeometricSeqPhenotype
     */
    public Phenotype mutate() {
        // Implementation:
        // Mutates a nucleotide in the sequence at a random index
        // If mutation creates a stop codon, then mutate at another index

        // this.nucleotideSequence is a char[] instead of
        //  - String since Strings are immutable. New Strings have to be created to mutate a sequence, which is messy and uses up heap space. Plus, Strings include an extra byte '\n'.
        //  - StringBuffer because it's just a Wrapper around a char[], and has functionality that we don't need such as resizing.

        String wildTypeAminoAcid = "", mutantAminoAcid = "";
        int mutationIndexSiteSequence = -1;
        char randomMutatedNucleotide = '@';

        while (mutantAminoAcid.equals("") || mutantAminoAcid.equals("STOP")) {
            mutationIndexSiteSequence = random(this.nucleotideSequence.length);

            int codonStartIndex = 3 * (mutationIndexSiteSequence / 3);
            String originalCodon = "" + this.nucleotideSequence[codonStartIndex] + this.nucleotideSequence[codonStartIndex + 1] + this.nucleotideSequence[codonStartIndex + 2];

            wildTypeAminoAcid = Simulation.codonMap.get(originalCodon);

            int mutationIndexSiteCodon = mutationIndexSiteSequence % 3;
            char originalNucleotideToMutate = this.nucleotideSequence[mutationIndexSiteSequence];
            double[] transitionTransversion = Simulation.transitionTranversionProbability.get(originalNucleotideToMutate);


            double randomNum = Math.random();
            int indexAlphabet = 0;
            for (int i = 0; i < 4; i++) {
                if (randomNum < transitionTransversion[i]) {
                    indexAlphabet = i;
                    break;
                }
            }

            randomMutatedNucleotide = Parameters.AlphabetType.NUCLEOTIDES.getValidCharacters().charAt(indexAlphabet);
            StringBuffer mutatedCodon = new StringBuffer(originalCodon);
            mutatedCodon.setCharAt(mutationIndexSiteCodon, randomMutatedNucleotide);
            mutantAminoAcid = Simulation.codonMap.get(mutatedCodon.toString());
        }

        this.nucleotideSequence[mutationIndexSiteSequence] = randomMutatedNucleotide;

        //

        int mSiteMutationVectors = Parameters.AlphabetType.AMINO_ACIDS.getValidCharacters().indexOf(wildTypeAminoAcid);
        int nSiteMutationVectors = Parameters.AlphabetType.AMINO_ACIDS.getValidCharacters().indexOf(mutantAminoAcid);

        double[] mutations = Simulation.siteMutationVectors.get(mutationIndexSiteSequence)[mSiteMutationVectors][nSiteMutationVectors];

        this.correspondingGeometricPhenotype = new GeometricPhenotype(this.correspondingGeometricPhenotype.getTraitA() + mutations[0], this.correspondingGeometricPhenotype.getTraitB() + mutations[1]);

        checkRep();

        return new GeometricSeqPhenotype(this.correspondingGeometricPhenotype.getTraitA(), this.correspondingGeometricPhenotype.getTraitB(), this.nucleotideSequence);
    }

    /**
     * Returns the sequence and the position of this GeometricSeqPhenotype (i.e. the String representation of the GeometricSeqPhenotype represented by this).
     * Valid example outputs include "ACG, 0.0, 0.0" and "ACGTGTACGTGT, 2.3, 8.9".
     *
     * @return the String representation of the GeometricSeqPhenotype represented by this.
     */
    public String toString() {
        String fullString = String.format("%s, %.4f,%.4f", String.valueOf(this.nucleotideSequence), this.correspondingGeometricPhenotype.getTraitA(),
                                          this.correspondingGeometricPhenotype.getTraitB());
        return fullString;
    }

    // Returns a random index within the bounds of String.length.()=maxRange
    private int random(int maxRange) {
        Random random = new Random();
        return random.nextInt(0, maxRange - 1);
    }

    /**
     * Throws an exception if the representation invariant is violated.
     */
    private void checkRep() {
        if (DEBUG)  {
            assert(this.nucleotideSequence.length == Parameters.startingSequence.length()) : "Nucleotide sequences must remain the same length throughout a Simulation";
            assert(this.correspondingGeometricPhenotype != null) : "All GeometricSeqPhenotypes must have a GeometricPhenotype";
            for (int i = 0; i < this.nucleotideSequence.length; i += 3) {
                String triplet = "" + this.nucleotideSequence[i] + this.nucleotideSequence[i + 1] + this.nucleotideSequence[i + 2];
                String translatedAminoAcid = Simulation.codonMap.get(triplet);

                assert(!translatedAminoAcid.equals("STOP")) : "There should not be a stop codon at site " + (i/3);
            }
        }
    }
}
