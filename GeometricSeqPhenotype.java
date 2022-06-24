import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

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
     * Run expensive tests iff DEBUG == true.
     */
    public static final boolean DEBUG = true;

    // Abstraction Function:
    // A SequencePhenotype, s, is null if s.sequence = null, otherwise s.sequence = sequence
    // for sequence.length() > 0
    // s1.sequence and s2.sequence are allowed to be equal for two given SequencePhenotypes

    // Representation invariant for every GeometricSeqPhenotype:
    // gsp.sequence != null && gsp.sequence.length() > 0 && gsp.sequence.length() % 3 == 0
    // for all indexSite such that s.sequence.charAt(i):
    //     indexSite is a character in ALPHABET

    /**
     * Constructor that creates a new GeometricSeqPhenotype.
     *
     * @spec.effects Constructs a new GeometricSeqPhenotype that has a random starting sequence of length Parameters.sequence.length()
     * and GeometricPhenotype
     */
    public GeometricSeqPhenotype() {
        this.correspondingGeometricPhenotype = new GeometricPhenotype();

        this.nucleotideSequence = startingSequenceGenerator();
        String aminoAcidSequence = translateSequenceToAminoAcids();

        while (aminoAcidSequence.contains("STOP")) {
            this.nucleotideSequence = startingSequenceGenerator();
            aminoAcidSequence = translateSequenceToAminoAcids();
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
        String aminoAcidSequence = translateSequenceToAminoAcids();

        while (aminoAcidSequence.contains("STOP")) {
            this.nucleotideSequence = startingSequenceGenerator();
            aminoAcidSequence = translateSequenceToAminoAcids();
        }
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

    private String translateSequenceToAminoAcids() {
        StringBuilder translatedSequence = new StringBuilder();
        for (int i = 0; i < this.nucleotideSequence.length(); i += 3) {
            String triplet = this.nucleotideSequence.substring(i, i+3);
            String translatedAminoAcid = Simulation.codonMap.get(triplet);

            translatedSequence.append(translatedAminoAcid);
        }

        return translatedSequence.toString();
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
        // requires: startingSequence to not include a stop codon
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
        // index of sequence to mutate
        int mutationIndexSiteSequence = random(this.nucleotideSequence.length());


        int codonStartIndex = 3 * (mutationIndexSiteSequence / 3);
        String originalCodon = this.nucleotideSequence.substring(codonStartIndex, codonStartIndex + 3);
        String originalAminoAcid = Simulation.codonMap.get(originalCodon);

        String mutatedAminoAcid = getRandomMutatedNucleotide(originalCodon, mutationIndexSiteSequence);

        int mSiteMutationVectors = Parameters.AlphabetType.AMINO_ACIDS.getValidCharacters().indexOf(originalAminoAcid);
        int nSiteMutationVectors = Parameters.AlphabetType.AMINO_ACIDS.getValidCharacters().indexOf(mutatedAminoAcid);

        double[] mutations = Simulation.siteMutationVectors.get(mutationIndexSiteSequence)[mSiteMutationVectors][nSiteMutationVectors];

        this.correspondingGeometricPhenotype = new GeometricPhenotype(this.correspondingGeometricPhenotype.getTraitA() + mutations[0], this.correspondingGeometricPhenotype.getTraitB() + mutations[1]);

        return new GeometricSeqPhenotype(this.correspondingGeometricPhenotype.getTraitA(), this.correspondingGeometricPhenotype.getTraitB(), this.nucleotideSequence);
    }

    public String toString() {
        String fullString = String.format("%s, %.4f,%.4f", this.nucleotideSequence, this.correspondingGeometricPhenotype.getTraitA(),
                                          this.correspondingGeometricPhenotype.getTraitB());
        return fullString;
    }

    private String getRandomMutatedNucleotide(String originalCodon, int mutationIndexSite) {
        String nucleotides = "AGTC";
        int mutationIndexSiteCodon = mutationIndexSite % 3;
        char originalNucleotideToMutate = originalCodon.charAt(mutationIndexSiteCodon);

        double[] transitionTransversion = Simulation.transitionTranversionProbability.get(originalNucleotideToMutate);

        double randomNum = Math.random();
        int indexAlphabet = 0;
        for (int i = 0; i < 4;i++) {
            if (randomNum < transitionTransversion[i]) {
                indexAlphabet = i;
                break;
            }
        }
        char randomMutatedNucleotide = nucleotides.charAt(indexAlphabet);

        StringBuilder mutatedCodon = new StringBuilder(originalCodon);

        mutatedCodon.setCharAt(mutationIndexSiteCodon, randomMutatedNucleotide);
        String translatedAminoAcidFromTriplet = Simulation.codonMap.get(mutatedCodon.toString());

        int count = 0;

        while (translatedAminoAcidFromTriplet.equals("STOP")) {
            randomNum = Math.random();
            indexAlphabet = 0;
            for (int i = 0; i < 4;i++) {
                if (randomNum < transitionTransversion[i]) {
                    indexAlphabet = i;
                    break;
                }
            }
             randomMutatedNucleotide = nucleotides.charAt(indexAlphabet);

            mutatedCodon.setCharAt(mutationIndexSiteCodon, randomMutatedNucleotide);
            translatedAminoAcidFromTriplet = Simulation.codonMap.get(mutatedCodon.toString());
            count++;
            if (count > 10) {
                //System.out.println(mutationIndexSiteCodon +  " -> " + originalCodon +  " -> " + randomNum + " -> " + randomMutatedNucleotide + " -> " + translatedAminoAcidFromTriplet);
            }
        }

        StringBuilder mutatedNucleotideSequence = new StringBuilder(this.nucleotideSequence);
        mutatedNucleotideSequence.setCharAt(mutationIndexSite, randomMutatedNucleotide);
        this.nucleotideSequence = mutatedNucleotideSequence.toString();

        Simulation.mutations.println("" + originalNucleotideToMutate + randomMutatedNucleotide);

        return translatedAminoAcidFromTriplet;
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
            assert(this.nucleotideSequence.length() == Parameters.startingSequence.length()) : "Nucleotide sequences must remain the same length throughout a Simulation";
            assert(this.correspondingGeometricPhenotype != null) : "All GeometricSeqPhenotypes must have a GeometricPhenotype";
            for (int i = 0; i < this.nucleotideSequence.length(); i += 3) {
                String triplet = this.nucleotideSequence.substring(i, i+3);
                String translatedAminoAcid = Simulation.codonMap.get(triplet);

                assert(!translatedAminoAcid.equals("STOP")) : "There should not be a stop codon at site " + i;
            }
        }
    }
}
