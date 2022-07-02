import java.util.Arrays;

/**
 * <b>GeometricSeqPhenotype</b> represents a phenotype.
 * Distance and cross-immunity can be calculated between two GeometricSeqPhenotypes.
 * Multiple Viruses can reference a single GeometricSeqPhenotype which is identified by its data contents.
 *
 * @author Thien Tran
 */
public class GeometricSeqPhenotype extends GeometricPhenotype {
    /**
     * The valid letters that make up the nucleotide sequence of this GeometricSeqPhenotype.
     *
     * Valid letters are {A, C, G, T}
     */
    private final char[] ALPHABET = Parameters.AlphabetType.NUCLEOTIDES.getValidCharacters().toCharArray();

    /**
     * The nucleotide sequence of this GeometricSeqPhenotype
     */
    private char[] nucleotideSequence;

    /**
     * The first parameter that determines the position of this GeometricPhenotype in Euclidean space
     */
    private double traitA;

    /**
     * The second parameter that determines the position of this GeometricPhenotype in Euclidean space
     */
    private double traitB;

    /**
     * The number of epitope mutations this GeometricPhenotype went through, beginning with Parameters.startingSequence
     */
    private int eMutation = 0;

    /**
     * The number of non-epitope mutations this GeometricPhenotype went through, beginning with Parameters.startingSequence
     */
    private int neMutation = 0;

    /**
     * Run expensive tests iff DEBUG == true.
     */
    public static final boolean DEBUG = false;

    // Abstraction Function:
    // A GeometricSeqPhenotype, gsp, is null if gsp.sequence = null,
    // otherwise gsp.sequence = sequence for sequence.length() > 0 and sequence.length() % 3 == 0.
    //
    // The location of gsp in Euclidean space is defined by (traitA, traitB) each gsp is associated with a sequence of nucleotides.
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
        traitA = tA;
        traitB = tB;

        startingSequenceGenerator();
        checkRep();
    }

    // Generates a random sequence of nucleotides of length Parameters.startingSequence.length()
    private void startingSequenceGenerator() {
        for (int siteIndex = 0; siteIndex < Parameters.startingSequence.length(); siteIndex += 3) {
            String aminoAcid = "";

            // Generate codon at siteIndex
            while (aminoAcid.equals("") || aminoAcid.equals("STOP")) {
                String codon = "";
                for (int codonIndex = 0; codonIndex < 3; codonIndex++) {
                    int indexAlphabet = random(this.ALPHABET.length);
                    char nucleotide = this.ALPHABET[indexAlphabet];
                    codon += nucleotide;
                    this.nucleotideSequence[siteIndex + codonIndex] = nucleotide;
                }
                aminoAcid = Biology.CodonMap.CODONS.getAminoAcid(codon);
            }
        }
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
        traitA = tA;
        traitB = tB;
        this.nucleotideSequence = startingSequence;
        checkRep();
    }

    public GeometricSeqPhenotype(double tA, double tB, char[] startingSequence, int e, int ne) {
        traitA = tA;
        traitB = tB;
        this.nucleotideSequence = startingSequence;
        this.eMutation = e;
        this.neMutation = ne;
        checkRep();
    }

    /**
     * Return x-coordinate of where this GeometricSeqPhenotype is in Euclidean space
     *
     * @return position of this GeometricSeqPhenotype along the x-axis
    */
    public double getTraitA() {
        return traitA;
    }

    /**
     * Return y-coordinate of where this GeometricSeqPhenotype is in Euclidean space
     *
     * @return position of this GeometricSeqPhenotype along the y-axis
     */
    public double getTraitB() {
        return traitB;
    }

    public String getSequence() {
        return String.valueOf(this.nucleotideSequence);
    }

    /**
     * Return the (Euclidean) distance between this GeometricSeqPhenotype and p in Euclidean space
     *
     * @param p  GeometricSeqPhenotype to calculate the distance between
     * @return distance between this GeometricSeqPhenotype and p
     */
    public double distance(Phenotype p) {
        GeometricSeqPhenotype p2d = (GeometricSeqPhenotype) p;
        double distA = (getTraitA() - p2d.getTraitA());
        double distB = (getTraitB() - p2d.getTraitB());
        double dist = (distA * distA) + (distB * distB);
        dist = Math.sqrt(dist);
        checkRep();
        return dist;
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
        //  - String since Strings are immutable. New Strings have to be created to mutate a sequence,
        //    which is messy and uses up extra heap space. Plus, Strings include an extra byte '\n'.
        //  - StringBuffer because it's just a Wrapper around a char[], and has functionality that we don't need such as resizing.

        String wildTypeAminoAcid = "", mutantAminoAcid = "";
        int mutationIndexSiteSequence = -1;
        int mutationIndexSite = -1;
        char mutatedNucleotide = '@', originalNucleotide = '@';

        // mutate at least once
        // keep mutating if mutation generates a stop codon
        while (mutantAminoAcid.equals("") || mutantAminoAcid.equals("STOP")) {
            // choose random index to mutate in this.nucleotideSequence
            mutationIndexSiteSequence = random(this.nucleotideSequence.length);

            originalNucleotide = this.nucleotideSequence[mutationIndexSiteSequence];

            int mutationIndexSiteCodon = mutationIndexSiteSequence % 3; // index of nucleotide being mutated in codon {0, 1, 2}
            mutationIndexSite = mutationIndexSiteSequence / 3; // site # where mutation is occurring {0, . . ., total number of sites - 1}
            int codonStartIndex = 3 * mutationIndexSite; // start index of nucleotide being mutated in codon {0, 3, 6, . . .}

            char originalNucleotideToMutate = this.nucleotideSequence[mutationIndexSiteSequence];
            String originalCodon = "" + this.nucleotideSequence[codonStartIndex] +
                                        this.nucleotideSequence[codonStartIndex + 1] +
                                        this.nucleotideSequence[codonStartIndex + 2];
            wildTypeAminoAcid = Biology.CodonMap.CODONS.getAminoAcid(originalCodon);

            // get mutant nucleotide (transition versus transversion)
            mutatedNucleotide = Biology.MutationType.MUTATION.getNucleotide(originalNucleotideToMutate);

            // get new codon after mutation occurs
            StringBuffer mutatedCodon = new StringBuffer(originalCodon);
            mutatedCodon.setCharAt(mutationIndexSiteCodon, mutatedNucleotide);
            mutantAminoAcid = Biology.CodonMap.CODONS.getAminoAcid(mutatedCodon.toString());
        }

        // get new nucleotide sequence after mutating at the random index
        char[] copyNucleotideSequence = Arrays.copyOf(this.nucleotideSequence, Parameters.startingSequence.length());
        copyNucleotideSequence[mutationIndexSiteSequence] = mutatedNucleotide;

        // get indices for the matrix based on the wild type and mutant amino acids
        // matrix i,j correspond with the String "ACDEFGHIKLMNPQRSTWYV"
        int mSiteMutationVectors = Parameters.AlphabetType.AMINO_ACIDS.getValidCharacters().indexOf(wildTypeAminoAcid);
        int nSiteMutationVectors = Parameters.AlphabetType.AMINO_ACIDS.getValidCharacters().indexOf(mutantAminoAcid);

        double[] mutations = Biology.SiteMutationVectors.VECTORS.getVector(mutationIndexSite, mSiteMutationVectors, nSiteMutationVectors); // use matrix at site # where mutation is occurring
        double mutA = getTraitA() + mutations[0]; // r * cos(theta)
        double mutB = getTraitB() + mutations[1]; // r * sin(theta)

        // count number of epitope and non-epitope mutations
        int eMutationNew = this.eMutation;
        int neMutationNew = this.neMutation;
        if (Biology.SiteMutationVectors.VECTORS.getEpitopeSites().contains(mutationIndexSite)) {
            eMutationNew += 1;
        } else {
            neMutationNew += 1;
        }

        // For testing transition versus transversion ratio
        // TestGeometricSeqPhenotype.mutations.println("" + originalNucleotide + mutatedNucleotide);

        checkRep();
        Phenotype mutP = new GeometricSeqPhenotype(mutA,mutB, copyNucleotideSequence, eMutationNew, neMutationNew);
        return mutP;
    }

    /**
     * Returns the sequence, the position of this GeometricSeqPhenotype, and number of epitope and non-epitope mutations
     * (i.e. the String representation of the GeometricSeqPhenotype represented by this).
     * Valid example outputs include "ACG, 0.0, 0.0, 0, 0" and "ACGTGTACGTGT, 2.3, 8.9, 40, 10".
     *
     * @return the String representation of the GeometricSeqPhenotype represented by this.
     */
    public String toString() {
        String fullString = String.format("%s, %.4f, %.4f, %d, %d", String.valueOf(this.nucleotideSequence), this.getTraitA(),
                                          this.getTraitB(), this.eMutation, this.neMutation);
        return fullString;
    }

    // Returns a random index within the bounds of String.length.()=maxRange
    private int random(int maxRange) {
        Random random = new Random();
        return random.nextInt(0, maxRange - 1);
    }

    // Throws an exception if the representation invariant is violated.
    private void checkRep() {
        if (DEBUG)  {
            assert(this.nucleotideSequence.length == Parameters.startingSequence.length()) : "Nucleotide sequences must remain the same length throughout a Simulation";
            for (int i = 0; i < this.nucleotideSequence.length; i += 3) {
                String triplet = "" + this.nucleotideSequence[i] + this.nucleotideSequence[i + 1] + this.nucleotideSequence[i + 2];
                String translatedAminoAcid = Biology.CodonMap.CODONS.getAminoAcid(triplet);

                assert(!translatedAminoAcid.equals("STOP")) : "There should not be a stop codon at site " + (i/3);
            }
        }
    }
}
