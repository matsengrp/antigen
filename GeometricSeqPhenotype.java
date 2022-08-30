import java.util.Arrays;

/**
 * <b>GeometricSeqPhenotype</b> stores a Virus's genetic sequence and antigenic phenotype,
 * and tracks the cumulative number of mutations in epitope and non-epitope sites.
 * Antigenic phenotype is defined as a position in 2D space, where the antigenic distance
 * between two viruses can be computed as the Euclidean distance between them.
 * Multiple Viruses can reference a single GeometricSeqPhenotype which is identified by its data contents.
 *
 * @author Thien Tran
 */
public class GeometricSeqPhenotype extends GeometricPhenotype {
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
     * The number of epitope mutations this GeometricPhenotype went through (counting from the startingSequence)
     */
    private int eMutation;

    /**
     * The number of non-epitope mutations this GeometricPhenotype went through (counting from the startingSequence)
     */
    private int neMutation;

    /**
     *
     */
    public static final boolean SANITY_TEST = true;

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
        this.traitA = 0.0;
        this.traitB = 0.0;
        this.nucleotideSequence = Parameters.startingSequence.toCharArray();
        this.eMutation = 0;
        this.neMutation = 0;
        checkRep();
    }

    /**
     * Constructor that creates a new GeometricSeqPhenotype.
     *
     * @param tA the x-coordinate of the new GeometricSeqPhenotype.
     * @param tB the y-coordinate of the new GeometricSeqPhenotype.
     * @effects Constructs a new GeometricSeqPhenotype that has a random starting sequence of length Parameters.sequence.length()
     * and GeometricPhenotype with the data content of the given parameters.
     */
    public GeometricSeqPhenotype(double tA, double tB) {
        this();
        this.traitA = tA;
        this.traitB = tB;
        checkRep();
    }

    /**
     * Constructor that creates a new GeometricSeqPhenotype.
     *
     * @param tA the x-coordinate of the new GeometricSeqPhenotype.
     * @param tB the y-coordinate of the new GeometricSeqPhenotype.
     * @param startingSequence the starting nucleotide sequence of the GeometricSeqPhenotype to be constructed.
     * @requires nucleotideSequence != null && nucleotideSequence.length() > 0 && nucleotideSequence.length() % 3 == 0
     * @effects Constructs a new GeometricSeqPhenotype that has a starting sequence
     * and GeometricPhenotype with the data content of the given parameters.
     */
    public GeometricSeqPhenotype(double tA, double tB, char[] startingSequence) {
        this();
        this.traitA = tA;
        this.traitB = tB;
        this.nucleotideSequence = startingSequence;
        checkRep();
    }

    /**
     * Constructor that creates a new GeometricSeqPhenotype.
     *
     * @param tA the x-coordinate of the new GeometricSeqPhenotype.
     * @param tB the y-coordinate of the new GeometricSeqPhenotype.
     * @param startingSequence the starting nucleotide sequence of the GeometricSeqPhenotype to be constructed.
     * @param e the number of epitope mutations this GeometricSeqPhenotype is from startingSequence
     * @param nE the number of non-epitope mutations from startingSequence is from startingSequence
     * @requires nucleotideSequence != null && nucleotideSequence.length() > 0 && nucleotideSequence.length() % 3 == 0
     * @effects Constructs a new GeometricSeqPhenotype that has a starting sequence
     * and GeometricPhenotype with the data content of the given parameters.
     */
    public GeometricSeqPhenotype(double tA, double tB, char[] startingSequence, int e, int nE) {
        this.traitA = tA;
        this.traitB = tB;
        this.nucleotideSequence = startingSequence;
        this.eMutation = e;
        this.neMutation = nE;
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

    /**
     * Returns the nucleotide sequence of this GeometricSeqPhenotype.
     * Valid example outputs include "ACG" and "ACGTGTACGTGT"
     *
     * @return the nucleotide sequence of the GeometricSeqPhenotype represented by this.
     */
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

        // this.nucleotideSequence is represented using a char[] instead of
        //  - String since Strings are immutable. New Strings have to be created to mutate a sequence,
        //    which is messy and uses up extra heap space. Plus, Strings include an extra byte '\n'.
        //  - StringBuffer because it's just a Wrapper around a char[], and has functionality that we don't need such as resizing.

        String wildTypeAminoAcid = "", mutantAminoAcid = "";
        int nucleotideMutationIndex = -1;
        char wildTypeNucleotide = ' ', mutantNucleotide = ' ';
        int cycle = 0;

        // Make a single nucleotide mutation to the sequence
        // If the mutation results in a stop codon, then throw that mutation away and try another one
        while (mutantAminoAcid.equals("") || mutantAminoAcid.equals("STOP")) {
            cycle++;
            // choose random index to mutate in this.nucleotideSequence
            nucleotideMutationIndex = Random.nextInt(0, this.nucleotideSequence.length - 1);

            wildTypeNucleotide = this.nucleotideSequence[nucleotideMutationIndex];

            // get mutant nucleotide (transition/transversion ratio)
            mutantNucleotide = Biology.MutationType.MUTATION.getNucleotide(wildTypeNucleotide);

            String[] wildTypeMutantAminoAcids = mutateHelper(nucleotideMutationIndex, mutantNucleotide);

            wildTypeAminoAcid = wildTypeMutantAminoAcids[0];
            mutantAminoAcid = wildTypeMutantAminoAcids[1];

            if (SANITY_TEST) {
                int nucleotideMutationCodonIndex = nucleotideMutationIndex % 3;
                int proteinMutationIndex = nucleotideMutationIndex / 3;
                int nucleotideMutationFirstCodonIndex = 3 * proteinMutationIndex;
                String wildTypeCodon = "" + this.nucleotideSequence[nucleotideMutationFirstCodonIndex] +
                        this.nucleotideSequence[nucleotideMutationFirstCodonIndex + 1] +
                        this.nucleotideSequence[nucleotideMutationFirstCodonIndex + 2];

                StringBuilder mutantCodon = new StringBuilder(wildTypeCodon);
                mutantCodon.setCharAt(nucleotideMutationCodonIndex, mutantNucleotide);

                // (proteinMutationIndex + 1) to show one-based numbering
                TestGeometricSeqPhenotype.mutations.println((nucleotideMutationIndex + 1) + "," + wildTypeNucleotide + mutantNucleotide + "," +
                        wildTypeAminoAcid + "," + mutantAminoAcid + "," + wildTypeCodon + "," + mutantCodon + "," +  cycle + "," + this.hashCode());
            }
        }

        int proteinMutationIndex = nucleotideMutationIndex / 3; // site # where mutation is occurring {0, . . ., total number of sites - 1}
        // Update the x and y coordinates of the virus in antigenic space
        double[] vectors = updateAntigenicPhenotype(proteinMutationIndex, wildTypeAminoAcid, mutantAminoAcid);
        double mutA = vectors[0];
        double mutB = vectors[1];

        // Update the virus's nucleotide sequence by introducing the mutation from above
        char[] copyNucleotideSequence = Arrays.copyOf(this.nucleotideSequence, Parameters.startingSequence.length());
        copyNucleotideSequence[nucleotideMutationIndex] = mutantNucleotide;

        // Determine whether the mutation occurred in an epitope or non-epitope site, and update
        // the counts of epitope and non-epitope mutations accordingly
        int eMutationNew = this.eMutation;
        int neMutationNew = this.neMutation;
        if (Biology.SiteMutationVectors.VECTORS.getEpitopeSites().contains(proteinMutationIndex)) {
            eMutationNew += 1;
        } else {
            neMutationNew += 1;
        }

        checkRep();
        Phenotype mutP = new GeometricSeqPhenotype(mutA, mutB, copyNucleotideSequence, eMutationNew, neMutationNew);
        return mutP;
    }

    // Updates a virus's location in antigenic space upon a mutation by taking the
    // vector giving the virus's current location and then summing it with a
    // precomputed vector that gives the antigenic effect of the mutation
    private double[] updateAntigenicPhenotype(int mutationIndexSite, String wildTypeAminoAcid, String mutantAminoAcid) {
        double mutA = 0.0; // virus's location in the x dimension
        double mutB = 0.0; // virus's location in the y dimension
        if (!wildTypeAminoAcid.equals(mutantAminoAcid)) {
            // get indices for the matrix based on the wild type and mutant amino acids
            // matrix i,j correspond with the String "ACDEFGHIKLMNPQRSTWYV"
            int mSiteMutationVectors = Parameters.AlphabetType.AMINO_ACIDS.getValidCharacters().indexOf(wildTypeAminoAcid);
            int nSiteMutationVectors = Parameters.AlphabetType.AMINO_ACIDS.getValidCharacters().indexOf(mutantAminoAcid);

            double[] mutations = Biology.SiteMutationVectors.VECTORS.getVector(mutationIndexSite, mSiteMutationVectors, nSiteMutationVectors); // use matrix at site # where mutation is occurring
            mutA = getTraitA() + mutations[0]; // r * cos(theta)
            mutB = getTraitB() + mutations[1]; // r * sin(theta)
        }

        return new double[]{mutA, mutB};
    }

    // Mutates nucleotide sequence at given nucleotideMutationIndex with given char mutantNucleotide.
    // Returns the wild type and mutant amino acid as a String[]
    //
    // package private helper method so that the updating of the nucleotide sequence can be tested in TestGeometricSeqPhenotype.java
    String[] mutateHelper(int nucleotideMutationIndex, char mutantNucleotide) {
        int nucleotideMutationCodonIndex = nucleotideMutationIndex % 3; // index of nucleotide being mutated in codon {0, 1, 2}
        int proteinMutationIndex = nucleotideMutationIndex / 3; // protein site # where mutation is occurring {0, . . ., total number of sites - 1}
        int nucleotideMutationFirstCodonIndex = 3 * proteinMutationIndex; // start index of nucleotide being mutated in codon {0, 3, 6, . . .}

        String wildTypeCodon = "" + this.nucleotideSequence[nucleotideMutationFirstCodonIndex] +
                this.nucleotideSequence[nucleotideMutationFirstCodonIndex + 1] +
                this.nucleotideSequence[nucleotideMutationFirstCodonIndex + 2];
        String wildTypeAminoAcid = Biology.CodonMap.CODONS.getAminoAcid(wildTypeCodon);

        // get new codon after mutation occurs
        StringBuilder mutantCodon = new StringBuilder(wildTypeCodon);
        mutantCodon.setCharAt(nucleotideMutationCodonIndex, mutantNucleotide);
        String mutantAminoAcid = Biology.CodonMap.CODONS.getAminoAcid(mutantCodon.toString());

        return new String[]{wildTypeAminoAcid, mutantAminoAcid};
    }

    /**
     * Returns the virus's sequence, position in antigenic space,
     * and cumulative number of epitope and non-epitope mutations in its evolutionary history
     * Valid example outputs include "ACG, 0.0, 0.0, 0, 0" and "ACGTGTACGTGT, 2.3, 8.9, 40, 10".
     *
     * @return the String representation of the GeometricSeqPhenotype represented by this.
     */
    public String toString() {
        String fullString = String.format("%s, %.4f, %.4f, %d, %d", String.valueOf(this.nucleotideSequence), this.getTraitA(),
                                          this.getTraitB(), this.eMutation, this.neMutation);
        return fullString;
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
