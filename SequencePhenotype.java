import EDU.oswego.cs.dl.util.concurrent.FJTask;

import java.lang.*;

import static java.lang.Math.exp;
import static java.lang.Math.pow;

/**
 * <b>SequencePhenotype</b> represents a phenotype.
 * SequencePhenotypes are identified by their data contents.
 * SequencePhenotype permits all data, except for <b>""</b> and <b>null</b>.
 *
 * @author Thien Tran
 */
public class SequencePhenotype implements Phenotype {

    /**
     * The valid letters that make up the sequence of this SequencePhenotype
     */
    public final String[] NUCLEOTIDES = new String[]{"A", "C", "G", "T"};

    /**
     * The sequence of this SequencePhenotype
     */
    private String sequence;

    /**
     * Run expensive tests iff DEBUG == true.
     */
    public static final boolean DEBUG = false;

    // Abstraction Function:
    // A SequencePhenotype, s, is null if s.sequence = null, otherwise s.sequence = sequence
    // for sequence.length() > 0
    // s1.sequence and s2.sequence are allowed to be equal for two given SequencePhenotypes

    // Representation invariant for every SequencePhenotype:
    // s.sequence != null && s.sequence.length() > 0
    // for all indexSite such that s.sequence.charAt(i):
    //     indexSite == 'A' || indexSite == 'G' || indexSite == 'T' || indexSite == 'C'

    /**
     * Constructor that creates a new SequencePhenotype.
     *
     * @spec.effects Constructs a new SequencePhenotype of length Parameters.sequence.length()
     */
    public SequencePhenotype() {
        this.sequence = "";
        for (int i = 0; i < Parameters.sequence.length(); i++) {
            int indexNucleotide = random(this.NUCLEOTIDES.length);
            this.sequence += this.NUCLEOTIDES[indexNucleotide];
        }
        checkRep();
    }

    /**
     * Constructor that creates a new SequencePhenotype.
     *
     * @param sequence the sequence of the new SequencePhenotype
     * @spec.requires sequence != null && sequence.length() > 0
     * @spec.effects Constructs a new SequencePhenotype with the data content of the given parameter.
     */
    public SequencePhenotype(String sequence) {
        sequence = sequence.toUpperCase();
        this.sequence = sequence;
        checkRep();
    }

    /**
     * Returns the sequence of this SequencePhenotype
     *
     * @return the sequence of the SequencePhenotype represented by this.
     */
    public String getSequence() {
        return this.sequence;
    }

    /**
     * Sets the sequence of this SequencePhenotype to the given sequence
     *
     * @param sequence the sequence to set the sequence of this SequencePhenotype to
     * @spec.requires sequence != null && sequence.length() > 0
     * @spec.effects Changes the sequence of this SequencePhenotype
     */
    public void setSequence(String sequence) {
        this.sequence = sequence;
        checkRep();
    }

    /**
     * Compare phenotypes and report antigenic distance
     *
     * @param p the SequencePhenotype to find the hamming distance from this SequencePhenotype
     * @spec.requires p.sequence != null && p.sequence.length() == this.sequence.length() && p instanceof SequencePhenotype
     * @return the hamming distance of this SequencePhenotype and (SequencePhenotype) p
     */
    public double distance(Phenotype p) {
        // String validation
        if (!(p instanceof SequencePhenotype)) {
            throw new IllegalArgumentException("Parameter Phenotype p should be an instance of SequencePhenotype!");
        }

        SequencePhenotype seqP = (SequencePhenotype) p;
        String seq2 = seqP.getSequence();

        // equal length validation
        if (this.sequence.length() != seq2.length()) {
            System.out.println("s1 " + this.sequence);
            System.out.println("s2 " + seq2);
            throw new IllegalArgumentException("Sequence lengths are not equal!");
        }

        // Calculates the hamming distance between this.sequence and seq2
        int hammingDistance = 0;
        for (int i = 0; i < this.sequence.length(); i++) {
            if (this.sequence.charAt(i) != seq2.charAt(i)) {
                hammingDistance++;
            }
        }
        checkRep();
        return hammingDistance;
    }

    /**
     * Provides the risk of infection (from 0 to 1) of a virus with this phenotype
     * when contacting a Host with a List of Phenotypes forming their immune history
     *
     * @param history the SequencePhenotypes that have contacted a Host
     * @spec.requires history[i] instanceof SequencePhenotype for all i in history.length
     * @return the cross immunity between a virus phenotype and a host's immune history
     *         here encoded more directly as risk of infection, which ranges from 0 to 1
     */
    public double riskOfInfection(Phenotype[] history) {
        double fullImmuneRisk = 1.0;

        int maxDistance = this.sequence.length(); // The largest number of the range input.

        for(Phenotype pHistory : history) {
            double inputDistance = this.distance(pHistory); // hamming distance of this and pHistory
            double localImmuneRisk;

            switch (Parameters.crossImmunity) {
                case "linear":
                    localImmuneRisk = pow(inputDistance / maxDistance, Parameters.crossImmunityStrength);
                    break;
                case "exponential":
                    localImmuneRisk = 1 - exp(-inputDistance / Parameters.crossImmunityStrength);
                    break;
            }
            fullImmuneRisk *= localImmuneRisk;
        }

        return fullImmuneRisk;
    }

    /**
     * Returns a mutated copy of this SequencePhenotype (point substitution), original SequencePhenotype is unharmed
     *
     * @return a mutated copy of this SequencePhenotype
     */
    public SequencePhenotype mutate() {
        int indexSite = random(this.sequence.length());
        int indexNucleotide = random(this.NUCLEOTIDES.length);

        // substitute a random index of sequence with a random nucleotide
        StringBuilder mutated = new StringBuilder(this.sequence);
        mutated.setCharAt(indexSite, this.NUCLEOTIDES[indexNucleotide].charAt(0));

        checkRep();
        return new SequencePhenotype(mutated.toString());
    }

    /**
     * Returns the sequence of this SequencePhenotype (i.e. the String representation of the SequencePhenotype represented by this).
     * Valid example outputs include "AGTC" and "AAGTCGTAGTCC" and "A" and "C."
     *
     * @return the String representation of the SequencePhenotype represented by this.
     */
    public String toString() {
        return this.sequence;
    }

    private int random(int length) {
        Random random = new Random();
        return random.nextInt(0, length - 1);
    }

    /**
     * Standard equality operation.
     *
     * @param obj the object to be compared for equality
     * @return true if and only if 'obj' is an instance of a Node and 'this' and 'obj' represent the same Node.
     */
    @Override
    public boolean equals(Object obj) {
        checkRep();
        if (this == obj) {
            return true;
        }

        if (obj instanceof SequencePhenotype) {
            SequencePhenotype other = (SequencePhenotype) obj;
            checkRep();
            return this.sequence.equals(other.sequence);
        } else {
            return false;
        }
    }

    /**
     * Throws an exception if the representation invariant is violated.
     */
    private void checkRep() {
        if (DEBUG)  {
            assert (this.sequence != null) : "sequence should never be null.";
            assert (this.sequence.length() > 0) : "sequence should never be empty.";

            for (int i = 0; i < this.sequence.length(); i++) {
                String sequenceChar = ("" + this.sequence.charAt(i));
                boolean contains = java.util.Arrays.asList(this.NUCLEOTIDES).contains(sequenceChar);
                assert (contains) : sequenceChar + " is not a valid nucleotide!";
            }
        }
    }
}
