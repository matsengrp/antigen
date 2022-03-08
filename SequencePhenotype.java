import java.util.*;
import java.lang.*;

import static java.lang.Math.exp;

public class SequencePhenotype implements Phenotype {

    // constants
    public final String[] NUCLEOTIDES = new String[]{"A", "C", "G", "T"};

    // fields
    private String sequence;

    // constructor
    public SequencePhenotype() {
        int indexSite = random(NUCLEOTIDES.length);
        this.sequence = NUCLEOTIDES[indexSite];
    }
    public SequencePhenotype(String sequence) {
        sequence = sequence.toUpperCase();
        for (int i = 0; i < sequence.length(); i++) {
            String sequenceChar = ("" + sequence.charAt(i));
            boolean contains = java.util.Arrays.asList(NUCLEOTIDES).contains(sequenceChar);
            if (!contains) {
                throw new IllegalArgumentException(sequenceChar + " is not a valid nucleotide!");
            }
        }
        this.sequence = sequence;
    }

    public String getSequence() {
        return this.sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public double distance(Phenotype p) {
        // String validation
        if (!(p instanceof SequencePhenotype)) {
            throw new IllegalArgumentException("Parameter Phenotype p should be an instance of SequencePhenotype!");
        }

        SequencePhenotype seqP = (SequencePhenotype) p;

        String seq2 = seqP.getSequence();

        // equal length validation
        if (this.sequence.length() != seq2.length()) {
            throw new IllegalArgumentException("Sequence lengths are not equal!");
        }

        // Calculates the hamming distance between this.sequence and seq2
        int hammingDistance = 0;
        for (int i = 0; i < this.sequence.length(); i++) {
            if (this.sequence.charAt(i) != seq2.charAt(i)) {
                hammingDistance++;
            }
        }

        return hammingDistance;
    }

    // cross immunity between a virus phenotype and a host's immune history
    // here encoded more directly as risk of infection, which ranges from 0 to 1
    public double riskOfInfection(Phenotype[] history) {
        double fullImmuneRisk = 1.0;

        int minDistance = 0; // The lowest number of the range input.
        int maxDistance = this.sequence.length(); // The largest number of the range input.
        int minOutput = 0; // The lowest number of the range output.
        int maxOutput = 1; // The largest number of the range output.

        double slope = (double) (maxOutput - minOutput) / (maxDistance - minDistance);

        for(Phenotype pHistory : history) {
            double input = this.distance(pHistory); // hamming distance of this and pHistory
            double output = minOutput + slope * (input - minDistance);

            switch (Parameters.crossImmunity) {
                case "linear":
                    fullImmuneRisk = fullImmuneRisk * output;
                    break;
                case "exponential":
                    fullImmuneRisk = fullImmuneRisk * (1 - exp(-output));
                    break;
            }
        }

        return fullImmuneRisk;
    }

    // returns a mutated copy, original SequencePhenotype is unharmed
    public Phenotype mutate() {
        int indexSite = random(this.sequence.length());
        int indexNucleotide = random(this.NUCLEOTIDES.length);

        // substitute a random index of sequence with a random nucleotide
        StringBuilder mutated = new StringBuilder(this.sequence);
        mutated.setCharAt(indexSite, this.NUCLEOTIDES[indexNucleotide].charAt(0));
        Phenotype mutatedP = new SequencePhenotype(mutated.toString());
        return mutatedP;

    }

    public String toString() {
        return this.sequence;
    }

    private int random(int length) {
        Random random = new Random();
        return random.nextInt(0, length - 1);
    }
}