import java.util.*;
import java.lang.*;
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
            boolean contains = Arrays.asList(NUCLEOTIDES).contains(sequenceChar);
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

        // find the closest phenotype in history
        double closestDistance = 100.0;
        if (history.length > 0) {
            for (Phenotype phenotype : history) {
                double thisDistance = distance(phenotype);
                if (thisDistance < closestDistance) {
                    closestDistance = thisDistance;
                }
                if (thisDistance < 0.01) {
                    break;
                }
            }
        }

        double risk = closestDistance * Parameters.smithConversion;
        double minRisk = 1.0 - Parameters.homologousImmunity;
        risk = Math.max(minRisk, risk);
        risk = Math.min(1.0, risk);

        return risk;

    }

    // returns a mutated copy, original SequencePhenotype is unharmed
    public Phenotype mutate() {
        int indexSite = random(this.sequence.length());
        int indexNucleotide = random(this.NUCLEOTIDES.length);

        // substitute a random index of sequence with a random nucleotide
        StringBuilder mutated = new StringBuilder(this.sequence);
        mutated.setCharAt(indexSite, this.NUCLEOTIDES[indexNucleotide].charAt(0));

        return new SequencePhenotype(mutated.toString());
    }

    public String toString() {
        return this.sequence;
    }

    private int random(int length) {
        Random random = new Random();
        return random.nextInt(0, length - 1);
    }
}