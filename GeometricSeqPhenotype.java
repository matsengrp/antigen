public class GeometricSeqPhenotype implements Phenotype {
    public final char[] ALPHABET = Parameters.alphabet.toCharArray();
    public GeometricPhenotype correspondingGeometricPhenotype;

    private String sequence;

    public GeometricSeqPhenotype() {
        this.correspondingGeometricPhenotype = new GeometricPhenotype();
        this.sequence = startingSequenceGenerator();
    }

    public GeometricSeqPhenotype(double tA, double tB) {
        this();
        this.correspondingGeometricPhenotype = new GeometricPhenotype(tA, tB);
    }

    public GeometricSeqPhenotype(double tA, double tB, String startingSequence) {
        this(tA, tB);
        startingSequence = startingSequence.toUpperCase();
        this.sequence = startingSequence;
    }

    /*
    public GeometricSeqPhenotype(String startingSequence, Phenotype correspondingGeometricPhenotype) {
        startingSequence = startingSequence.toUpperCase();
        this.sequence = startingSequence;
        this.correspondingGeometricPhenotype = (GeometricPhenotype) correspondingGeometricPhenotype;
    }
    */

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
        return this.sequence;
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
        int indexSite = random(this.sequence.length());
        int indexAlphabet = random(this.ALPHABET.length);

        // substitute a random index of sequence with a random character in the ALPHABET
        StringBuilder mutated = new StringBuilder(this.sequence);
        mutated.setCharAt(indexSite, this.ALPHABET[indexAlphabet]);

        GeometricPhenotype mutatedCorrespondingGeometricPhenotype = (GeometricPhenotype) correspondingGeometricPhenotype.mutate();

        // TODO: GeometricSeqPhenotype(String startingSequence, Phenotype correspondingGeometricPhenotype) instead?
        return new GeometricSeqPhenotype(mutatedCorrespondingGeometricPhenotype.getTraitA(),
                                         mutatedCorrespondingGeometricPhenotype.getTraitB(), mutated.toString());
    }

    private int random(int maxRange) {
        Random random = new Random();
        return random.nextInt(0, maxRange - 1);
    }

    public String toString() {
        return this.sequence;
    }

}
