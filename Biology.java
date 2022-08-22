import cern.colt.Arrays;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Scanner;

/**
 * A class that allows antigen to model a virus's genetic sequence,
 * and how changes in sequence give rise to changes in both:
 * i) antigenic phenotype
 * ii) viral fitness in the absence of immune selection (e.g., quantified by DMS data)
 *
 * @author Thien Tran
 */
public class Biology {
    public enum DMSData {
        DMS_DATA();
        public final double[][] aminoAcidPreference;

        DMSData() {
            int numberOfAminoAcidSites = Parameters.startingSequence.length() / 3;
            int numberOfAminoAcids = Parameters.AlphabetType.AMINO_ACIDS.getValidCharacters().length();
            aminoAcidPreference = new double[numberOfAminoAcidSites][numberOfAminoAcids];

            try {
                Scanner dms = new Scanner(new File(Parameters.DMSFile));
                dms.nextLine(); // read header

                for (int i = 0; i < numberOfAminoAcidSites; i++) {
                    Scanner currentSite = new Scanner(dms.nextLine());
                    currentSite.useDelimiter(",");

                    currentSite.next(); // ignore site number

                    for (int j = 0; j < numberOfAminoAcids; j++) {
                        double probability = currentSite.nextDouble();
                        aminoAcidPreference[i][j] = probability;
                    }
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        public double[][] getAminoAcidPreference() {
            return aminoAcidPreference;
        }
    }

    /**
     * Transition/transversion ratio used to determine what a given wild type nucleotide should be mutated to.
     */
    public enum MutationType {
        MUTATION();
        public final Map<Character, double[]> transitionTranversionProbability;

        /**
         * Constructor that creates a probability distribution for each nucleotide,
         * where each possible mutation is weighted by a pre-defined transition/transversion ratio.
         */
        MutationType() {
            this.transitionTranversionProbability = new HashMap<>()  {{
                double transition = 0.5 / (Parameters.transitionTransversionRatio + 1.0);
                double transversion = Parameters.transitionTransversionRatio / (Parameters.transitionTransversionRatio + 1.0);
                // key: nucleotide
                // value: int[] that gives the probability of being < to the index which corresponds to [A, C, G, T]

                // Example:
                // Parameters.transitionTransversionRatio: 5.0
                // A: {0, 0.5/6, 5.5/6.0, 1.0}
                // C: {0.5/6.0, 0.5/6.0, 1.0/6.0, 1.0}
                // G: {5.0/6.0, 5.5/6.0, 5.5/6.0, 1.0}
                // T: {0.5/6.0, 5.5/6.0, 1.0, 1.0}

                put('A', new double[]{0, transition, transition + transversion, 1.0});
                put('C', new double[]{transition, transition, transition + transition, 1.0});
                put('G', new double[]{transversion, transversion + transition, transversion + transition, 1.0});
                put('T', new double[]{transition, transition + transversion, 1.0, 1.0});
            }};
        }

        /**
         * Returns a mutant nucleotide using the transition/transversion ratio.
         *
         * @param originalNucleotideToMutate the original nucleotide to mutate
         * @return the new nucleotide to change the original nucleotide to
         */
        public char getNucleotide(char originalNucleotideToMutate) {
            double[] transitionTransversion = this.transitionTranversionProbability.get(originalNucleotideToMutate);

            // choose a random number between 0-1
            double randomNum = Math.random();
            int indexAlphabet = 0;

            for (int i = 0; i < 4; i++) {
                if (randomNum < transitionTransversion[i]) {
                    indexAlphabet = i;
                    break;
                }
            }

            return Parameters.AlphabetType.NUCLEOTIDES.getValidCharacters().charAt(indexAlphabet);
        }
    }

    /**
     * DNA codon table
     */
    public enum CodonMap {
        CODONS();
        public final Map<String, String> codonMap;

        CodonMap() {
            this.codonMap = new HashMap<>() {{
                put("TTT", "F");
                put("TTC", "F");
                put("TTA", "L");
                put("TTG", "L");

                put("CTT", "L");
                put("CTC", "L");
                put("CTA", "L");
                put("CTG", "L");

                put("ATT", "I");
                put("ATC", "I");
                put("ATA", "I");
                put("ATG", "M");

                put("GTT", "V");
                put("GTC", "V");
                put("GTA", "V");
                put("GTG", "V");


                put("TCT", "S");
                put("TCC", "S");
                put("TCA", "S");
                put("TCG", "S");

                put("CCT", "P");
                put("CCC", "P");
                put("CCA", "P");
                put("CCG", "P");

                put("ACT", "T");
                put("ACC", "T");
                put("ACA", "T");
                put("ACG", "T");

                put("GCT", "A");
                put("GCC", "A");
                put("GCA", "A");
                put("GCG", "A");


                put("TAT", "Y");
                put("TAC", "Y");
                put("TAA", "STOP");
                put("TAG", "STOP");

                put("CAT", "H");
                put("CAC", "H");
                put("CAA", "Q");
                put("CAG", "Q");

                put("AAT", "N");
                put("AAC", "N");
                put("AAA", "K");
                put("AAG", "K");

                put("GAT", "D");
                put("GAC", "D");
                put("GAA", "E");
                put("GAG", "E");


                put("TGT", "C");
                put("TGC", "C");
                put("TGA", "STOP");
                put("TGG", "W");

                put("CGT", "R");
                put("CGC", "R");
                put("CGA", "R");
                put("CGG", "R");

                put("AGT", "S");
                put("AGC", "S");
                put("AGA", "R");
                put("AGG", "R");

                put("GGT", "G");
                put("GGC", "G");
                put("GGA", "G");
                put("GGG", "G");
            }};
        }

        /**
         * Translates and returns the amino acid specified by the given codon
         *
         * @param codon the sequence of three nucleotides to translate
         * @return the corresponding amino acid
         */
        public String getAminoAcid(String codon) {
            return this.codonMap.get(codon);
        }
    }

    /**
     * Data type used to determine how changes in sequence lead to changes in to a virus's location in antigenic space.
     * Specifically, for each site, a vector for each possible amino-acid mutations is pre-computed.
     * Vectors angles are drawn from a uniform distribution, while vector magnitudes are drawn from a gamma distribution.
     * Gamma distributions differ for epitope and non-epitope sites.
     *
     * Each site is drawn from the gamma distribution
     */
    public enum SiteMutationVectors  {
        VECTORS();

        private final int totalSites = Parameters.startingSequence.length() / 3;

        public final Map<Integer, double[][][]> siteMutationVectors;
        public final HashSet<Integer> epitopeSites;
        public final String[] stringOutputCSV = new String[totalSites]; // String[] to create CSV from

        SiteMutationVectors() {
            // SiteMutationVectors is represented using a Map that maps each protein site to a 2D array of double[].
            // Rows in the 2D array represent wild type amino acids, while columns represent mutant amino acids.
            // Each entry is a double[] of size 2 where the first index represents the x-coordinate
            // and the second index represents the y-coordinate of the corresponding vector of the wild type and mutant amino acids.
            // Example:
            //          A               C               D           . . .
            // A    [0.0, 0.0]     [-1.0, -2.0]    [-1.0, -1.0]
            // C    [1.0, 2.0]      [0.0, 0.0]      [3.0, 0.5]
            // D    [1.0, 1.0]     [-3.0, -0.5]     [0.0, 0.0]
            // .
            // .
            // .

            HashSet<Integer> epitopeSitesSet = new HashSet<>();
            for (int i = 0; i < Parameters.epitopeSites.length; i++) {
                epitopeSitesSet.add(Parameters.epitopeSites[i] - 1); // Allow users to index starting at 1, but store values starting at 0
            }
            this.epitopeSites = epitopeSitesSet;

            // create a mapping of site # to 2D array of vectors
            // mutations: 2D Array
            // i: index of wild type amino acid
            // j: index of mutant amino acid
            // ACDEFGHIKLMNPQRSTWYV 0 -> A, 1 -> C, 2 -> D, etc.
            int matrixSize = Parameters.AlphabetType.AMINO_ACIDS.getValidCharacters().length();
            Map<Integer, double[][][]> currentSiteMutationVectors  = new HashMap<>();

            // Cycle over each site in the protein sequence and create a 2D array
            for (int nucleotideSiteNumber = 0; nucleotideSiteNumber < Parameters.startingSequence.length() / 3; nucleotideSiteNumber += 1) {
                String currentStringOutputCSV = "mutation,r,theta\n";

                double[][][] currentSiteMutationMatrix = new double[matrixSize][matrixSize][];

                // Create 2D where rows are wild type mino acids and columns are mutant amino acids
                for (int wildTypeIndex = 0; wildTypeIndex < matrixSize; wildTypeIndex++) {

                    for (int mutationIndex = 0; mutationIndex < matrixSize; mutationIndex++) {
                        if (mutationIndex < wildTypeIndex) { // update lower and upper triangle in this branch.
                            // direction of mutation
                            double theta;
                            if (Parameters.mut2D) {
                                theta = Random.nextDouble(0, 2 * Math.PI);
                            } else {
                                if (Random.nextBoolean(0.5)) {
                                    theta = 0;
                                } else {
                                    theta = Math.PI;
                                }
                            }

                            // size of mutation
                            double r;
                            double alpha;
                            double beta;

                            double meanStep;
                            double sdStep;
                            if (this.epitopeSites.contains(nucleotideSiteNumber)) {
                                // epitope sites
                                meanStep = Parameters.meanStepEpitope;
                                sdStep = Parameters.sdStepEpitope;
                            } else {
                                // non-epitope sites
                                meanStep = Parameters.meanStep;
                                sdStep = Parameters.sdStep;
                            }

                            r = meanStep;
                            if (!Parameters.fixedStep) {
                                alpha = (meanStep * meanStep) / (sdStep * sdStep);
                                beta = (sdStep * sdStep) / meanStep;
                                r = Random.nextGamma(alpha, beta);
                            }

                            // create antigenic phenotype by computing the x and y coordinates (mutA and mutB) of the virus in antigenic space
                            double mutA = r * Math.cos(theta);
                            double mutB = r * Math.sin(theta);

                            // add vector=(mutA, mutB) to mutations_i,j, where
                            // mutA is the x-coordinate
                            // mutB is the y-coordinate
                            // mutations is the 2D array that corresponds to the current amino acid site
                            currentSiteMutationMatrix[wildTypeIndex][mutationIndex] = new double[]{mutA, mutB};

                            // add vector=(-mutA, -mutB) to mutations_j,i
                            currentSiteMutationMatrix[mutationIndex][wildTypeIndex] = new double[]{-1 * mutA, -1 * mutB};

                            // amino acid mutation notation
                            // wild type amino acid + site # + mutant amino acid
                            currentStringOutputCSV += "" + Parameters.AlphabetType.AMINO_ACIDS.getValidCharacters().charAt(wildTypeIndex) +
                                      nucleotideSiteNumber + Parameters.AlphabetType.AMINO_ACIDS.getValidCharacters().charAt(mutationIndex) +
                                      "," + r + "," + theta + '\n';
                        }
                    }
                }
                stringOutputCSV[nucleotideSiteNumber] = currentStringOutputCSV;
                currentSiteMutationVectors.put(nucleotideSiteNumber, currentSiteMutationMatrix);
            }
            this.siteMutationVectors = currentSiteMutationVectors;
        }

        /**
         * Retrieves the matrix associated with the (protein level) site mutation and
         * returns the vector that represents where the virus will move in Euclidean space based on the wild type and mutant amino acid.
         * Note, it's possible for wild type == mutant amino acid.
         *
         * @param mutationIndexSite the (protein level) site where a nucleotide mutated
         * @param mSiteMutationVectors the index (within Parameters.AMINO_ACIDS) that represent the wild type amino acid
         * @param nSiteMutationVectors the index (within Parameters.AMINO_ACIDS) that represent the mutant amino acid
         * @return a vector to represent the change in Euclidean space
         */
        public double[] getVector(int mutationIndexSite, int mSiteMutationVectors, int nSiteMutationVectors) {
            return this.siteMutationVectors.get(mutationIndexSite)[mSiteMutationVectors][nSiteMutationVectors];
        }

        /**
         * Returns a representation of a csv formatted as amino acid mutation notation (e.g., A4F), r, theta
         * that represents all the entries from each site's matrix of vectors.
         *
         * @return a String representing vectors in each site's matrix
         */
        public String[] getStringOutputCSV() {
            return this.stringOutputCSV;
        }

        /**
         * Return each site's matrix of vectors.
         *
         * @return site's matrices
         */
        public Map<Integer, double[][][]> getMatrices() {
            return this.siteMutationVectors;
        }

        /**
         * Returns a set of indices that represent epitope sites, which are fixed within a simulation run.
         *
         * @return epitope site indices
         */
        public HashSet getEpitopeSites() {
            return this.epitopeSites;
        }
    }
}
