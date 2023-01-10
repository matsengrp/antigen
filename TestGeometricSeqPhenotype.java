import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

import static org.junit.Assert.*;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

/**
 * A class for testing the methods in GeometricSeqPhenotype.java
 *
 * @author Thien Tran
 */
public class TestGeometricSeqPhenotype {
    // Private fields
    private GeometricSeqPhenotype emptyPheno; // Default constructor
    private GeometricSeqPhenotype simplePheno; // Give a specific sequence
    private GeometricSeqPhenotype simplePheno2; // Give a specific sequence

    GeometricSeqPhenotype [] history;

    /**
     * Define variables needed for multiple tests.
     */
    @Before
    public void setUp() {
        // Use default static parameters from Parameters.
        // Update startingSequence and homologousImmunity only.
        Parameters.startingSequence = "TGCATC";
        Parameters.homologousImmunity = 0.95;

        emptyPheno = new GeometricSeqPhenotype();
        simplePheno = new GeometricSeqPhenotype(0.0, 0.0, new char[]{'T', 'G', 'C', 'A', 'T', 'C'});

        simplePheno2 = new GeometricSeqPhenotype(20.0, 10.0, new char[]{'G', 'C', 'T', 'G', 'G', 'A'});

        // Create a small history of GeometricSeq phenotypes.
        history = new GeometricSeqPhenotype[3];

        history[0] = new GeometricSeqPhenotype(1.0, 2.0, new char[]{'T', 'G', 'C', 'G', 'C', 'C'});
        history[1] = new GeometricSeqPhenotype(3.0, 2.0, new char[]{'T', 'G', 'C', 'G', 'C', 'T'});
        history[2] = new GeometricSeqPhenotype(5.0, 8.0, new char[]{'C', 'T', 'C', 'G', 'C', 'T'});
    }

    /**
     * Test constructors for the SequencePhenotype class.
     */
    @Test
    public void testConstructors() {
        // Check that empty constructor works and data contents are valid.
        // Nucleotide sequence must be a multiple of 3.
        assertEquals(0, emptyPheno.getSequence().length() % 3);
        // Nucleotide sequence must be the same length as the startingSequence.
        assertEquals(Parameters.startingSequence.length(), emptyPheno.getSequence().length());

        // Check that non-empty constructor works and data contents are equal to the given parameters.
        assertEquals("GCTGGA", simplePheno2.getSequence());
        assertEquals(20.0, simplePheno2.getTraitA(), 0.0);
        assertEquals(10.0, simplePheno2.getTraitB(), 0.0);
    }

    /**
     * Test getter method.
     */
    @Test
    public void testGetTraits() {
        assertEquals(0.0, simplePheno.getTraitA(), 0);
        assertEquals(0.0, simplePheno.getTraitB(), 0);

        assertEquals(20.0, simplePheno2.getTraitA(), 0);
        assertEquals(10.0, simplePheno2.getTraitB(), 0);
    }

    /**
     * Test getter method.
     */
    @Test
    public void testGetSequence() {
        assertEquals("TGCATC", simplePheno.getSequence());
    }

    /**
     * Test distance method.
     */
    @Test
    public void testDistance() {
        // Calculates the distance between this GeometricSeqPhenotype and p in Euclidean space.

        // simplePheno: traitA=0.0, traitB=0.0
        // history[0]: traitA=1.0, traitB=2.0
        // double distA = (0.0 - 1.0)
        // double distB = (0.0 - 2.0)
        // double dist = (-1.0 * -1.0) + (-2.0 * -2.0) = 5.0
        // dist = Math.sqrt(5.0) = 2.2360679775
        assertEquals(2.236068, simplePheno.distance(history[0]), 0.000001);
        assertEquals(3.605551, simplePheno.distance(history[1]), 0.000001);
        assertEquals(9.433981, simplePheno.distance(history[2]), 0.000001);

        // history[0]: traitA=1.0, traitB=2.0
        // history[2]: traitA=5.0, traitB=8.0
        // double distA = (1.0 - 5.0)
        // double distB = (2.0 - 8.0)
        // double dist = (-4.0 * -4.0) + (-6.0 * -6.0) = 52.0
        // dist = Math.sqrt(52.0) = 2.7.21110255093
        assertEquals(7.211103, history[2].distance(history[0]), 0.000001);
        assertEquals(7.211103, history[0].distance(history[2]), 0.000001);
    }

    /**
     * Test mutate function.
     */
    @Test
    public void testMutate() {
        GeometricSeqPhenotype mutantPheno = (GeometricSeqPhenotype) simplePheno.mutate();

        String simplePhenoSeq = simplePheno.getSequence();
        String mutantPhenoSeq = mutantPheno.getSequence();

        // Make sure that the Hamming distance between the mutated and un-mutated nucleotide sequences equals one.
        assertEquals(simplePhenoSeq.length(), mutantPhenoSeq.length());
        int counter = 0;
        for (int i = 0; i < simplePhenoSeq.length(); i++) {
            if (simplePhenoSeq.charAt(i) != mutantPhenoSeq.charAt(i)) {
                counter++;
            }
        }
        assertEquals(1, counter);

        assertEquals("TGCATC", simplePhenoSeq);
        assertNotEquals("TGCATC", mutantPhenoSeq);

        // Mutation in first codon and changes protein sequence
        String[] simpleWildTypeMutantAminoAcids1 = simplePheno.mutateHelper(0, 'A');
        assertArrayEquals(new String[]{"C", "S"}, simpleWildTypeMutantAminoAcids1);

        // Mutation in first codon and doesn't change protein sequence
        String[] simple2WildTypeMutantAminoAcids1 = simplePheno2.mutateHelper(1, 'C');
        assertArrayEquals(new String[]{"A", "A"}, simple2WildTypeMutantAminoAcids1);

        // Mutation in second codon and changes protein sequence
        String[] simple2WildTypeMutantAminoAcids2 = simplePheno2.mutateHelper(3, 'A');
        assertArrayEquals(new String[]{"G", "R"}, simple2WildTypeMutantAminoAcids2);

        // Mutation in second codon and doesn't change protein sequence
        String[] simpleWildTypeMutantAminoAcids2 = simplePheno.mutateHelper(5, 'T');
        assertArrayEquals(new String[]{"I", "I"}, simpleWildTypeMutantAminoAcids2);
    }

    /**
     * Test riskOfInfection calculations.
     */
    @Test
    public void testRiskOfInfection() {
        // closestDistance * Parameters.smithConversion = 2.236068 * 0.1 = 0.2236068
        // minRisk = 1.0 - Parameters.homologousImmunity = 1.0 - 0.95 = 0.05
        // risk = Math.max(minRisk, risk) = max(0.05, 0.2236068) = 0.2236068
        // risk = Math.min(1.0, risk) = 0.2236068
        assertEquals(0.2236068, simplePheno.riskOfInfection(history), 0.00001);

        // closestDistance * Parameters.smithConversion = 15.132746 * 0.1 = 1.5132746
        // minRisk = 1.0 - Parameters.homologousImmunity = 1.0 - 0.95 = 0.05
        // risk = Math.max(minRisk, risk) = max(0.05, 1.5132746) = 1.5132746
        // risk = Math.min(1.0, risk) = 1.0
        assertEquals(1.0, simplePheno2.riskOfInfection(history), 0.0);
    }

    /**
     * A redundant test for returning the sequence, but with toString().
     */
    @Test
    public void testToString() {
        // This may change as the class develops (i.e., new enhancements).
        assertEquals("TGCATC, 0.0000, 0.0000, 0, 0", simplePheno.toString());

        // Check that String representation changes after mutation.
        assertNotEquals("TGCATC, 0.0000, 0.0000, 0, 0", simplePheno.mutate().toString());
        // traitA and traitB are not deterministic.
        // E should be 0 and nE should be 1
        assertTrue(simplePheno.mutate().toString().contains(", 0, 1"));
    }

    /**
     * Test that all 64 codons are accounted for and corresponds to the correct amino acid.
     */
    @Test
    public void testCodonMap() throws FileNotFoundException {
        // Initialize static parameters to create CodonMap in Biology
        Parameters.load();
        Parameters.initialize();

        // Create a new CodonMap from a codon table txt file,
        // which will be used to validate CodonMap in Biology
        Map<String, String> codonMapTest = new HashMap<>();

        Scanner codonTable = new Scanner(new File("codon_table.txt"));
        codonTable.nextLine();

        while (codonTable.hasNextLine()) {
            Scanner codonLine = new Scanner(codonTable.nextLine());

            String codon = codonLine.next(); // codon
            codonLine.next(); // throw out second column (amino acid abbreviation)
            String aminoAcid = codonLine.next(); // amino acid codon codes for

            // txt file saved from GitHub uses "O" instead of "STOP" to represent stop codons
            if (aminoAcid.equals("O")) {
                aminoAcid = "STOP";
            }

            codonMapTest.put(codon, aminoAcid);
        }

        assertEquals(codonMapTest, Biology.CodonMap.CODONS.codonMap);
    }

    /**
     * Test that DMS data stored in Antigen reflects the csv file read from parameters.yml.
     */
    @Test
    public void testDMSData() {
        //  The number of rows in the DMS file is checked in Parameters.

        // Check if each row (array) in the DMS data sums up to 1.0.
        int numberOfAminoAcidSites = Parameters.startingSequence.length() / 3;
        int numberOfAminoAcids = Biology.AlphabetType.AMINO_ACIDS.getValidCharacters().length();

        // Cycle through each array, and sum up each element in the nexted array.
        for (int i = 0; i < numberOfAminoAcidSites; i++) {
            double currentSiteSum = 0.0;
            for (int j = 0; j < numberOfAminoAcids; j++) {
                currentSiteSum += Biology.DMSData.DMS_DATA.getAminoAcidPreference(i)[j];
            }
            assertEquals(1.0, currentSiteSum, 0.0001);
        }
    }

    /**
     * Creates a csv file for each amino acid site's matrix of vectors in a directory, test/valuesGammaDistribution.
     */
    @Test
    public void testGammaDistribution() throws IOException {
        // Initialize static parameters to test different distributions
        Parameters.load();
        Parameters.initialize();
        new File("testGeometricSeqPhenotype/valuesGammaDistribution").mkdirs();

        String[] values = Biology.SiteMutationVectors.VECTORS.getStringOutputCSV();
        Map<Integer, double[][][]> matrices = Biology.SiteMutationVectors.VECTORS.getMatrices();

        for (int i = 0; i < values.length; i++) {
            String value = values[i];
            PrintStream output = new PrintStream("testGeometricSeqPhenotype/valuesGammaDistribution/0_site" + i + ".csv");
            output.println(value);

            double[][][] matrix = matrices.get(i);
            System.out.println("Matrix " + i);
            for(int j = 0; j < Biology.AlphabetType.AMINO_ACIDS.getValidCharacters().length(); j++) {
                for(int k = 0; k < Biology.AlphabetType.AMINO_ACIDS.getValidCharacters().length(); k++) {
                    // Nulls are along the diagonal
                    System.out.print(Arrays.toString(matrix[j][k]));
                }
                System.out.println();
            }
        }

        // Run: python testGammaDistribution.py for each i such that "0_site" + i + ".csv"
    }

    /**
     * PrintStream (mutations.csv) to print wild type and mutant nucleotide pairs to.
     */
    public static PrintStream mutations;

    static {
        try {
            if (GeometricSeqPhenotype.SANITY_TEST) {
                mutations = new PrintStream("testGeometricSeqPhenotype/mutations.csv");
                mutations.println("siteN,wildCodon,mutantCodon,pairWildMutantN,wildAA,mutantAA,cycle,id");
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // Run: python testTransitionTransversion.py
    }
}
