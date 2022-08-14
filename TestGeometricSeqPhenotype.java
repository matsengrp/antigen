import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Map;
import java.util.Random;
import java.util.Scanner;

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
        // initialize static parameters
        Parameters.load();
        Parameters.initialize();

        Random r = new Random();
        long seed = 5;
        r.setSeed(seed);
        emptyPheno = new GeometricSeqPhenotype();
        simplePheno = new GeometricSeqPhenotype(0.0, 0.0, new char[]{'T', 'G', 'C', 'A', 'T', 'C'});

        simplePheno2 = new GeometricSeqPhenotype(20.0, 10.0, new char[]{'T', 'G', 'C', 'A', 'T', 'C'});

        // Create a small history of GeometricSeq phenotypes
        history = new GeometricSeqPhenotype[3];

        history[0] = new GeometricSeqPhenotype(1.0, 2.0, new char[]{'T', 'G', 'C', 'A', 'T', 'C'});
        history[1] = new GeometricSeqPhenotype(3.0, 2.0, new char[]{'T', 'G', 'C', 'A', 'T', 'C'});
        history[2] = new GeometricSeqPhenotype(5.0, 8.0, new char[]{'T', 'G', 'C', 'A', 'T', 'C'});
    }

    /**
     * Test constructors for the SequencePhenotype class.
     */
    @Test
    public void testConstructors() {
        assertEquals(0, emptyPheno.getSequence().length() % 3);
        assertEquals(Parameters.startingSequence.length(), emptyPheno.getSequence().length());

        assertEquals("TGCATC", simplePheno.getSequence());
    }

    /**
     * Test getter method.
     */
    @Test
    public void testGetTraits() {
        assertEquals(0.0, simplePheno.getTraitA(), 0);
        assertEquals(0.0, simplePheno.getTraitB(), 0);
    }

    /**
     * Test getter method.
     */
    @Test
    public void testGetSequence() {
        assertEquals("TGCATC", simplePheno.getSequence());
    }

    @Test
    public void testDistance() {
        assertEquals(2.236068, simplePheno.distance(history[0]), 0.000001);
        assertEquals(3.605551, simplePheno.distance(history[1]), 0.000001);
        assertEquals(9.433981, simplePheno.distance(history[2]), 0.000001);

        assertEquals(7.211103, history[2].distance(history[0]), 0.000001);
        assertEquals(7.211103, history[0].distance(history[2]), 0.000001);
    }

    /**
     * Test mutate function.
     */
    @Test
    public void testMutate() {
        GeometricSeqPhenotype mutantPheno = (GeometricSeqPhenotype) simplePheno.mutate();

        assertEquals("TGCATC", simplePheno.getSequence());
        assertNotEquals("TGCATC", mutantPheno.getSequence());

        String[] wildTypeMutantAminoAcids = simplePheno.mutateHelper(0, 'A');
        assertEquals(new String[]{"C", "S"}, wildTypeMutantAminoAcids);
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
        assertNotEquals("TGCATC, 0.0000, 0.0000, 0, 0", simplePheno.mutate().toString());
    }

    /**
     * Test that DMS data stored in Antigen reflects the data passed in.
     */
    @Test
    public void testDMSData() throws FileNotFoundException {
        int numberOfAminoAcidSites = Parameters.startingSequence.length() / 3;
        int numberOfAminoAcids = Parameters.AlphabetType.AMINO_ACIDS.getValidCharacters().length();

        int numberOfSites = 0;

        for (int i = 0; i < numberOfAminoAcidSites; i++) {
            numberOfSites++;
            double currentSiteSum = 0.0;
            for (int j = 0; j < numberOfAminoAcids; j++) {
                currentSiteSum += Biology.DMSData.DMS_DATA.getAminoAcidPreference()[i][j];
            }
            assertEquals(1.0, currentSiteSum, 0.0001);
        }

        Scanner DMSData = new Scanner(new File(Parameters.DMSData));

        int dmsDataLineCount = 0;
        DMSData.nextLine(); // ignore the header
        while (DMSData.hasNextLine()) {
            dmsDataLineCount++;
            DMSData.nextLine();
        }

        assertEquals(dmsDataLineCount, numberOfSites);
    }

    /**
     * Creates a csv file for each amino acid site's matrix of vectors in a directory, test/valuesGammaDistribution.
     */
    @Test
    public void testGammaDistribution() throws IOException {
        new File("testGeometricSeqPhenotype/valuesGammaDistribution").mkdirs();

        String[] values = Biology.SiteMutationVectors.VECTORS.getStringOutputCSV();

        for (int i = 0; i < values.length; i++) {
            String value = values[i];
            PrintStream output = new PrintStream("testGeometricSeqPhenotype/valuesGammaDistribution/0_site" + i + ".csv");
            output.println(value);

            Map<Integer, double[][][]> matrices = Biology.SiteMutationVectors.VECTORS.getMatrices();
        }

        // run python testGammaDistribution.py for each i such that "0_site" + i + ".csv"
    }

    /**
     * PrintStream (mutations.csv) to print wild type and mutant nucleotide pairs to.
     */
    public static PrintStream mutations;

    static {
        try {
            mutations = new PrintStream("testGeometricSeqPhenotype/mutations.csv");
            mutations.println(Parameters.startingSequence);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // run python testTransitionTransversion.py
    }
}
