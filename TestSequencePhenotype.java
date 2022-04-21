import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import java.util.Random;

/**
 * @author zorian15
 * @date 2.14.2022
 *
 * A class for testing the methods in SequencePhenotype.java
 */
public class TestSequencePhenotype {

	// Private fields
	private SequencePhenotype emptyPheno; // Default constructor
	private SequencePhenotype simplePheno; // Give a specific sequence
	

	/**
	 * Define variables needed for multiple tests.
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {
		Random r = new Random();
		long seed = 5;
		r.setSeed(seed);
		emptyPheno = new SequencePhenotype();

		// Create a small history of Sequence phenotypes
		SequencePhenotype [] history = new SequencePhenotype[3];

		switch (Parameters.alphabet) {
			case "nucleotides":
				simplePheno = new SequencePhenotype("ACGT");

				history[0] = new SequencePhenotype("ATGT");
				history[1] = new SequencePhenotype("ATCT");
				history[2] = new SequencePhenotype("ATCA");
				break;
			case "aminoAcids":
				simplePheno = new SequencePhenotype("EQGZ");

				history[0] = new SequencePhenotype("ESGZ");
				history[1] = new SequencePhenotype("ESVZ");
				history[2] = new SequencePhenotype("ESVY");
				break;
		}
	}


	/**
	 * Test constructors for the SequencePhenotype class.
	 */
	@Test
	public void testConstructors() {
		// For now, assert that an emptyConstructor returns a single nucleotide
		assertEquals(Parameters.startingSequence.length(), emptyPheno.getSequence().length());

		// Make sure simple constructor is fine
		switch (Parameters.alphabet) {
			case "nucleotides":
				assertEquals("ACGT", simplePheno.getSequence());
				break;
			case "aminoAcids":
				assertEquals("EQGZ", simplePheno.getSequence());
				break;
		}
	}

	/**
	 * Test hamming distance calculations and exception throws.
	 */
	@Test
	public void testDistance() {
		SequencePhenotype testPheno;
		SequencePhenotype diffPheno;
		Exception exception;
		SequencePhenotype shortPheno;

		switch (Parameters.alphabet) {
			case "nucleotides":
				// Hamming distance of same sequence
				testPheno = new SequencePhenotype("ACGT");
				assertEquals(0, simplePheno.distance(testPheno), 0.0);

				// Hamming distance of differing sequences
				diffPheno = new SequencePhenotype("ATGT");
				assertEquals(1, simplePheno.distance(diffPheno), 0.0);

				// Hamming distance of sequences of different length
				shortPheno = new SequencePhenotype("AT");
				exception = assertThrows(Exception.class, () -> testPheno.distance(shortPheno));
				assertEquals("Sequence lengths are not equal!", exception.getMessage());
				break;
			case "aminoAcids":
				// Hamming distance of same sequence
				testPheno = new SequencePhenotype("EQGZ");
				assertEquals(0, simplePheno.distance(testPheno), 0.0);

				// Hamming distance of differing sequences
				diffPheno = new SequencePhenotype("ESGZ");
				assertEquals(1, simplePheno.distance(diffPheno), 0.0);

				// Hamming distance of sequences of different length
				shortPheno = new SequencePhenotype("EQ");
				exception = assertThrows(Exception.class, () -> testPheno.distance(shortPheno));
				assertEquals("Sequence lengths are not equal!", exception.getMessage());
				break;
		}
	}

	/**
	 * Test getter method.
	 */
	@Test
	public void testGetSequence() {
		// Assert the sequence was updated

		switch (Parameters.alphabet) {
			case "nucleotides":
				assertEquals("ACGT", simplePheno.getSequence());
				break;
			case "aminoAcids":
				assertEquals("EQGZ", simplePheno.getSequence());
				break;
		}
	}

	/**
	 * Test mutate function.
	 * Original implementation adds a random base to the end of the sequence.
	 * Another version will mutate a base at random.
	 */
	@Test
	public void testMutate() {
		SequencePhenotype originalPheno;
		switch (Parameters.alphabet) {
			case "nucleotides":
				originalPheno = new SequencePhenotype("ACGT");
				// Test new substitution mutate.
				originalPheno = originalPheno.mutate();
				assertNotEquals(simplePheno.getSequence(), originalPheno.getSequence());
				assertEquals(simplePheno.getSequence().length(), originalPheno.getSequence().length());
				break;
			case "aminoAcids":
				originalPheno = new SequencePhenotype("EQGZ");
				// Test new substitution mutate.
				originalPheno = originalPheno.mutate();
				assertNotEquals(simplePheno.getSequence(), originalPheno.getSequence());
				assertEquals(simplePheno.getSequence().length(), originalPheno.getSequence().length());
				break;
		}
	}

	/**
	 * Test riskOfInfection calculations.
	 */
	@Test
	public void testRiskOfInfection() {
		// Add this test after implementing a new (or old) CI function.
	}

	/**
	 * A redundant test for returning the sequence, but with toString().
	 * This may change a bit as the class develops (i.e., specific epitope sequences)
	 */
	@Test
	public void testToString() {
		switch (Parameters.alphabet) {
			case "nucleotides":
				assertEquals("ACGT", simplePheno.toString());
				break;
			case "aminoAcids":
				assertEquals("EGQZ", simplePheno.toString());
				break;
		}
	}

	/**
	 * Test equals()
	 */
	@Test
	public void testEquals() {
		// Same SequencePhenotype objects are equal
		assertEquals(emptyPheno, emptyPheno);

		SequencePhenotype simplePhenoSame;
		SequencePhenotype simplePhenoDifferent;

		switch (Parameters.alphabet) {
			case "nucleotides":
				simplePhenoSame = new SequencePhenotype("ACGT");
				simplePhenoDifferent = new SequencePhenotype("CCGT");

				// SequencePhenotype objects with the same sequence are equal
				assertEquals(simplePheno, simplePhenoSame);
				// SequencePhenotype objects with different sequences are not equal
				assertNotEquals(simplePheno, simplePhenoDifferent);
				break;
			case "aminoAcids":
				simplePhenoSame = new SequencePhenotype("EGQZ");
				simplePhenoDifferent = new SequencePhenotype("ESVY");

				// SequencePhenotype objects with the same sequence are equal
				assertEquals(simplePheno, simplePhenoSame);
				// SequencePhenotype objects with different sequences are not equal
				assertNotEquals(simplePheno, simplePhenoDifferent);
				break;
		}
	}

}
