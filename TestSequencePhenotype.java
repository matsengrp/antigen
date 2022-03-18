import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

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
		emptyPheno = new SequencePhenotype();
		simplePheno = new SequencePhenotype("ACGT");
		// Create a small history of Sequence phenotypes
		SequencePhenotype [] history = new SequencePhenotype[3];
		history[0] = new SequencePhenotype("ATGT");
		history[1] = new SequencePhenotype("ATCT");
		history[2] = new SequencePhenotype("ATCA");
	}


	/**
	 * Test constructors for the SequencePhenotype class.
	 */
	@Test
	public void testConstructors() {
		// For now, assert that an emptyConstructor returns a single nucleotide
		assertEquals(1, emptyPheno.getSequence().length());

		// Make sure simple constructor is fine
		assertEquals("ACGT", simplePheno.getSequence());
	}

	/**
	 * Test hamming distance calculations and exception throws.
	 */
	@Test
	public void testDistance() {
		// Hamming distance of same sequence
		SequencePhenotype testPheno = new SequencePhenotype("ACGT");
		assertEquals(0, simplePheno.distance(testPheno), 0.0);

		// Hamming distance of differing sequences
		SequencePhenotype diffPheno = new SequencePhenotype("ATGT");
		assertEquals(1, simplePheno.distance(diffPheno), 0.0);

		// Hamming distance of sequences of different length
		SequencePhenotype shortPheno = new SequencePhenotype("AT");
		Exception exception = assertThrows(Exception.class, () -> testPheno.distance(shortPheno));
		assertEquals("Sequence lengths are not equal!", exception.getMessage());
	}

	/**
	 * Test getter and setter methods.
	 */
	@Test
	public void testSetAndGetSequence() {
		// Set sequence of testPheno to something different
		simplePheno.setSequence("TCGA");

		// Assert the sequence was updated
		assertEquals("TCGA", simplePheno.getSequence());
	}

	/**
	 * Test mutate function.
	 * Original implementation adds a random base to the end of the sequence.
	 * Another version will mutate a base at random.
	 */
	@Test
	public void testMutate() {
		// For now, test Thien's original version.
		SequencePhenotype mutantPheno = simplePheno.mutate();
		assertNotEquals(simplePheno.getSequence(), mutantPheno.getSequence());
		assertNotEquals(simplePheno.getSequence().length(), mutantPheno.getSequence().length());

		SequencePhenotype subPheno = simplePheno.mutate();
		assertNotEquals(simplePheno.getSequence(), subPheno.getSequence());

		// Test that will fail until the substitution is put in.
		assertEquals(simplePheno.getSequence().length(), subPheno.getSequence().length());
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
		assertEquals("ACGT", simplePheno.toString());
	}

	/**
	 * Test equals()
	 */
	@Test
	public void testEquals() {
		// Same SequencePhenotype objects are equal
		assertTrue(emptyPheno.equals(emptyPheno));

		SequencePhenotype simplePhenoTest = new SequencePhenotype("ACGT");

		// SequencePhenotype objects with the same sequence are equal
		assertTrue(simplePheno.equals(simplePhenoTest));
	}

}
