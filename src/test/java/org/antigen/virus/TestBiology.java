package org.antigen.virus;

import static org.junit.Assert.*;

import org.antigen.core.Parameters;
import org.junit.Before;
import org.junit.Test;

/** Test class for Biology.java, specifically the K80 DNA evolution model */
public class TestBiology {

  @Before
  public void setUp() {
    // Set a known transition/transversion ratio for testing
    Parameters.transitionTransversionRatio = 5.0;
  }

  /** Test that K80 model probabilities are calculated correctly */
  @Test
  public void testK80ModelProbabilities() {
    Biology.K80DNAEvolutionModel model = Biology.K80DNAEvolutionModel.MUTATION;

    // With ratio = 5.0:
    // transversionProb = 1/(2+5) = 1/7 ≈ 0.1428571
    // transitionProb = 5/7 ≈ 0.7142857

    double expectedTransversionProb = 1.0 / 7.0;
    double expectedTransitionProb = 5.0 / 7.0;

    // Test nucleotide A (transitions to G, transversions to C and T)
    double[] probsA = model.transitionTranversionProbability.get('A');
    assertEquals(0.0, probsA[0], 0.0001); // A->A (no change)
    assertEquals(expectedTransversionProb, probsA[1], 0.0001); // A->C (transversion)
    assertEquals(
        expectedTransversionProb + expectedTransitionProb, probsA[2], 0.0001); // A->G (transition)
    assertEquals(1.0, probsA[3], 0.0001); // Cumulative total

    // Test nucleotide C (transitions to T, transversions to A and G)
    double[] probsC = model.transitionTranversionProbability.get('C');
    assertEquals(expectedTransversionProb, probsC[0], 0.0001); // C->A (transversion)
    assertEquals(expectedTransversionProb, probsC[1], 0.0001); // C->C (no change)
    assertEquals(
        expectedTransversionProb + expectedTransitionProb,
        probsC[2],
        0.0001); // C->G (transversion)
    assertEquals(1.0, probsC[3], 0.0001); // C->T (transition)

    // Test nucleotide G (transitions to A, transversions to C and T)
    double[] probsG = model.transitionTranversionProbability.get('G');
    assertEquals(expectedTransitionProb, probsG[0], 0.0001); // G->A (transition)
    assertEquals(
        expectedTransitionProb + expectedTransversionProb,
        probsG[1],
        0.0001); // G->C (transversion)
    assertEquals(
        expectedTransitionProb + expectedTransversionProb, probsG[2], 0.0001); // G->G (no change)
    assertEquals(1.0, probsG[3], 0.0001); // G->T (transversion)

    // Test nucleotide T (transitions to C, transversions to A and G)
    double[] probsT = model.transitionTranversionProbability.get('T');
    assertEquals(expectedTransversionProb, probsT[0], 0.0001); // T->A (transversion)
    assertEquals(
        expectedTransversionProb + expectedTransitionProb, probsT[1], 0.0001); // T->C (transition)
    assertEquals(1.0, probsT[2], 0.0001); // T->G (transversion)
    assertEquals(1.0, probsT[3], 0.0001); // T->T (no change)
  }

  /** Test that mutation sampling respects the transition/transversion ratio */
  @Test
  public void testMutationSampling() {
    Biology.K80DNAEvolutionModel model = Biology.K80DNAEvolutionModel.MUTATION;

    // Sample many mutations and check that the ratio is approximately correct
    int numSamples = 10000;
    int transitionsFromA = 0;
    int transversionsFromA = 0;

    // Set seed for reproducibility
    java.util.Random rand = new java.util.Random(12345);

    for (int i = 0; i < numSamples; i++) {
      // Manually sample using the probability array
      double randomValue = rand.nextDouble();
      double[] probs = model.transitionTranversionProbability.get('A');

      char mutant = 'A';
      if (randomValue < probs[1]) {
        mutant = 'C'; // transversion
        transversionsFromA++;
      } else if (randomValue < probs[2]) {
        mutant = 'G'; // transition
        transitionsFromA++;
      } else {
        mutant = 'T'; // transversion
        transversionsFromA++;
      }
    }

    // The ratio should be approximately 5.0
    // But since we have 2 transversions and 1 transition,
    // the observed ratio should be transition/(transversions/2)
    double observedRatio = (double) transitionsFromA / (transversionsFromA / 2.0);
    assertEquals(5.0, observedRatio, 0.5); // Allow some statistical variance
  }
}
