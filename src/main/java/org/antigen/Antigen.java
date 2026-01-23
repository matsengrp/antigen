package org.antigen;

/* Implements an individual-based model in which the infection's genealogical history is tracked through time */

import java.io.FileNotFoundException;
import org.antigen.core.*;

public class Antigen {
  public static void main(String[] args) throws FileNotFoundException {

    // initialize random number generator
    cern.jet.random.AbstractDistribution.makeDefaultGenerator();

    // initialize static parameters
    Parameters.load();
    Parameters.initialize();

    // run simulation
    Simulation sim = new Simulation();
    sim.run();
  }
}
