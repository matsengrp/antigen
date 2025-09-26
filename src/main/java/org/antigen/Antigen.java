package org.antigen;

/* Implements an individual-based model in which the infection's genealogical history is tracked through time */

import java.io.FileNotFoundException;
import org.antigen.core.*;

public class Antigen {
    public static void main(String[] args) throws FileNotFoundException {

		// initialize random number generator
		cern.jet.random.AbstractDistribution.makeDefaultGenerator();
		
		// initialize static parameters
		if (args.length > 0) {
			// Use command line argument for parameter file
			System.out.println("Loading parameters from command line argument: " + args[0]);
			Parameters.load(args[0]);
		} else {
			// Fall back to default embedded parameters
			System.out.println("No parameter file specified, using default embedded parameters.yml");
			Parameters.load();
		}
		Parameters.initialize();
		
		// run simulation
		Simulation sim = new Simulation();
		sim.run();	
		
	}
   	
}

