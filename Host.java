/* A human individual that harbors viruses and immunity */

import java.util.*;
import java.io.*;
import java.util.regex.*;

public class Host {

	// fields
	private Virus infection;												
	private Phenotype[] immuneHistory = new Phenotype[0];
	private double birth;
	
	// naive host
	public Host() {
		double lifeExpectancy = (1.0 / Parameters.birthRate) / 365.0;
		birth = -1.0 * (Parameters.burnin / 365.0) - Random.nextExponential(lifeExpectancy);		
		initializeHistory();	
	}
	
	// initial infected host
	public Host(Virus v) {
		double lifeExpectancy = (1.0 / Parameters.birthRate) / 365.0;
		birth = -1.0 * (Parameters.burnin / 365.0) - Random.nextExponential(lifeExpectancy);
		infection = v;
		initializeHistory();
	}
	
	// checkpointed host
	public Host(int d, String sVirus, String sHist) {
		double lifeExpectancy = (1.0 / Parameters.birthRate) / 365.0;
		birth = -1.0 * (Parameters.burnin / 365.0) - Random.nextExponential(lifeExpectancy);
		if (!sVirus.equals("n")) {
			Pattern rc = Pattern.compile(",");
    		String[] traitList = rc.split(sVirus);
    		double x = Double.parseDouble(traitList[0]);
    		double y = Double.parseDouble(traitList[1]);
			Phenotype p = PhenotypeFactory.makeArbitaryPhenotype(x,y);
			infection = new Virus(Parameters.urVirus,d,p);
		}
		if (!sHist.equals("n")) {
			Pattern rsc = Pattern.compile(";");
    		String[] phenotypeList = rsc.split(sHist);
    		for (int i = 0; i < phenotypeList.length; i++) {
				Pattern rc = Pattern.compile(",");
    			String[] traitList = rc.split(phenotypeList[i]);
    			double x = Double.parseDouble(traitList[0]);
    			double y = Double.parseDouble(traitList[1]);
				Phenotype p = PhenotypeFactory.makeArbitaryPhenotype(x,y);
				addToHistory(p);
			}
		}		
	}
	
	// sometimes start with immunity	
	public void initializeHistory() {
		double chanceOfSuccess = Parameters.initialPrR;
		if (Random.nextBoolean(chanceOfSuccess)) {	
			Phenotype p = Parameters.urImmunity;
			addToHistory(p);
		}	
	}
	
	public void addToHistory(Phenotype p) {
		Phenotype[] newHistory = new Phenotype[immuneHistory.length + 1];
		for (int i = 0; i < immuneHistory.length; i++) {
			newHistory[i] = immuneHistory[i];
		}
		newHistory[immuneHistory.length] = p;
		immuneHistory = newHistory;
	}
	
	// infection methods
	public void reset() {
		infection = null;
		immuneHistory = new Phenotype[0];
		birth = Parameters.getDate();
	}
	public double getBirth() {
		return birth;
	}
	public double getAge() {
		return Parameters.getDate() - birth;
	}
	public boolean isAdult() {
		boolean adult = false;
		if (getAge() >= 15) {
			adult = true;
		}
		return adult;
	}
	public boolean isInfected() {
		boolean infected = false;
		if (infection != null) {
			infected = true;
		}
		return infected;
	}
	public Virus getInfection() {
		return infection;
	}
	public void infect(Virus pV, int d) {
		Virus nV = new Virus(pV, d);
		infection = nV;
	}
	public void clearInfection() {
		Phenotype p = infection.getPhenotype();
		addToHistory(p);
		infection = null;
	}
	public int getHistoryLength() {
		return immuneHistory.length;
	}
	
	// make a new virus with the mutated phenotype
	public void mutate() {
		Virus mutV = infection.mutate();
		infection = mutV;
	}
	
	// history methods
	public Phenotype[] getHistory() {
		return immuneHistory;
	}	
	
	public void printHistory() {
		for (int i = 0; i < immuneHistory.length; i++) {
			System.out.println(immuneHistory[i]);
		}
	}

	public void printInfection(PrintStream stream) {
		if (infection != null) {
			stream.print(infection.getPhenotype());
		}
		else {
			stream.print("n");
		}
	}
	
	public void printHistory(PrintStream stream) {
		if (immuneHistory.length > 0) {
			stream.print(immuneHistory[0]);
			for (int i = 1; i < immuneHistory.length; i++) {
				stream.print(";" + immuneHistory[i]);
			}
		}
		else {
			stream.print("n");
		}
	}	
		
	public String toString() {
		return Integer.toHexString(this.hashCode());
	}	
	
}