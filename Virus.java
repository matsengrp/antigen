/* Virus infection that has genotype, phenotype and ancestry */

import java.util.*;

public class Virus {

	// simulation fields
	private Virus parent;
	private Phenotype phenotype;
	private double birth;		// measured in years relative to burnin
	private int deme;
	private double fitness; 	// seasonal fitness (averageInfectionRisk * seasonality * probSusceptible)
	private double averageInfectionRisk; 	// raw average risk of infection
	private double probSusceptible; 	// fraction of susceptible hosts in deme
	private double demeSeasonality; 	// seasonality value for deme
	
	// additional reconstruction fields
	private boolean marked;
	private boolean trunk;	// fill this at the end of the simulation
	private List<Virus> children = new ArrayList<>(0);	// will be void until simulation ends
	private double layout;
	private int coverage;		// how many times this Virus has been covered in tracing the tree backwards
	
	// initialization
	public Virus() {
		phenotype = PhenotypeFactory.makeVirusPhenotype();
		birth = Parameters.getDate();
	}
		
	// replication, copies the virus, but remembers the ancestry
	public Virus(Virus v, int d) {
		parent = v;
		phenotype = v.getPhenotype();
		birth = Parameters.getDate();
		deme = d;
	}
	
	public Virus(Virus v, int d, Phenotype p) {
		parent = v;
		phenotype = p;
		birth = Parameters.getDate();
		deme = d;
	}	
	
	public Virus(int d, Phenotype p) {
		parent = null;
		phenotype = p;
		birth = Parameters.getDate();
		deme = d;
	}		
	
	// methods
	public Phenotype getPhenotype() {
		return phenotype;
	}
	public void setPhenotype(Phenotype p) {
		phenotype = p;
	}	
	public double getBirth() {
		return birth;
	}
	public Virus getParent() {
		return parent;
	}
	public void setParent(Virus v) {
		parent = v;
	}
	public boolean isTrunk() {
		return trunk; 
	}
	public void makeTrunk() {
		trunk = true;
	}
	public void mark() {
		marked = true;
	}
	public boolean isMarked() {
		return marked;
	}
	public int getDeme() {
		return deme;
	}	
	public double getLayout() {
		return layout;
	}
	public void setLayout(double y) {
		layout = y;
	}
	public int getCoverage() {
		return coverage;
	}
	public void incrementCoverage() {
		coverage++;
	}
	public double getFitness(){
		return fitness;
	}
	public void setFitness(double f){
		fitness = f;
	}
	public double getAverageInfectionRisk(){
		return averageInfectionRisk;
	}
	public void setAverageInfectionRisk(double risk){
		averageInfectionRisk = risk;
	}
	public double getProbSusceptible(){
		return probSusceptible;
	}
	public void setProbSusceptible(double prob){
		probSusceptible = prob;
	}
	public double getDemeSeasonality(){
		return demeSeasonality;
	}
	public void setDemeSeasonality(double seasonality){
		demeSeasonality = seasonality;
	}
	
	// add virus node as child if does not already exist
	public void addChild(Virus v) {
		if (!children.contains(v)) {
			children.add(v);
		}
	}		
	public int getNumberOfChildren() {
		return children.size();
	}
	public List<Virus> getChildren() {
		return children;
	}	
	public boolean isTip() {
		return getNumberOfChildren() == 0;
	}
	
	// returns a mutated copy, original virus left intact
	public Virus mutate() {
	
		Phenotype mutP = phenotype.mutate();			// mutated copy
		Virus mutV = new Virus(this,deme,mutP);
		return mutV;
		
	}
	
	public Virus commonAncestor(Virus virusB) {
		// Algorithm proceeds by recursively visiting parents.
		// It terminates when either lineage arrives at a node already in the ancestry set,
		// which must have been already visited by the other lineage and thus represents
		// a common ancestor.
		
		assert(virusB != null);
		if(virusB == this) {
			return this;
		}
		
		Virus lineageA = this;
		Virus lineageB = virusB;
		Set<Virus> ancestrySet = new HashSet<>();
		ancestrySet.add(lineageA);
		ancestrySet.add(lineageB);
		while(lineageA != null || lineageB != null) {
			if(lineageA != null) {
				lineageA = lineageA.getParent();
				if(lineageA != null && !ancestrySet.add(lineageA)) {
					return lineageA;
				}
			}
			if(lineageB != null) {
				lineageB = lineageB.getParent();
				if(lineageB != null && !ancestrySet.add(lineageB)) {
					return lineageB;
				}
			}
		}
		return null;
	}
	
	public double distance(Virus virusB) {
		Virus ancestor = commonAncestor(virusB);
		if (ancestor != null) {
			double distA = getBirth() - ancestor.getBirth();
			double distB = virusB.getBirth() - ancestor.getBirth();
			return distA + distB;
		}
		else {
			return 0;
		}
	}
	
	public double antigenicDistance(Virus virusB) {
		return phenotype.distance(virusB.getPhenotype());
	}	
	
	// is there a coalescence event within x amount of time? (measured in years)
	public double coalescence(Virus virusB, double windowTime) {

		Virus lineageA = this;
		Virus lineageB = virusB;
		Set<Virus> ancestry = new HashSet<>();
		double success = 0.0;
		
		double startTime = lineageA.getBirth();
		double time = startTime;
		while (time > startTime - windowTime) {
			if (lineageA.getParent() != null) {		
				lineageA = lineageA.getParent();
				time = lineageA.getBirth();
				ancestry.add(lineageA);
			}
			else {
				break;
			}
		}
		
		startTime = lineageB.getBirth();
		time = startTime;
		while (time > startTime - windowTime) {
			if (lineageB.getParent() != null) {		
				lineageB = lineageB.getParent();
				time = lineageB.getBirth();				
				if (!ancestry.add(lineageB)) { 
					success = 1.0;
					break; 
				}
			}
			else {
				break;
			}			
		}
		
		return success;	

	}	
	
	// this is the interval from this virus's birth back to its parent's birth
	public double serialInterval() {
		Virus p = getParent();
		return getBirth() - p.getBirth();
	}
	
	public String toString() {
		return Integer.toHexString(this.hashCode());
	}

}