/* Virus infection that has genotype, phenotype and ancestry */

import java.util.*;

public class Virus {

	// simulation fields
	private Virus parent;
	private Phenotype phenotype;
	private double birth;	// measured in years relative to burnin
	private boolean trunk;	// fill this at the end of the simulation
	private int deme;

	// additional reconstruction fields
	private List<Virus> children = new ArrayList<Virus>();	// will be void until simulation ends	
	private double layout;
	
	// initialization
	public Virus() {
		phenotype = PhenotypeFactory.makeVirusPhenotype();
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
	public boolean isTrunk() {
		return trunk; 
	}
	public void makeTrunk() {
		trunk = true;
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
		return getNumberOfChildren() == 0 ? true : false;
	}
	
	// returns a mutated copy, original virus left intact
	public Virus mutate() {
	
		Phenotype mutP = phenotype.mutate();			// mutated copy
		Virus mutV = new Virus(this,deme,mutP);
		return mutV;
		
	}
	
	public Virus commonAncestor(Virus virusB) {
		
		// go through current virus's history and add to a set
		Virus lineage = this;
		Set<Virus> ancestry = new HashSet<Virus>();
		while (lineage.getParent() != null) {
			lineage = lineage.getParent();
			ancestry.add(lineage);
		}
		
		// go through other virus's history and add to this set, stop when duplicate is encountered and return this duplicate
		lineage = virusB;
		while (lineage.getParent() != null) {
			lineage = lineage.getParent();
			if (!ancestry.add(lineage)) {
				break;
			}
		}		
		
		
		return lineage;
		
	}
	
	public double distance(Virus virusB) {
		Virus ancestor = commonAncestor(virusB);
		double distA = getBirth() - ancestor.getBirth();
		double distB = virusB.getBirth() - ancestor.getBirth();
		return distA + distB;
	}
	
	public String toString() {
		return Integer.toHexString(this.hashCode());
	}

}