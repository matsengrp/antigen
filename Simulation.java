/* Simulation functions, holds the host population */

import java.util.*;
import java.io.*;

import com.javamex.classmexer.*;

public class Simulation {
	public static PrintStream mutations;

	static {
		try {
			mutations = new PrintStream("mutation.csv");
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	public static Map<Character, double[]> transitionTranversionProbability = new HashMap<>()  {{
		// key: nucleotide
		// value: int array of probability that key will mutate to the letter in the corresponding index
		// {A, G, T, C}
		put('A', new double[]{0, 0.5, 0.75, 1.0});
		put('G', new double[]{0.5, 0.5, 0.75, 1.0});
		put('T', new double[]{0.25, 0.5, 0.5, 1.0});
		put('C', new double[]{0.25, 0.5, 1.0, 1.0});
	}};

	public static Map<String, String> codonMap = new HashMap<>() {{
		put("TTT", "F");
		put("TTC", "F");
		put("TTA", "L");
		put("TTG", "L");

		put("CTT", "L");
		put("CTC", "L");
		put("CTA", "L");
		put("CTG", "L");

		put("ATT", "I");
		put("ATC", "I");
		put("ATA", "I");
		put("ATG", "M");

		put("GTT", "V");
		put("GTC", "V");
		put("GTA", "V");
		put("GTG", "V");


		put("TCT", "S");
		put("TCC", "S");
		put("TCA", "S");
		put("TCG", "S");

		put("CCT", "P");
		put("CCC", "P");
		put("CCA", "P");
		put("CCG", "P");

		put("ACT", "T");
		put("ACC", "T");
		put("ACA", "T");
		put("ACG", "T");

		put("GCT", "A");
		put("GCC", "A");
		put("GCA", "A");
		put("GCG", "A");


		put("TAT", "Y");
		put("TAC", "Y");
		put("TAA", "STOP");
		put("TAG", "STOP");

		put("CAT", "H");
		put("CAC", "H");
		put("CAA", "Q");
		put("CAG", "Q");

		put("AAT", "N");
		put("AAC", "N");
		put("AAA", "K");
		put("AAG", "K");

		put("GAT", "D");
		put("GAC", "D");
		put("GAA", "E");
		put("GAG", "E");


		put("TGT", "C");
		put("TGC", "C");
		put("TGA", "STOP");
		put("TGG", "W");

		put("CGT", "R");
		put("CGC", "R");
		put("CGA", "R");
		put("CGG", "R");

		put("AGT", "S");
		put("AGC", "S");
		put("AGA", "T");
		put("AGG", "T");

		put("GGT", "G");
		put("GGC", "G");
		put("GGA", "G");
		put("GGG", "G");
	}};

	// fields
	private List<HostPopulation> demes = new ArrayList<>();
	private double diversity;
	private double tmrca;
	private double netau;
	private double serialInterval;
	private double antigenicDiversity;

	private List<Double> diversityList = new ArrayList<>();
	private List<Double> tmrcaList = new ArrayList<>();
	private List<Double> netauList = new ArrayList<>();
	private List<Double> serialIntervalList = new ArrayList<>();
	private List<Double> antigenicDiversityList = new ArrayList<>();
	private List<Double> nList = new ArrayList<>();
	private List<Double> sList = new ArrayList<>();
	private List<Double> iList = new ArrayList<>();
	private List<Double> rList = new ArrayList<>();
	private List<Double> casesList = new ArrayList<>();

	public static int matrixSize = Parameters.AlphabetType.AMINO_ACIDS.getValidCharacters().length();
	public static Map<Integer, double[][][]> siteMutationVectors = new HashMap<>();

	// constructor
	public Simulation() {
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp;
			if (Parameters.restartFromCheckpoint) {
				hp = new HostPopulation(i, true);
			} else {
				hp = new HostPopulation(i);
			}
			demes.add(hp);
		}

		initializeSiteMutationVectors();

	}

	private void initializeSiteMutationVectors() {
		for (int siteNumber = 0; siteNumber < Parameters.startingSequence.length(); siteNumber++) {
			double[][][] currentSiteMutationMatrix = new double[matrixSize][matrixSize][];
			for (int wildTypeIndex = 0; wildTypeIndex < matrixSize; wildTypeIndex++) {
				for (int mutationIndex = 0; mutationIndex < matrixSize; mutationIndex++) {
					if (mutationIndex < wildTypeIndex) { // update lower and upper triangle in this branch.
						// direction of mutation
						double theta;
						if (Parameters.mut2D) {
							theta = Random.nextDouble(0, 2 * Math.PI);
						} else {
							if (Random.nextBoolean(0.5)) {
								theta = 0;
							} else {
								theta = Math.PI;
							}
						}

						// size of mutation
						double r = Parameters.meanStep;
						if (!Parameters.fixedStep) {
							double alpha = (Parameters.meanStep * Parameters.meanStep) / (Parameters.sdStep * Parameters.sdStep);
							double beta = (Parameters.sdStep * Parameters.sdStep) / Parameters.meanStep;
							r = Random.nextGamma(alpha, beta);
						}

						// create phenotype
						double mutA = r * Math.cos(theta);
						double mutB = r * Math.sin(theta);

						double[] mutations = new double[]{mutA, mutB};

						currentSiteMutationMatrix[wildTypeIndex][mutationIndex] = mutations;

						mutations = new double[]{-1 * mutA, -1 * mutB};

						currentSiteMutationMatrix[mutationIndex][wildTypeIndex] = mutations;
					} else if (mutationIndex == wildTypeIndex) {
						currentSiteMutationMatrix[wildTypeIndex][mutationIndex] = new double[]{0.0, 0.0};
					}
				}
			}

			siteMutationVectors.put(siteNumber, currentSiteMutationMatrix);
		}
	}

	// methods

	public int getN() {
		int count = 0;
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getN();
		}
		return count;
	}

	public int getS() {
		int count = 0;
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getS();
		}
		return count;
	}

	public int getI() {
		int count = 0;
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getI();
		}
		return count;
	}

	public int getR() {
		int count = 0;
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getR();
		}
		return count;
	}

	public int getCases() {
		int count = 0;
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getCases();
		}
		return count;
	}

	public double getDiversity() {
		return diversity;
	}

	public double getNetau() {
		return netau;
	}

	public double getTmrca() {
		return tmrca;
	}

	public double getSerialInterval() {
		return serialInterval;
	}

	public double getAntigenicDiversity() {
		return antigenicDiversity;
	}

	// proportional to infecteds in each deme
	public int getRandomDeme() {
		int n = Random.nextInt(0,getN()-1);
		int d = 0;
		int target = (demes.get(0)).getN();
		while (n < target) {
			d += 1;
			target += (demes.get(d)).getN();
		}
		return d;
	}

	// return random virus proportional to worldwide prevalence
	public Virus getRandomInfection() {

		Virus v = null;

		if (getI() > 0) {

			// get deme proportional to prevalence
			int n = Random.nextInt(0,getI()-1);
			int d = 0;
			int target = (demes.get(0)).getI();
			while (d < Parameters.demeCount) {
				if (n < target) {
					break;
				} else {
					d++;
					target += (demes.get(d)).getI();
				}
			}
			HostPopulation hp = demes.get(d);

			// return random infection from this deme
			if (hp.getI()>0) {
				Host h = hp.getRandomHostI();
				v = h.getInfection();
			}

		}

		return v;

	}

	// return random host from random deme
	public Host getRandomHost() {
		int d = Random.nextInt(0,Parameters.demeCount-1);
		HostPopulation hp = demes.get(d);
		return hp.getRandomHost();
	}

	public double getAverageRisk(Phenotype p) {

		double averageRisk = 0;
		for (int i = 0; i < 10000; i++) {
			Host h = getRandomHost();
			Phenotype[] history = h.getHistory();
			averageRisk += p.riskOfInfection(history);
		}
		averageRisk /= 10000.0;
		return averageRisk;

	}

	public void printImmunity() {

		try {
			File immunityFile = new File("out.immunity");
			immunityFile.delete();
			immunityFile.createNewFile();
			PrintStream immunityStream = new PrintStream(immunityFile);

			for (double x = VirusTree.xMin; x <= VirusTree.xMax; x += 0.5) {
				for (double y = VirusTree.yMin; y <= VirusTree.yMax; y += 0.5) {

					Phenotype p = PhenotypeFactory.makeArbitaryPhenotype(x,y);
					double risk = getAverageRisk(p);
					immunityStream.printf("%.4f,", risk);

				}
				immunityStream.println();
			}

			immunityStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file");
			System.exit(0);
		}

	}

	public void printHostPopulation() {

		try {
			File hostFile = new File("out.hosts");
			hostFile.delete();
			hostFile.createNewFile();
			PrintStream hostStream = new PrintStream(hostFile);
			for (int i = 0; i < Parameters.demeCount; i++) {
				HostPopulation hp = demes.get(i);
				hp.printHostPopulation(hostStream);
			}
			hostStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file");
			System.exit(0);
		}

	}

	public void makeTrunk() {
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.makeTrunk();
		}
	}

	public void printState() {

		System.out.printf("%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%d\t%d\n", (int) Parameters.day, getDiversity(), getTmrca(),  getNetau(), getSerialInterval(), getAntigenicDiversity(), getN(), getS(), getI(), getR(), getCases());

		if (Parameters.memoryProfiling && Parameters.day % 10 == 0) {
			long noBytes = MemoryUtil.deepMemoryUsageOf(this);
			System.out.println("Total: " + noBytes);
			HostPopulation hp = demes.get(1);
			noBytes = MemoryUtil.deepMemoryUsageOf(hp);
			System.out.println("One host population: " + noBytes);
			Host h = hp.getRandomHostS();
			noBytes = MemoryUtil.deepMemoryUsageOf(h);
			System.out.println("One susceptible host with " +  h.getHistoryLength() + " previous infection: " + noBytes);
			//h.printHistory();
			if (getI() > 0) {
				Virus v = getRandomInfection();
				noBytes = MemoryUtil.memoryUsageOf(v);
				System.out.println("One virus: " + noBytes);
				noBytes = MemoryUtil.deepMemoryUsageOf(VirusTree.getTips());
				System.out.println("Virus tree: " + noBytes);
			}
		}

	}

	public void printHeader(PrintStream stream) {
		stream.print("date\tdiversity\ttmrca\tnetau\tserialInterval\tantigenicDiversity\ttotalN\ttotalS\ttotalI\ttotalR\ttotalCases");
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.printHeader(stream);
		}
		stream.println();
	}

	public void printState(PrintStream stream) {
		stream.printf("%.4f\t%.4f\t%.4f\t%.4f\t%.5f\t%.4f\t%d\t%d\t%d\t%d\t%d", Parameters.getDate(), getDiversity(), getTmrca(), getNetau(), getSerialInterval(), getAntigenicDiversity(), getN(), getS(), getI(), getR(), getCases());
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.printState(stream);
		}
		stream.println();
	}

	public void printSummary() {

		try {
			File summaryFile = new File("out.summary");
			summaryFile.delete();
			summaryFile.createNewFile();
			PrintStream summaryStream = new PrintStream(summaryFile);
			summaryStream.printf("parameter\tfull\n");
			summaryStream.printf("endDate\t%.4f\n", Parameters.getDate());
			summaryStream.printf("diversity\t%.4f\n", mean(diversityList));
			summaryStream.printf("tmrca\t%.4f\n", mean(tmrcaList));
			summaryStream.printf("netau\t%.4f\n", mean(netauList));
			summaryStream.printf("serialInterval\t%.5f\n", mean(serialIntervalList));
			summaryStream.printf("antigenicDiversity\t%.4f\n", mean(antigenicDiversityList));
			summaryStream.printf("N\t%.4f\n", mean(nList));
			summaryStream.printf("S\t%.4f\n", mean(sList));
			summaryStream.printf("I\t%.4f\n", mean(iList));
			summaryStream.printf("R\t%.4f\n", mean(rList));
			summaryStream.printf("cases\t%.4f\n", mean(casesList));

			summaryStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file");
			System.exit(0);
		}

	}

	private double mean(List<Double> list) {
		double mean = 0;
		if(!list.isEmpty()) {
			for (Double item : list) {
				mean += item;
			}
			mean /= list.size();
		}
		return mean;
	}

	public void updateDiversity() {

		diversity = 0.0;
		tmrca = 0.0;
		antigenicDiversity = 0.0;
		netau = 0.0;
		serialInterval = 0.0;

		double coalCount = 0.0;
		double coalOpp = 0.0;
		double coalWindow = Parameters.netauWindow / 365.0;
		int sampleCount = Parameters.diversitySamplingCount;

		for (int i = 0; i < sampleCount; i++) {
			Virus vA = getRandomInfection();
			Virus vB = getRandomInfection();
			if (vA != null && vB != null) {
				double dist = vA.distance(vB);
				diversity += dist;
				if (dist > tmrca) {
					tmrca = dist;
				}
				antigenicDiversity += vA.antigenicDistance(vB);
				coalOpp += coalWindow;
				coalCount += vA.coalescence(vB, coalWindow);
				serialInterval += vA.serialInterval();
			}
		}

		diversity /= sampleCount;
		tmrca /= 2.0;
		netau = coalOpp / coalCount;
		serialInterval /= sampleCount;
		antigenicDiversity /= sampleCount;

	}

	public void pushLists() {
		diversityList.add(diversity);
		tmrcaList.add(tmrca);
		netauList.add(netau);
		serialIntervalList.add(serialInterval);
		antigenicDiversityList.add(antigenicDiversity);
		nList.add((double) getN());
		sList.add((double) getS());
		iList.add((double) getI());
		rList.add((double) getR());
		casesList.add((double) getCases());
	}

	public void resetCases() {
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.resetCases();
		}
	}

	public void stepForward() {

		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.stepForward();
			for (int j = 0; j < Parameters.demeCount; j++) {
				if (i != j) {
					HostPopulation hpOther = demes.get(j);
					hp.betweenDemeContact(hpOther);
				}
			}
		}

		Parameters.day += Parameters.deltaT;

	}

	public void run() {

		try {

			File seriesFile = new File("out.timeseries");
			seriesFile.delete();
			seriesFile.createNewFile();
			PrintStream seriesStream = new PrintStream(seriesFile);
			System.out.println("day\tdiversity\ttmrca\tnetau\tserialInterval\tantigenicDiversity\tN\tS\tI\tR\tcases");
			printHeader(seriesStream);

			while (Parameters.day < (double) Parameters.endDay) {

				if (Parameters.day % (double) Parameters.printStep < Parameters.deltaT) {
					updateDiversity();
					printState();
					if (Parameters.day > Parameters.burnin) {
						printState(seriesStream);
						pushLists();
					}
					resetCases();
				}

				if (getI()==0) {
					if (Parameters.repeatSim) {
						reset();
						seriesFile.delete();
						seriesFile.createNewFile();
						seriesStream = new PrintStream(seriesFile);
						printHeader(seriesStream);
					} else {
						break;
					}
				}

				stepForward();

			}

			seriesStream.close();

			writeDataCSV();
		} catch(IOException ex) {
			System.out.println("Could not write to file");
			System.exit(0);
		}

		// tree reduction
		VirusTree.pruneTips();
		VirusTree.markTips();
		VirusTree.reroot();

		// tree prep
		makeTrunk();
		VirusTree.fillBackward();
		VirusTree.sortChildrenByDescendants();
		VirusTree.setLayoutByDescendants();
		VirusTree.streamline();

		// rotation
		if (Parameters.pcaSamples) {
			VirusTree.rotate();
			VirusTree.flip();
		}


		// Summary
		printSummary();
		VirusTree.printMKSummary();		// appends to out.summary


		if (!Parameters.reducedOutput) {

			// tip and tree output
			VirusTree.printTips();
			VirusTree.printBranches();
			VirusTree.printNewick();

			// immunity output
			if (Parameters.phenotypeSpace.equals("geometric") || Parameters.phenotypeSpace.equals("geometricSeq")) {
				VirusTree.updateRange();
				VirusTree.printRange();
				if (Parameters.immunityReconstruction) {
					printImmunity();
				}
			}

			// detailed output
			if (Parameters.detailedOutput) {
				printHostPopulation();
			}

		}
	}

	private void writeDataCSV() throws FileNotFoundException {
		// Creates csv file from the most recent out.timeseries (i.e., not from example/out.timeseries)
		Scanner input = new Scanner(new File("out.timeseries"));
		PrintStream output = new PrintStream("out_timeseries.csv");

		// Check for next line
		while(input.hasNextLine()) {
			String line = input.nextLine();
			Scanner lineScan = new Scanner(line);

			// Check for next token
			while(lineScan.hasNext()) {
				String token = lineScan.next();
				if (lineScan.hasNext()) {
					output.print(token + ",");
				} else {
					output.print(token);
				}
			}
			output.println();
		}
	}

	public void reset() {
		Parameters.day = 0;
		diversity = 0;
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.reset();
		}
		VirusTree.clear();
	}

}