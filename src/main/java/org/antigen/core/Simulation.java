package org.antigen.core;

/* Simulation functions, holds the host population */

import java.util.*;
import java.io.*;

import com.javamex.classmexer.*;
import org.antigen.host.*;
import org.antigen.virus.*;
import org.antigen.phenotype.*;
import org.antigen.analysis.*;

public class Simulation {
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

	// Get average infection risk of a phenotype amongst a given sample size



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

	public void printPopulationImmunityCentroids(PrintStream historyStream, int day) {
		// Write CSV header on first call
		if (day == 0) {
			historyStream.println("year,deme,ag1,ag2,naive_fraction,experienced_hosts");
		}
		
		double year = (double) day / 365.0;
		
		// Calculate and output individual deme centroids
		double globalSumAg1 = 0.0;
		double globalSumAg2 = 0.0;
		int globalExperiencedHosts = 0;
		int globalTotalSamples = 0;
		
		for (int i = 0; i < Parameters.demeCount; i++) {
			int nSamples = Parameters.hostImmunitySamplesPerDeme[i];
			ImmunitySummary summary = demes.get(i).getPopulationImmunitySummary(nSamples);
			
			// Output deme-specific data
			if (summary.hasValidCentroid()) {
				historyStream.printf("%.4f,%s,%.6f,%.6f,%.4f,%d%n", 
					year, 
					Parameters.demeNames[i], 
					summary.getCentroid()[0], 
					summary.getCentroid()[1],
					summary.getNaiveFraction(),
					summary.getExperiencedHosts());
			} else {
				// All hosts are naive in this deme
				historyStream.printf("%.4f,%s,NaN,NaN,%.4f,%d%n", 
					year, 
					Parameters.demeNames[i], 
					summary.getNaiveFraction(),
					summary.getExperiencedHosts());
			}
			
			// Accumulate for global centroid (only if valid)
			if (summary.hasValidCentroid()) {
				globalSumAg1 += summary.getCentroid()[0] * summary.getExperiencedHosts();
				globalSumAg2 += summary.getCentroid()[1] * summary.getExperiencedHosts();
			}
			globalExperiencedHosts += summary.getExperiencedHosts();
			globalTotalSamples += summary.getTotalSampled();
		}
		
		// Calculate and output global centroid
		if (globalExperiencedHosts > 0) {
			double globalAg1 = globalSumAg1 / globalExperiencedHosts;
			double globalAg2 = globalSumAg2 / globalExperiencedHosts;
			double globalNaiveFraction = (double) (globalTotalSamples - globalExperiencedHosts) / globalTotalSamples;
			
			historyStream.printf("%.4f,total,%.6f,%.6f,%.4f,%d%n", 
				year, 
				globalAg1, 
				globalAg2,
				globalNaiveFraction,
				globalExperiencedHosts);
		} else {
			// All sampled hosts are naive globally
			historyStream.printf("%.4f,total,NaN,NaN,1.0000,%d%n", 
				year, 
				globalExperiencedHosts);
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

			File outDirs = new File(Parameters.outPath);
			outDirs.mkdirs();
			File historyFile = new File("out.histories.csv");
			historyFile.delete();
			historyFile.createNewFile();
			PrintStream historyStream = new PrintStream(historyFile);
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

				// print immunity if needed
				if (Parameters.sampleHostImmunity && Parameters.day % (double) Parameters.printHostImmunityStep < Parameters.deltaT) {
					printPopulationImmunityCentroids(historyStream, (int) Parameters.day);
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
			historyStream.close();

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
			System.out.println("Writing tips file...");
			VirusTree.printTips();
			System.out.println("Writing branches file...");
			VirusTree.printBranches();
			System.out.println("Writing FASTA file...");
			VirusTree.printFASTA();
			System.out.println("Writing newick tree file...");
			VirusTree.printNewick();

			// immunity output
			if (Parameters.phenotypeSpace.equals("geometric") || Parameters.phenotypeSpace.equals("geometricSeq")) {
				VirusTree.updateRange();
				VirusTree.printRange();
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