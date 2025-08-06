## Introduction

Antigen implements an SIR epidemiological model where hosts in a population are infected with
viruses that have distinct antigenic phenotypes.  Hosts make contacts transmitting viruses and also
recover from infection.  After recovery, a host remembers the antigenic phenotype it was infected
with as part of its immune history.  The risk of infection after contact depends on comparing the
infecting virus's phenotype to the phenotypes in the host immune history.

Antigenic phenotype is currently implemented as a *n*-d continuous vector.  Virus mutations move
phenotype randomly in this *n*-d Euclidean space.  However, other phenotype models may be
implemented through the Phenotype interface.

Additionally, population structure is implemented in terms of discrete demes.  Contacts within a
deme occur through standard mass action, while contacts between demes occur at some fraction of the
rate of within deme contact.

-------------------------------------------

## Running

### Setting up Java on Fred Hutch cluster


To load Java 11, run these commands:

	source /app/lmod/lmod/init/profile
	ml Java/11;



The program can be compiled with Maven:

	mvn clean compile

To run tests:

	mvn test

A transportable jar file with all dependencies can be created with:

	mvn package

This creates two JAR files:
- `target/antigen-prime.jar` - Complete executable with all dependencies (recommended)
- `target/antigen-prime-no-dependencies.jar` - Classes only, requires classpath setup
	
Then to run from the main jar:

	java -jar target/antigen-prime.jar -XX:+UseSerialGC -Xmx1G

For development, you can compile and run directly:

	mvn compile exec:java -Dexec.mainClass="org.antigen.Antigen" -Dexec.args="-XX:+UseSerialGC -Xmx32G"
	
This requires Java 16 or higher to compile and run. The `-Xmx1G` option is used to increase memory allocation. 
This may need to be increased further with larger host population sizes. The `-XX:+UseSerialGC` option 
swaps the default Java garbage collector to something that works much more efficiently for Antigen.
	
## Parameters
	
Parameter defaults can be seen in [`src/main/java/org/antigen/core/Parameters.java`](src/main/java/org/antigen/core/Parameters.java).  When run, the program looks 
for the file [`src/main/resources/parameters.yml`](src/main/resources/parameters.yml) and dynamically loads these in, overwriting defaults.

## Output

The simulation will output a timeseries of region-specific prevalence and incidence to
[`out.timeseries`](example/out.timeseries).  It will also sample viruses periodically and output their 
geographic and antigenic locations to [`out.tips`](example/out.tips) and a tree connecting these samples 
to [`out.branches`](example/out.branches).  This file contains pairs of viruses, child and parent, 
representing nodes in a genealogy.  Average values are output to [`out.summary`](example/out.summary).

If you have Mathematica, you can generate a number of figures from this output by running the
notebook [`antigen-analysis.nb`](example/antigen-analysis.nb).

## Memory

The `-Xmx1G` is required, because as an individual-based model, the memory requirements are
typically quite large. Each host requires a minimum of 40 bytes of memory, plus 8 bytes per
Phenotype recorded in its immune history.  If the yearly attack rate is 15% and the host life span
is 30 years, at equilibrium the average size of the immune history will be 4.5 references.  This
gives memory usage of: population size x 76 bytes.  With 7.5 million hosts (used in the default
parameters), the equals 570MB.

In addition to hosts and immune histories, the simulation tracks the virus genealogy through
VirusTree.  This is harder to profile, and will continually grow in memory usage throughout the
simulation.  With the default parameters, VirusTree takes 5.5 MB at the end of a simulated year and
may up to 110 MB at the end of the default 20 simulated years.

Memory can be easily profiled by calling `jmap -histo <PID>`.

## Manual Compilation (Alternative)

If Maven is not available, the program can be compiled manually with:

	javac -cp "src/main/java:src/main/resources" src/main/java/org/antigen/*.java src/main/java/org/antigen/*/*.java

Or if you need to specify JAR dependencies explicitly:

	javac -classpath "src/main/java:src/main/resources:~/.m2/repository/colt/colt/1.2.0/colt-1.2.0.jar:~/.m2/repository/org/yaml/snakeyaml/2.0/snakeyaml-2.0.jar" src/main/java/org/antigen/*.java src/main/java/org/antigen/*/*.java

Then to run:

	java -cp "src/main/java:src/main/resources:~/.m2/repository/colt/colt/1.2.0/colt-1.2.0.jar:~/.m2/repository/org/yaml/snakeyaml/2.0/snakeyaml-2.0.jar" -XX:+UseSerialGC -Xmx32G org.antigen.Antigen

Note: Maven is the recommended build method as it handles all dependencies automatically.

-------------------------------------------

Copyright Trevor Bedford 2010-2014. Distributed under the GPL v3.
