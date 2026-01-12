# antigen-prime: Epidemiological Simulation Framework

Antigen implements an SIR epidemiological model where hosts in a population are infected with
viruses that have distinct antigenic phenotypes.  Hosts make contacts transmitting viruses and also
recover from infection.  After recovery, a host remembers the antigenic phenotype it was infected
with as part of its immune history.  The risk of infection after contact depends on comparing the
infecting virus's phenotype to the phenotypes in the host immune history.

## 🦠 What Antigen Does

- **Models virus evolution**: Simulate how viruses evolve over time with antigenic drift and immune escape
- **Population dynamics**: Track infections across structured populations with multiple demes/regions
- **Immune history**: Model host immunity based on past infection history and cross-immunity
- **Phylogenetic tracking**: Generate virus genealogies and sample phylogenetic trees
- **Flexible phenotypes**: Support multiple phenotype models (geometric, sequence-based)

## 🚀 Quick Start

### Prerequisites
- Java 16 or higher
- Maven (recommended) or manual Java compilation

### Installation & Basic Run
```bash
# Clone and navigate to directory
git clone https://github.com/matsengrp/antigen-prime.git
cd antigen-prime

# Compile with Maven
mvn clean compile package

# Run simulation (basic)
java -jar target/antigen-prime.jar

# Run with more memory for larger simulations
java -jar target/antigen-prime.jar -XX:+UseSerialGC -Xmx8G
```

### Your First Simulation
The simulation will run with default parameters and create output files in the `output/` directory:
- `out.timeseries` - Epidemic curves and prevalence over time
- `out.tips` - Sampled virus genetic and antigenic data  
- `out.branches` - Phylogenetic relationships between samples
- `out.summary` - Summary statistics

## 📚 Documentation

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

## 🏥 Citing Antigen

If you use Antigen in your research, please cite:
```
TBD...
```

## 🤝 Contributing & Support

- **Issues & Questions**: [GitHub Issues](https://github.com/matsengrp/antigen-prime/issues)
- **Feature Requests**: Open an issue with the "enhancement" label
- **Contributing**: See [Development Guide](docs/development/contributing.md)

## 📄 License

Copyright Trevor Bedford 2010-2024. Distributed under the GPL v3.

---


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
