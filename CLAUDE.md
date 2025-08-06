# CLAUDE.md - Project Context for Antigen

## Project Overview
Antigen is an epidemiological simulation framework written in Java that models virus evolution and population dynamics. The project simulates SIR (Susceptible-Infected-Recovered) models with genetic/phenotypic evolution of pathogens.

## Development Workflow

When implementing new features or changes:
1. **Draft Phase**: Create a markdown file to outline the implementation plan
2. **Issue Creation**: Once the plan is finalized, create a GitHub issue using the markdown content
3. **Branch Naming**: Create a new branch following the convention: `<issue-number>-<brief-feature-description>`
   - Example: `42-add-new-visualization`
   - Example: `15-fix-data-pipeline-bug`

This workflow ensures proper documentation and tracking of all changes.

## Maven Project Structure
```
/ (root directory)
├── pom.xml - Maven configuration
├── src/
│   ├── main/
│   │   ├── java/org/antigen/
│   │   │   ├── Antigen.java - Main entry point
│   │   │   ├── core/
│   │   │   │   ├── Simulation.java - Main simulation logic
│   │   │   │   ├── Parameters.java - Configuration parameters
│   │   │   │   └── Random.java - Random number generation
│   │   │   ├── host/
│   │   │   │   ├── Host.java - Individual host modeling
│   │   │   │   └── HostPopulation.java - Population dynamics
│   │   │   ├── virus/
│   │   │   │   ├── Virus.java - Virus entities
│   │   │   │   ├── VirusTree.java - Phylogenetic tree tracking
│   │   │   │   └── Biology.java - Biological constants and utilities
│   │   │   ├── phenotype/
│   │   │   │   ├── Phenotype.java - Interface for phenotype models
│   │   │   │   ├── PhenotypeFactory.java - Creates phenotype instances
│   │   │   │   ├── GeometricPhenotype.java
│   │   │   │   ├── GeometricPhenotype3D.java
│   │   │   │   ├── GeometricPhenotype10D.java
│   │   │   │   └── GeometricSeqPhenotype.java
│   │   │   └── analysis/
│   │   │       └── SimplePCA.java - PCA analysis utilities
│   │   └── resources/
│   │       ├── parameters.yml - Main parameter file
│   │       ├── codon_table.txt - Genetic code reference
│   │       ├── input/ - Additional parameter files
│   │       └── lib/
│   │           └── classmexer.jar - Memory profiling library
│   └── test/
│       └── java/org/antigen/phenotype/
│           └── TestGeometricSeqPhenotype.java - JUnit tests
├── target/ - Maven build output
├── output/, example/ - Simulation output directories
└── Python scripts - Analysis utilities
```

## Build and Run Commands
```bash
# Compile and package using Maven
mvn clean compile

# Run tests
mvn test

# Create JAR with dependencies
mvn package

# Run simulation (from root directory)
java -jar target/antigen-prime-jar-with-dependencies.jar

# With memory allocation
java -Xmx10G -jar target/antigen-prime-jar-with-dependencies.jar

# Run with specific parameters file
java -Xmx10G -jar target/antigen-prime-jar-with-dependencies.jar parameters.yml

# For development - compile and run directly
mvn compile exec:java -Dexec.mainClass="org.antigen.Antigen"
```

## Key Configuration
- `src/main/resources/parameters.yml` - Main parameter file (YAML format)
- `src/main/resources/codon_table.txt` - Genetic code reference
- `src/main/resources/input/` - Additional parameter files
- `pom.xml` - Maven build configuration and dependencies

## Output Files Generated
- `*.trees` - Phylogenetic tree structure
- `*.tips` - Tree tip information
- `*.branches` - Branch information
- `*.immunity` - Population immunity states
- `*.histories` - Host infection histories
- `*.timeseries` - Epidemiological dynamics
- `out.summary` - Summary statistics

## External Dependencies (Maven managed)
- SnakeYAML - YAML parsing (managed by Maven)
- CERN Colt library - Scientific computing (managed by Maven) 
- JUnit - Unit testing framework (managed by Maven)
- Classmexer - Memory profiling (bundled in `src/main/resources/lib/` and installed as local Maven dependency)

## Python Scripts
- `clustering.py` - Analysis scripts
- `output_csv.py` - Convert output to CSV
- `setup.py` - Setup utilities

## Project Status
- ✅ Proper Maven directory structure implemented
- ✅ Source and compiled files separated
- ✅ Maven build system configured
- ✅ Dependencies managed via Maven
- ✅ Clear package organization by functionality
- ✅ JUnit test framework integrated

## Documentation Status
- Currently restructuring documentation (see `documentation-restructuring-plan.md`)
- Main documentation in `README.md`
- Additional description in `description.md`

## Notes
- Project uses phenotype abstraction for different evolution models
- Memory management is critical for large simulations
- JUnit test framework configured - tests located in `src/test/java/`
- Java 16 compatible for deployment on remote systems