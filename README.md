# Antigen: Epidemiological Simulation Framework

Antigen is a powerful simulation framework for studying virus evolution and population dynamics in epidemiological contexts. Built for researchers in infectious disease epidemiology, it implements SIR (Susceptible-Infected-Recovered) models with genetic and phenotypic evolution of pathogens.

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
git clone https://github.com/matsengrp/antigen.git
cd antigen

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

### For Epidemiologists & Researchers
- **[Getting Started Guide](docs/user-guide/running-simulations.md)** - Step-by-step simulation workflow
- **[Epidemiological Model](docs/user-guide/epidemiological-model.md)** - Understanding the SIR framework  
- **[Parameter Guide](docs/user-guide/parameters-reference.md)** - Configure simulations for your research
- **[Output Analysis](docs/user-guide/output-analysis.md)** - Interpret and analyze results
- **[Example Simulations](docs/examples/)** - Common research scenarios

### Technical Resources
- **[Installation Details](docs/installation/)** - Compilation, troubleshooting, system requirements
- **[API Reference](docs/api-reference/)** - Class documentation for developers
- **[Development Guide](docs/development/)** - Extending Antigen and contributing

<!-- ## 🔬 Research Applications

Antigen is designed for epidemiological research questions such as:
- **Seasonal influenza dynamics**: Model antigenic drift and vaccine effectiveness
- **Multi-strain competition**: Study how different virus strains compete and evolve
- **Population structure effects**: Analyze how geographic structure affects spread
- **Immune history impact**: Understand how past infections shape future disease risk
- **Vaccination strategies**: Evaluate intervention timing and coverage -->

## 💡 Key Features

### Epidemiological Model
- **SIR dynamics** with recovered immunity based on antigenic similarity
- **Multi-deme populations** with configurable migration rates
- **Seasonal transmission** patterns with region-specific parameters
- **Host demographics** including birth, death, and aging processes

### Virus Evolution
- **Antigenic phenotypes** in continuous n-dimensional space or sequence-based
- **Mutation models** with epitope/non-epitope site differences
- **Phylogenetic tracking** of virus genealogies throughout simulation
- **Cross-immunity** based on antigenic distance between strains

## 📊 Example Output

Running the default simulation generates epidemiological timeseries showing:
- Regional infection prevalence and incidence
- Virus sampling and genetic diversity
- Phylogenetic relationships between circulating strains
- Summary statistics of key epidemiological parameters

## 🏥 Citing Antigen

If you use Antigen in your research, please cite:
```
TBD...
```

## 🤝 Contributing & Support

- **Issues & Questions**: [GitHub Issues](https://github.com/matsengrp/antigen/issues)
- **Feature Requests**: Open an issue with the "enhancement" label
- **Contributing**: See [Development Guide](docs/development/contributing.md)

## 📄 License

Copyright Trevor Bedford 2010-2024. Distributed under the GPL v3.

---

**Quick Links**: [Installation](docs/installation/) | [User Guide](docs/user-guide/) | [Examples](docs/examples/) | [API Reference](docs/api-reference/)