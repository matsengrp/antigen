# antigen-prime

antigen-prime is an epidemiological simulation framework written in Java that models virus evolution and population dynamics. It simulates SIR (Susceptible-Infected-Recovered) models with genetic/phenotypic evolution of pathogens.

## Quick Start

```bash
# Clone and build
git clone https://github.com/matsengrp/antigen-prime.git
cd antigen-prime
mvn clean package

# Run simulation
java -jar target/antigen-prime.jar -XX:+UseSerialGC -Xmx8G
```

## Documentation

- [Requirements](installation/requirements.md) - System requirements and dependencies
- [Compilation](installation/compilation.md) - Build instructions
- [Mutation Model](user-guide/mutation-model.md) - K80 mutation and antigenic evolution
- [Immunity Centroids](user-guide/immunity-centroids.md) - Host population immunity sampling
- [Parameters Reference](user-guide/parameters-reference.md) - Configuration options
- [Output Analysis](user-guide/output-analysis.md) - Understanding simulation output

## Key Features

- **SIR Epidemiology**: Full susceptible-infected-recovered dynamics with demographics
- **Metapopulation Structure**: Multiple demes with migration and seasonal forcing
- **Phenotype Evolution**: Multiple evolution models (geometric, sequence-based, DMS-informed)
- **Phylogenetic Tracking**: Complete genealogy of virus lineages
- **Flexible Sampling**: Configurable sampling schemes for downstream analysis

## Building Documentation

```bash
pip install -r docs/requirements.txt
mkdocs serve  # preview at http://127.0.0.1:8000
mkdocs build  # build static site to site/
```
