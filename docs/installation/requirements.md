# System Requirements

## Java Requirements

antigen-prime requires **Java 16 or higher** to compile and run. This is due to specific language features and performance optimizations used in the codebase.

### Checking Your Java Version

```bash
java -version
```

You should see something like:
```
openjdk version "16.0.1" 2021-04-20
```

### Installing Java

#### macOS
```bash
# Using Homebrew
brew install openjdk@17

# Using SDKMAN
sdk install java 17.0.1-open
```

#### Linux (Ubuntu/Debian)
```bash
sudo apt update
sudo apt install openjdk-17-jdk
```

#### Linux (CentOS/RHEL)
```bash
sudo yum install java-17-openjdk-devel
```

#### Fred Hutch Cluster (Specific Environment - Skip if Not Applicable)
```bash
source /app/lmod/lmod/init/profile
ml Java/11
```
*Note: This is only for users working on the Fred Hutch computational cluster. 
Most users can skip this section and install Java through the standard methods above. Java 11 is the minimum for this specific cluster, but Java 16+ is still recommended for optimal performance.*

## Build Tools

### Maven (Recommended)
Maven handles all dependencies automatically and is the preferred build method.

#### Installation
- **macOS**: `brew install maven`
- **Linux**: `sudo apt install maven` (Ubuntu) or `sudo yum install maven` (CentOS)
- **Windows**: Download from [Apache Maven](https://maven.apache.org/download.cgi)

### Alternative: Manual Compilation
If Maven is unavailable, manual compilation is supported but requires managing dependencies manually.

## Memory Requirements

### Minimum System Memory
- **Small simulations** (< 100,000 hosts): 2GB RAM
- **Medium simulations** (100,000 - 1M hosts): 8GB RAM  
- **Large simulations** (> 1M hosts): 16GB+ RAM

### JVM Memory Allocation
The `-Xmx` flag controls Java memory allocation:
```bash
# Small simulation
java -jar antigen-prime.jar -Xmx2G

# Large simulation  
java -jar antigen-prime.jar -Xmx16G
```

### Memory Usage Calculation
Each host requires approximately:
- **40 bytes** base memory
- **8 bytes per immune history entry**

With default parameters (15% yearly attack rate, 30-year lifespan):
- Average immune history: 4.5 entries per host
- Memory per host: 40 + (4.5 × 8) = 76 bytes
- For 7.5M hosts: ~570MB

Additional memory is needed for virus genealogy tracking, which grows throughout the simulation.

## Storage Requirements

### Input Files
- Minimal: Configuration files only (~1KB)
- With sequence data: Up to 100MB for large genomes and DMS data

### Output Files
- **Small simulations**: 10-100MB
- **Medium simulations**: 100MB-1GB
- **Large simulations**: 1GB-10GB+

Output size depends on:
- Simulation duration
- Sampling frequency
- Number of demes
- Whether detailed output is enabled

## Performance Considerations

### CPU Requirements
- **Minimum**: Single-core 2GHz processor
- **Recommended**: Multi-core processor (antigen-prime can utilize multiple cores for certain operations)

### Garbage Collection
For best performance, use the serial garbage collector:
```bash
java -jar antigen-prime.jar -XX:+UseSerialGC -Xmx8G
```

This is much more efficient for antigen-prime's memory allocation patterns than the default parallel collector.

## Development Tools (Optional)

For developers and contributors:
- **Git**: For version control and contributing
- **IDE**: IntelliJ IDEA, Eclipse, or VS Code with Java support
- **Testing**: JUnit (included in Maven dependencies)