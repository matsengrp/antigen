# Compilation Guide

## Maven Compilation (Recommended)

Maven is the recommended build method as it automatically handles all dependencies and provides a standardized build process.

### Basic Compilation

```bash
# Clean and compile
mvn clean compile

# Run tests
mvn test

# Create JAR file with all dependencies
mvn package
```

### Output Files

Maven creates two JAR files in the `target/` directory:

- **`target/antigen-prime.jar`** - Complete executable with all dependencies (recommended)
- **`target/antigen-prime-no-dependencies.jar`** - Classes only, requires classpath setup

### Running the Application

#### Using the Complete JAR (Recommended)
```bash
java -jar target/antigen-prime.jar -XX:+UseSerialGC -Xmx8G
```

#### Using Maven Exec Plugin (Development)
```bash
mvn compile exec:java -Dexec.mainClass="org.antigen.Antigen" -Dexec.args="-XX:+UseSerialGC -Xmx8G"
```

### Common Maven Commands

```bash
# Clean build artifacts
mvn clean

# Compile without running tests
mvn compile -DskipTests

# Run only specific test
mvn test -Dtest=TestGeometricSeqPhenotype

# Generate project documentation
mvn javadoc:javadoc

# Check for dependency updates
mvn versions:display-dependency-updates
```

## Manual Compilation (Alternative)

If Maven is not available, you can compile manually. However, this requires managing dependencies manually and is not recommended for most users.

### Prerequisites for Manual Compilation

You'll need to have the following JAR files available (typically in `~/.m2/repository/` after running Maven once):

- `snakeyaml-2.0.jar` - YAML parsing
- `colt-1.2.0.jar` - Scientific computing library
- `classmexer.jar` - Memory profiling (optional)

### Manual Compilation Commands

#### Basic Compilation
```bash
javac -cp "src/main/java:src/main/resources" \
  src/main/java/org/antigen/*.java \
  src/main/java/org/antigen/*/*.java
```

#### With Explicit Dependencies
```bash
javac -classpath "src/main/java:src/main/resources:~/.m2/repository/colt/colt/1.2.0/colt-1.2.0.jar:~/.m2/repository/org/yaml/snakeyaml/2.0/snakeyaml-2.0.jar" \
  src/main/java/org/antigen/*.java \
  src/main/java/org/antigen/*/*.java
```

### Running Manually Compiled Code

```bash
java -cp "src/main/java:src/main/resources:~/.m2/repository/colt/colt/1.2.0/colt-1.2.0.jar:~/.m2/repository/org/yaml/snakeyaml/2.0/snakeyaml-2.0.jar" \
  -XX:+UseSerialGC -Xmx8G \
  org.antigen.Antigen
```

## JVM Options Explained

### Memory Management
- **`-Xmx8G`**: Sets maximum heap size to 8GB (adjust based on your simulation size)
- **`-XX:+UseSerialGC`**: Uses serial garbage collector (more efficient for antigen-prime's memory patterns)

### Performance Tuning
```bash
# For large simulations
java -jar target/antigen-prime.jar -XX:+UseSerialGC -Xmx16G -server

# For memory-constrained systems
java -jar target/antigen-prime.jar -XX:+UseSerialGC -Xmx2G -Xms1G
```

### Memory Profiling (Optional)
If you have `classmexer.jar` available:
```bash
java -javaagent:classmexer.jar -jar target/antigen-prime.jar -XX:+UseSerialGC -Xmx8G
```

## Troubleshooting Compilation

### Common Issues

#### "JAVA_HOME not set"
```bash
# macOS
export JAVA_HOME=$(/usr/libexec/java_home)

# Linux
export JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64

# Make permanent by adding to ~/.bashrc or ~/.bash_profile
echo 'export JAVA_HOME=$(/usr/libexec/java_home)' >> ~/.bash_profile
```

#### "Maven command not found"
Install Maven using your system's package manager:
```bash
# macOS
brew install maven

# Ubuntu/Debian
sudo apt install maven

# CentOS/RHEL
sudo yum install maven
```

#### "Out of memory" during compilation
Increase Maven's memory:
```bash
export MAVEN_OPTS="-Xmx2G"
mvn clean package
```

#### Dependency download failures
Clear Maven cache and retry:
```bash
rm -rf ~/.m2/repository
mvn clean package
```

## Build Verification

After successful compilation, verify the build:

```bash
# Check JAR file exists
ls -la target/antigen-prime.jar

# Quick test run (should show help/parameter info)
java -jar target/antigen-prime.jar --help

# Memory test
java -XX:+UseSerialGC -Xmx1G -jar target/antigen-prime.jar
```

A successful build will create a functional JAR file that can run simulations according to the configured parameters.