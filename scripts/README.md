# Running multiple antigen simulations

## Running simulations in a loop with `run-sims-serial.sh`.

This first script automates running many simulations in antigen with different simulation parameters, and storing the results in seperate directories.
This script takes two yaml files as positional arguments and stores simulation outputs from antigen in a directory called `simulations`.

### Positional Arguments
- Input file: path of file holding orginal parameters or "base" paramaters.
- Edits file: path to a yaml file containing the name of the simulation paramaters you want to change.
This tool will create a cartesian product of the parameters specified in this file, so only list what you want to change.

Below is an examples of an `edits.yml` file:
```
endDay: [100, 150]
meanStep: [0.1, 0.3]
```

Run from main directory: `./scripts/run-sims-serial.sh parameters.yml edits.yml`
