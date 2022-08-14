#!/bin/bash
set -eu

########### Setup variables ###########
jar=$(find $PWD -name "antigen.jar");
PYTHON_DIRS_SCRIPT=$(find $PWD -name "setup_dirs.py");
IN_FILE=$(find $PWD -name $1 | head -n 1);
EDIT_FILE=$(find $PWD -name $2);

########### Create directory hierarchy ###########
mkdir -p simulations;
cd simulations;

echo "LOG: Running python script to create directory hierachy.";
python $PYTHON_DIRS_SCRIPT $IN_FILE $EDIT_FILE;
echo "LOG: Simulation directory hierarchy creation done.";

########### Run simulations ###########
echo "LOG: Starting simulation...";

for YML in $(find $PWD -d -name "*.yml")
do
    cd $(dirname $YML)
    echo "\n*** Starting $YML"
      java -jar $jar -XX:+UseSerialGC -Xmx1G Antigen
    echo "*** Completed $YML"
done

echo "LOG: Simulations complete!";
