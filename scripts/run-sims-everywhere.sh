#!/bin/bash
set -eu
SCRIPT_PATH="/path/to/script.sh"

# Move to desired ouput directory
cd ../python/results

for YML in $(find $PWD -d -name "*.yml")
do
    cd $(dirname $YML)
    echo "\n*** Starting $YML"
      java -jar ../../../antigen.jar -XX:+UseSerialGC -Xmx1G Antigen
    echo "*** Completed $YML"
done
