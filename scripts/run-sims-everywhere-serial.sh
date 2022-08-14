#!/bin/bash
set -eu

jar=$(find $PWD -name "antigen.jar");

# Move to desired ouput directory
cd ../python/results

for YML in $(find $PWD -d -name "*.yml")
do
    cd $(dirname $YML)
    echo "\n*** Starting $YML"
      java -jar $jar -XX:+UseSerialGC -Xmx1G Antigen
    echo "*** Completed $YML"
done
