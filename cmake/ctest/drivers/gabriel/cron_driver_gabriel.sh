#!/bin/bash

echo
echo "Starting nightly Trilinos testing on gabriel: `date`"
echo

export TDD_PARALLEL_LEVEL=2
export TDD_HTTP_PROXY="http://wwwproxy.sandia.gov:80/"
time ../cron_driver.py

echo
echo "Ending nightly Trilinos testing on gabriel: `date`"
echo

/home/rabartl/mailmsg.py "Finished nightly Trilinos CMake tests gabriel: http://trilinos-dev.sandia.gov/cdash/index.php?project=Trilinos"
