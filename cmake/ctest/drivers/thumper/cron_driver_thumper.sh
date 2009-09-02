#!/bin/bash

# Source the 
cd $HOME
source .bash_profile

cd /scratch/rabartl/dashboards/Trilinos/cmake/ctest
cvs -q update -dP

ctest -S /scratch/rabartl/dashboards/Trilinos/cmake/ctest/ctest_linux_nightly_serial_debug_thumper.cmake

$HOME/mailmsg.py "Trilinos serial debug on thumper: http://trilinos.sandia.gov/cdash/index.php?project=Trilinos"
