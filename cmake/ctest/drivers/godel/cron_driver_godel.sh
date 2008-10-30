#!/bin/bash

# Source the 
cd $HOME
source .bash_profile

cd /home/rabartl/PROJECTS/dashboards/Trilinos/SERIAL_DEBUG/Trilinos/cmake/ctest
cvs -q update -dP

ctest -S /home/rabartl/PROJECTS/dashboards/Trilinos/SERIAL_DEBUG/Trilinos/cmake/ctest/ctest_linux_nightly_serial_debug_godel.cmake

$HOME/mailmsg.py "Trilinos serial debug on thumper: http://trilinos.sandia.gov/cdash"
