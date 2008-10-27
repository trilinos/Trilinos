#!/bin/bash

cd $HOME
source .bash_profile

ctest -S /scratch/rabartl/dashboards/Trilinos/cmake/ctest/ctest_linux_nightly_serial_debug_thumper.cmake

$HOME/mailmsg.py "Trilinos serial debug on thumper: http://trilinos.sandia.gov/cdash/index.php?project=Trilinos"
