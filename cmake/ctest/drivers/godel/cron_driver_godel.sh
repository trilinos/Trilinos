#!/bin/bash

# Source the 
cd $HOME
source .bash_profile

echo
echo "A) Checking out just the drivers"
echo

BASEDIR=/home/rabartl/PROJECTS/dashboards/Trilinos

cd $BASEDIR
cvs -q -d :ext:@software.sandia.gov:/space/CVS co -d scripts Trilinos/cmake/ctest

echo
echo "Doing dependency checking-only build"
echo

time /usr/local/bin/ctest -S $BASEDIR/scripts/ctest_linux_nightly_package_deps_godel.cmake -VV

/home/rabartl/mailmsg.py "Trilinos dependency checking on godel: http://trilinos.sandia.gov/cdash/index.php?project=Trilinos"

echo
echo "Doing serial debug build"
echo

time /usr/local/bin/ctest -S $BASEDIR/scripts/ctest_linux_nightly_serial_debug_godel.cmake -VV

/home/rabartl/mailmsg.py "Trilinos serial debug on godel: http://trilinos.sandia.gov/cdash/index.php?project=Trilinos"

echo
echo "Doing mpi optimized build"
echo

time /usr/local/bin/ctest -S $BASEDIR/scripts/ctest_linux_nightly_mpi_optimized_godel.cmake -VV

/home/rabartl/mailmsg.py "Trilinos mpi opt on godel: http://trilinos.sandia.gov/cdash/index.php?project=Trilinos"

echo
echo "Doing serial performance build"
echo

time /usr/local/bin/ctest -S $BASEDIR/scripts/ctest_linux_nightly_serial_performance_godel.cmake -VV

/home/rabartl/mailmsg.py "Trilinos serial performance on godel: http://trilinos.sandia.gov/cdash/index.php?project=Trilinos"
