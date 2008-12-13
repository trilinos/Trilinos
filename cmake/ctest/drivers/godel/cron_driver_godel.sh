#!/bin/bash

# Source the 
cd $HOME
source .bash_profile

CTEST_EXE=/usr/local/bin/ctest

echo
echo "Starting time: `date`"
echo

echo
echo "A) Checking out just the drivers"
echo

BASEDIR=/home/rabartl/PROJECTS/dashboards/Trilinos

cd $BASEDIR
cvs -q -d :ext:@software.sandia.gov:/space/CVS co -d scripts Trilinos/cmake/ctest

echo
echo "Doing dependency checking-only build"
echo

time ${CTEST_EXE} -S $BASEDIR/scripts/ctest_linux_nightly_package_deps_godel.cmake -VV

/home/rabartl/mailmsg.py "Trilinos dependency checking finished on godel: http://trilinos.sandia.gov/cdash/index.php?project=Trilinos"

echo
echo "Doing serial performance build"
echo

time ${CTEST_EXE} -S $BASEDIR/scripts/ctest_linux_nightly_serial_performance_godel.cmake -VV

/home/rabartl/mailmsg.py "Trilinos serial performance finished on godel: http://trilinos.sandia.gov/cdash/index.php?project=Trilinos"

echo
echo "Doing mpi optimized build"
echo

time ${CTEST_EXE} -S $BASEDIR/scripts/ctest_linux_nightly_mpi_optimized_godel.cmake -VV

/home/rabartl/mailmsg.py "Trilinos mpi opt finished on godel: http://trilinos.sandia.gov/cdash/index.php?project=Trilinos"

echo
echo "Doing serial debug build"
echo

time ${CTEST_EXE} -S $BASEDIR/scripts/ctest_linux_nightly_serial_debug_godel.cmake -VV

/home/rabartl/mailmsg.py "Trilinos serial debug finished on godel: http://trilinos.sandia.gov/cdash/index.php?project=Trilinos"

echo
echo "Doing mpi optimized shared library build"
echo

time ${CTEST_EXE} -S $BASEDIR/scripts/ctest_linux_nightly_mpi_optimized_shared_godel.cmake -VV

/home/rabartl/mailmsg.py "Trilinos mpi opt shared finished on godel: http://trilinos.sandia.gov/cdash/index.php?project=Trilinos"

echo
echo "Ending time: `date`"
echo
