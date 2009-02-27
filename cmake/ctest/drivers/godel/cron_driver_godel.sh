#!/bin/bash

# Source the 
cd $HOME
source .bash_profile

#CTEST_EXE=/usr/local/bin/ctest
CTEST_EXE=/home/rabartl/install/bin/ctest

echo
echo "Starting nightly Trilinos testing on godel: `date`"
echo

echo
echo "Checking out just the skeleton cmake/ctest code: `date`"
echo

BASEDIR=/home/rabartl/PROJECTS/dashboards/Trilinos.base

cd $BASEDIR
cvs -q -d :ext:software:/space/CVS co Trilinos/cmake Trilinos/CTestConfig.cmake

#echo
#echo "Doing dependency checking-only build: `date`"
#echo
#
#time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/ctest_linux_nightly_package_deps_godel.cmake -VV

echo
echo "Doing serial performance build: `date`"
echo

time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/ctest_linux_nightly_serial_performance_godel.cmake -VV

echo
echo "Doing mpi optimized build: `date`"
echo

time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/ctest_linux_nightly_mpi_optimized_godel.cmake -VV

echo
echo "Doing serial debug build: `date`"
echo

time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/ctest_linux_nightly_serial_debug_godel.cmake -VV

echo
echo "Doing mpi optimized shared library build: `date`"
echo

time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/ctest_linux_nightly_mpi_optimized_shared_godel.cmake -VV

echo
echo "Doing serial optimized implicit instantiation build: `date`"
echo

time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/ctest_linux_nightly_serial_opt_impl_instant_godel.cmake -VV

echo
echo "Doing mpi optimized zoltan c-only build: `date`"
echo

time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/ctest_linux_nightly_mpi_opt_zoltan_c_godel.cmake -VV

echo
echo "Doing serial debug intel build: `date`"
echo

time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/ctest_linux_nightly_serial_debug_icpc_godel.cmake -VV

echo
echo "Doing serial debug memcheck build: `date`"
echo

time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/ctest_linux_nightly_serial_debug_memcheck_godel.cmake -VV

echo
echo "Doing mpi debug memcheck build: `date`"
echo

time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/ctest_linux_nightly_mpi_debug_memcheck_godel.cmake -VV

echo
echo "Ending nightly Trilinos testing on godel: `date`"
echo

/home/rabartl/mailmsg.py "Finished nightly Trilinos tests godel: http://trilinos.sandia.gov/cdash/index.php?project=Trilinos"
