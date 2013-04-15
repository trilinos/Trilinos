#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on trilinos-test: `date`"
echo

# TrilinosDriver settings:
#
export TDD_GIT_EXE=/home/trilinos/git/bin/git
export TDD_PARALLEL_LEVEL=2
export TDD_HTTP_PROXY="http://wwwproxy.sandia.gov:80/"

# Trilinos settings:
#
export TDD_CTEST_TEST_TYPE=Nightly

# Machine specific environment:
#
export PATH=/home/trilinos/cmake/bin:/home/trilinos/git/bin:/home/trilinos/tpl/gcc4.1.2/qt-4.5.2/bin:/usr/kerberos/bin:/usr/local/bin:/bin:/usr/bin

# BASEDIR is the parent directory of this script's Trilinos source tree...
BASEDIR=`cd "\`dirname \"$0\"\`/../../../../..";pwd`
echo BASEDIR=$BASEDIR
BASEDATADIR=$BASEDIR

export CVS_RSH=ssh
export LD_LIBRARY_PATH="$BASEDIR/MPI_OPT_DEV_SHARED/BUILD/packages/PyTrilinos/src"
export PYTHONPATH="/home/trilinos/tpl/gcc4.1.2/mpi4py-1.3/lib64/python:/home/trilinos/tpl/gcc4.1.2/numpy1.4.0/lib64/python2.4/site-packages"
export TRILINOSDATADIRECTORY="$BASEDATADIR/TrilinosData"

pushd "$BASEDATADIR"
cvs -q -d :ext:software.sandia.gov:/space/CVS co TrilinosData
popd

# Machine independent cron_driver:
#
SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on trilinos-test: `date`"
echo

