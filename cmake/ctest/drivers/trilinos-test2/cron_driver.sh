#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on trilinos-test2: `date`"
echo

# TrilinosDriver settings:
#
export TDD_GIT_EXE=/home/trilinos/bin/git
export TDD_PARALLEL_LEVEL=2
export TDD_CTEST_TEST_TYPE=Nightly

# Trilinos settings:
#
#export CTEST_TEST_TYPE=Experimental
#export Trilinos_PACKAGES=Teuchos

# Machine specific environment:
#
export PATH=/home/trilinos/tpl/gcc4.1.2/swig-3.0.4/bin:/home/trilinos/python-2.7.9/bin:/home/trilinos/bin:/usr/kerberos/bin:/usr/local/bin:/bin:/usr/bin:/home/jmwille/bin:/home/trilinos/tpl/gcc4.1.2/qt-4.5.2/bin

# BASEDIR is the parent directory of this script's Trilinos source tree...
BASEDIR=`cd "\`dirname \"$0\"\`/../../../../..";pwd`
echo BASEDIR=$BASEDIR
BASEDATADIR=$BASEDIR

export CVS_RSH=ssh
export LD_LIBRARY_PATH="/home/trilinos/lib:$BASEDIR/MPI_OPT_DEV_SHARED/BUILD/packages/PyTrilinos/src:/home/trilinos/gcc4.7.2/base/lib64"
export PYTHONPATH="/home/trilinos/tpl/gcc4.1.2/numpy-1.9.1/lib/python2.7/site-packages"
export TRILINOSDATADIRECTORY="$BASEDATADIR/TrilinosData"

pushd "$BASEDATADIR"
cvs -q -d :ext:software.sandia.gov:/space/CVS co TrilinosData
popd

# Machine independent cron_driver:
#
SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on trilinos-test2: `date`"
echo

