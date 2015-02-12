#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on standard rhel6-x86_64: `date`"
echo

# TrilinosDriver settings:
#
#update this after git is installed on the NFS mount
export TDD_GIT_EXE=/usr/bin/git
export TDD_PARALLEL_LEVEL=2
export TDD_HTTP_PROXY="http://wwwproxy.sandia.gov:80/"

# Trilinos settings:
#
export TDD_CTEST_TEST_TYPE=Nightly

# Machine specific environment:
#
#export PATH=/usr/local/bin:/bin:/usr/bin/$PATH:/home/trilinos/tpl/gcc4.4.4/qt-4.5.2/bin

# BASEDIR is the parent directory of this script's Trilinos source tree...
BASEDIR=`cd "\`dirname \"$0\"\`/../../../../..";pwd`
echo BASEDIR=$BASEDIR
BASEDATADIR=$BASEDIR

export CVS_RSH=ssh

#export LD_LIBRARY_PATH=/home/trilinos/compilers/gcc/support_libs/mpc-1.0.1/lib:/home/trilinos/compilers/gcc/support_libs/mpfr-3.1.2/lib:/home/trilinos/compilers/gcc/support_libs/gmp-5.1.1/lib:/home/trilinos/compilers/gcc/4.7.2/lib64:$LD_LIBRARY_PATH
#export PYTHONPATH="/home/trilinos/tpl/gcc4.4.4/mpi4py-1.3/lib64/python2.6:/home/trilinos/tpl/gcc4.4.4/numpy-1.6.2/lib64/python2.6/site-packages/"
export TRILINOSDATADIRECTORY="$BASEDATADIR/TrilinosData"

pushd "$BASEDATADIR"
cvs -q -d :ext:software.sandia.gov:/space/CVS co TrilinosData
popd

# Machine independent cron_driver:
#
SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on standard rhel6-x86_64: `date`"
echo

