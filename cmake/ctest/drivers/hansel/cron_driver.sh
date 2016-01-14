#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on hansel: `date`"
echo

# TrilinosDriver settings:
#
export TDD_GIT_EXE=/home/trilinos/git/bin/git
export TDD_PARALLEL_LEVEL=2
export TDD_HTTP_PROXY="http://wwwproxy.sandia.gov:80/"

# Trilinos settings:
#
#export TDD_CTEST_TEST_TYPE=Nightly

# Machine specific environment:
#
#export PATH=/home/trilinos/gcc4.7.2/base/bin:/home/trilinos/cmake/bin:/home/trilinos/git/bin:/home/trilinos/tpl/gcc4.1.2/qt-4.5.2/bin:/usr/kerberos/bin:/usr/local/bin:/bin:/usr/bin

# BASEDIR is the parent directory of this script's Trilinos source tree...
BASEDIR=`cd "\`dirname \"$0\"\`/../../../../..";pwd`
echo BASEDIR=$BASEDIR

# Machine independent cron_driver:
#
SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on hansel: `date`"
echo

