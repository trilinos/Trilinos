#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on trilinos-test: `date`"
echo

#
# Is mpd running
#

mpdhost=$(mpdtrace)

if [ "$mpdhost" != "s861036" ] ; then
  echo "Starting mpd"
  mpd &
fi

# TrilinosDriver settings:
#
export TDD_GIT_EXE=/usr/local/bin/git
export TDD_PARALLEL_LEVEL=2
export TDD_HTTP_PROXY="http://wwwproxy.sandia.gov:80/"
export TDD_CTEST_TEST_TYPE=Nightly

#
# just testing
#
export TDD_IN_TESTING_MODE=ON
export CTEST_DO_SUBMIT=FALSE
export CTEST_DO_UPDATES=FALSE



# Machine specific environment:
#
export PATH=.:/usr/local/bin:/usr/kerberos/bin:/usr/bin:/bin:/usr/X11R6/bin:/sbin:/home/lriesen/bin

# BASEDIR is the parent directory of this script's Trilinos source tree...
BASEDIR=`cd "\`dirname \"$0\"\`/../../../../..";pwd`
echo BASEDIR=$BASEDIR
BASEDATADIR=$BASEDIR

export CVS_RSH=ssh
export LD_LIBRARY_PATH="$BASEDIR/MPI_OPT_DEV_SHARED/BUILD/packages/PyTrilinos/src"
export PYTHONPATH="/usr/local/lib/python2.7:/usr/local/lib/python2.7/site-packages:/usr/local/lib/python2.7/lib-dynload"

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

