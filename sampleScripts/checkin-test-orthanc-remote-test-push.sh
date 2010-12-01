#!/bin/bash

#
# This is the script that I use for remote test/push on
# godel.sandia.gov pulling from gabriel.sandia.gov.  You can copy this
# script and adapt it to your own machine.
#
# To do a remote test/push by pulling from gabriel.sandia.gov, just pass in the options:
#
#   --do-all --push
#
# then just wait for the return email
#

EXTRA_ARGS=$@

source ~/.bashrc

#
# Set up basic environment options
#

echo "
" > COMMON.config

echo "
" > MPI_DEBUG.config


echo "
-DCMAKE_C_COMPILER=icc
-DCMAKE_CXX_COMPILER=icpc
-DCMAKE_Fortran_COMPILER=ifort
" > SERIAL_RELEASE.config


#
# Run the checkin-test.py script with more arguments
#

../../Trilinos/checkin-test.py \
--send-email-to=bakercg@ornl.gov \
--make-options=\"-j4\" \
--ctest-options=\"-j1\" \
--ctest-timeout=180 \
--commit-msg-header-file=checkin_message \
$EXTRA_ARGS


# NOTE: The above remote 'gabriel' was created from:
#
#   $ eg remote add gabriel gabriel:~/PROJECTS/Trilinos.base/Trilinos
#
# This allows the shorthand 'gabriel' in lots of eg/git operations.

# NOTE: Even though godel has 8 cores, I will only use 6 of them so that I
# don't dominate the machine.

# Above, I left in the --commit-msg-header option in case I need to
# fix problems and add more commits.  I like editing a file before I
# commit.
