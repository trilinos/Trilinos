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
-DBUILD_SHARED_LIBS:BOOL=ON
-DTrilinos_ENABLE_Fortran:BOOL=OFF
-DTrilinos_EXTRA_LINK_FLAGS:STRING=\"-Wl,-rpath,/opt/gcc/gcc-4.3.2/lib64\"
" > COMMON.config


echo "
-DMPI_BASE_DIR:PATH=/usr/local/openmpi/gcc-4.3.2/openmpi-1.3.2
-DTrilinos_WARNINGS_AS_ERRORS_FLAGS:STRING=\"\"
" > MPI_DEBUG.config


echo "
-DCMAKE_CXX_COMPILER=/opt/gcc/gcc-4.3.2/bin/g++
-DCMAKE_C_COMPILER=/opt/gcc/gcc-4.3.2/bin/gcc
" > SERIAL_RELEASE.config


#
# Run the checkin-test.py script with more arguments
#

../../Trilinos/checkin-test.py \
--send-email-to=bakercg@ornl.gov \
--make-options=\"-j12\" \
--ctest-options=\"-j2\" \
--ctest-timeout=180 \
--commit-msg-header-file=checkin_message \
--extra-pull-from=\"zan master\" \
--do-all \
$EXTRA_ARGS


# NOTE: The above remote 'gabriel' was created from:
#
#   $ eg remote add gabriel gabriel:~/PROJECTS/Trilinos.base/Trilinos
#
# This allows the shorthand 'gabriel' in lots of eg/git operations.

# NOTE: If you want to pull from a different remote you have to quote
# the --extra-pull-from option to this shell script as:
#
#    "--extra-pull-from=\"otherrepo master\""
#
# Otherwise, the space will mess things up!

# NOTE: Even though godel has 8 cores, I will only use 6 of them so that I
# don't dominate the machine.

# Above, I left in the --commit-msg-header option in case I need to
# fix problems and add more commits.  I like editing a file before I
# commit.
