#!/bin/bash

#
# This is the script that I used to checkin to Trilinos on my Windows
# Vista laptop using Cygwin.
#

EXTRA_ARGS=$@

#
# Set up basic environment options
#


# I have to disable MPI in the MPI_DEBUG build (which becomes a
# SERIAL_DEBUG build) because I don't have MPI on this system.
# However, this is better than nothing.
echo "-DTPL_ENABLE_MPI:BOOL=OFF" > COMMON.config

# I also have to turn off warnings as errors because there is a
# warning coming from the STL headers with gcc 3.4.4 which comes
# standard on Cygwin.
echo "-DTrilinos_WARNINGS_AS_ERRORS_FLAGS:STRING=" >> COMMON.config


#
# Run the checkin-test.py script with more arguments
#

/cygdrive/c/_mystuff/PROJECTS/Trilinos.base/Trilinos/checkin-test.py \
--make-options="-j2" \
--ctest-options="-j2" \
--ctest-timeout=180 \
--send-email-to= \
--commit-msg-header-file=checkin_message \
$EXTRA_ARGS  

# NOTE: Above I have to turn off sending email because there is no
# such thing as mail or mailx on cygwin.
