#!/bin/bash

#
# This is the script that I used to checkin to Trilinos on gabriel.sandia.gov.
# You can copy this script and adapt it to your own machine.
#

#
# Allow command-line arguments to pass through to cmake configure!
#

EXTRA_ARGS=$@

#
# Set up some specific options that we will not want to change through the
# command-line
#

echo "-DBUILD_SHARED_LIBS:BOOL=ON" > COMMON.config

# 2009/08/28: Sundance is not building for this build case.  Kevin long said
# to turn it off for now.  Plese don't do this in your scripts unless you
# absolutely have to!
echo "-DTrilinos_ENABLE_Sundance:BOOL=OFF" > SERIAL_RELEASE.config

#
# Run the standard checkin testing script with my specializations
#

/home/rabartl/PROJECTS/Trilinos.base/Trilinos/cmake/python/checkin-test.py \
--make-options="-j4" \
--ctest-options="-j4" \
--commit-msg-header-file=checkin_message \
$EXTRA_ARGS  


#
# NOTES:
#
# (*) Enabling shared libaries makes relinks go *much* faster and massively
# speeds up the checkin testing process when rebulding!
#
# (*) Be sure to set --make-options="-jN" to speed up building.  It makes a
# *big* difference.
#
# (*) Passing -jN to ctest with --ctest-optioins="-jN" can speed up running
# the tests but I have not seen very good speedup in general and some people
# have reported no speedup at all.  You should experiment with ctest -jN to
# see what number N works well on your machine.
#
