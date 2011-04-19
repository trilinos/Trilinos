#!/bin/bash

# Simple script I used to drive checkin-test.py from within the base Trilinos
# source repo, which uses a CHECKIN subdir to store the output files.  This
# makes it easy to do commands like --pull with extra repos to see what the
# status of things are.  To use this script just do:
#
#  $ mkdir CHECKIN
#  $ ln -s sampleScripts/checkin-test-local.sh .
#
# then run it, for example, as:
#
#  $ ./checkin-test-local.sh --extra-repos=preCopyrightTrilinos --push
#

EXTRA_ARGS=$@

cd CHECKIN && ../checkin-test.py $EXTRA_ARGS
