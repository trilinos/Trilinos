#!/bin/bash

#
# Grab just the top lines of each top output to see what processes are taking
# the most.  This script should be run on the output from the monitor-top.sh
# command.  Therefore, if monitor-top.sh was run as:
#
#   $ monitor-top.sh &> SOME_BASE_DIR/top.log
#
# then you would run this script as:
#
#   $ summarize-monitor-top.sh SOME_BASE_DIR/top.log [NUM_LINES_TO_TAKE]
#

# The top log file
TOP_LOG_INPUT_FILE=$1

# Number of line to take.  Note, this will include the 'top' header so you
# will get several less of the actual process lines.
NUMBER_OF_LINES_TO_TAKE=$2

if [ "$NUMBER_OF_LINES_TO_TAKE" == "" ] ; then
  NUMBER_OF_LINES_TO_TAKE=15
fi

grep -A $NUMBER_OF_LINES_TO_TAKE '^top - ' $TOP_LOG_INPUT_FILE
