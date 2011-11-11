#!/bin/bash

# Run top in batch mode (-b) to update every 15 minutes ('-b 900' seconds)
# repreated -n 98 times so it will cover an entire day (i.e. 24*60/15=100)
#
# You can set up a cron entry to log top by doing:
#
#   00 00  *  *  *   SRC_BASE_DIR/Trilinos/cmake/ctest/tools/monitor-top.sh &> SOME_BASE_DIR/top.out
#
# You can then summarize the results with the script summarize-monitor-top.sh

top -b -d 900 -n 98
