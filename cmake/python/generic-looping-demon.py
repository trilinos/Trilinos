#!/bin/env python

# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER


# TODO: Put a time limit on the command so that it will not run past the
# RUNTILL time.  You could do this by running a cmake -P driver script that
# calls EXECUTE_PROCESS(...).  There is also some python code that exists out
# there that could also implement a time limit but it is complex.


#
# Read in the command-line arguments
#

usageHelp = r"""generic-looping-demon.py [OPTIONS]

This simple program takes a command as input and runs it over and over again
(pausing in-between iterations for a given time) and then stops at the given
abolute time.

The reason that the script takes an absolute time instead of a relative time
is that this script is desiged to drive continuous itegration (CI) processes
where the CI process should shut down at some point.

NOTE: The last iteration will not start later than --run-till=RUNTILL but the
last command could run until significantly after that time.
"""

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--command", dest="command", type="string", default="",
  help="The shell command that will get run in the loop" )

clp.add_option(
  "--loop-interval", dest="loopInterval", type="string", default="",
  help="Input to the standard unix/linux 'sleep' command" \
  +" (e.g. '60s') to space out iterations of the script.")

dateFormat = "%Y-%m-%d"
timeFormat = "%H:%M:%S"
dateTimeFormat = dateFormat+" "+timeFormat
clp.add_option(
  "--run-till", dest="runTill", type="string", default="",
  help="The absolute time the script will run iterations till." \
  +" This takes the format "+dateTimeFormat+" (e.g. 2011-06-09 18:00:00).")

clp.add_option(
  "--today-run-till", dest="todayRunTill", type="string", default="",
  help="The time today the script will run iterations till." \
  +" This takes the format "+timeFormat+" (e.g. 18:00:00).")

(options, args) = clp.parse_args()


#
# Check the input arguments
#

import sys

if not options.command:
  print "\nError, you must set the --command argument!"
  sys.exit(1)

if not options.loopInterval:
  print "\nError, you must set the --loop-interval argument!"
  sys.exit(1)

if not (options.runTill or options.todayRunTill):
  print "\nError, you must set either the --run-till or --today-run-till argument!"
  sys.exit(1)


#
# Echo the command-line
#

print ""
print "**************************************************************************"
print "Script: generic-looping-demon.py \\"

print "  --command='"+options.command+"' \\"
print "  --loop-interval='"+options.loopInterval+"' \\"
if options.runTill:
  print "  --run-till='"+options.runTill+"' \\"
if options.todayRunTill:
  print "  --run-till='"+options.todayRunTill+"' \\"


#
# Helper functions
#

# Parses "%Y-%m-%d %H:%M:%S" into a datetime.datetime object
def parseDateTimeString(dateTimeStr):
  (dateStr, timeStr) = dateTimeStr.split(" ")
  (year, month, day) = dateStr.split("-")
  (hour, minute, second) = timeStr.split(":")
  return datetime.datetime(int(year), int(month), int(day),
    int(hour), int(minute), int(second))


def formatDateTime(dateTimeObj):
  return datetime.datetime.strftime(dateTimeObj, dateTimeFormat)


#
# Executable statements
#

import os
import datetime

scriptsDir = os.path.abspath(os.path.dirname(sys.argv[0]))+"/cmake/python"
sys.path.insert(0, scriptsDir)

from GeneralScriptSupport import *

if options.runTill:
  finalDateTime = parseDateTimeString(options.runTill)
elif options.todayRunTill:
  todayDate = datetime.datetime.today()
  todayDateStr = datetime.datetime.strftime(todayDate, dateFormat)
  finalDateTime = parseDateTimeString(todayDateStr+" "+options.todayRunTill)

print "\nThe script will run iterations till = " + formatDateTime(finalDateTime) + "\n"

currentTime = datetime.datetime.now()
iteration = 0

while currentTime < finalDateTime:
  print "*****************************************************************************\n" \
    + str(iteration) + ":" \
    + " current time = " + formatDateTime(currentTime) \
    + ", final time = " + formatDateTime(finalDateTime)
  echoRunSysCmnd(options.command, throwExcept=False, timeCmnd=True)
  echoRunSysCmnd("sleep "+options.loopInterval)
  currentTime = datetime.datetime.now()
  iteration += 1
