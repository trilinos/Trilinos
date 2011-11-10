#!/usr/bin/env python

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

"""
Script for doing checkin testing of Trilinos.  Please run checkin-test.py -h
for details
"""


#
# Import commands
#


import sys
import os
import traceback

scriptsDir = os.path.abspath(os.path.dirname(sys.argv[0]))+"/cmake/python"
sys.path.insert(0, scriptsDir)

from GeneralScriptSupport import *

#
# Read in the commandline arguments
#

#print "sys.argv:", sys.argv

# Create a deep copy of the commandline arguments
cmndLineArgs = []
cmndLineArgs.extend(sys.argv)

# See if the help option is set or not
helpOpt = len( set(cmndLineArgs) & set(("--help", "-h")) ) > 0

# See if --show-defaults was set or not
showDefaultsOpt = len( set(cmndLineArgs) & set(("--show-defaults", "dummy")) ) > 0

#
# Forward the options but tee the output
#

if (not helpOpt) and (not showDefaultsOpt):
  logFileName = "checkin-test.out"
else:
  logFileName = ""

cmnd = scriptsDir+"/checkin-test-impl.py " + requoteCmndLineArgs(sys.argv[1:])

if logFileName:
  cmnd = cmnd + " 2>&1 | tee "+logFileName

# This return value is always 0 even if it fails?
rtnVal = echoRunSysCmnd(cmnd, throwExcept=False)

# Grep the output to determine success or failure
success = True
if logFileName and getCmndOutput("grep 'REQUESTED ACTIONS: PASSED' "+logFileName, True, False)=="":
  success = False
  
if success:
  sys.exit(0)
else:
  sys.exit(1)
