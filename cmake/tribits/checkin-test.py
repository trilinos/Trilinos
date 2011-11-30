#!/usr/bin/env python

# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Test, and Integrate System
#                 Copyright (2011) Sandia Corporation
#
#
# Copyright (2011) Sandia Corporation. Under the terms of Contract
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
Script for doing checkin testing of a TriBITS built project.  Please run
checkin-test.py -h for details
"""


#
# Import commands
#


import sys
import os
import traceback

from optparse import OptionParser

# Get the location of the scripts directory whether from a sym link or the
# actual

checkinTestFilePath = os.path.abspath(sys.argv[0])
checkinTestFileRealPath = os.path.realpath(checkinTestFilePath)
scriptsDir = os.path.dirname(checkinTestFileRealPath)+"/python"
print "scriptsDir='"+scriptsDir+"'"
sys.path.insert(0, scriptsDir)

from GeneralScriptSupport import *
from CheckinTestImpl import *

class WritableTee(object):
  """
  Object that directs all calls to its write method to stdout as well
  as a file. This is to be used as a simple replacement for the Unix
  tee command.
  """
  def __init__(self, outputfile):
    """ Constructor takes a file-like object to write output to."""
    self._realstdout = sys.stdout
    self._outputfile = outputfile

  def _safe_outputfile_method(self, methodname, *args):
    """
    Calls the method specified by methodname with the given args on
    the internal file object if it is non-null.
    """
    if self._outputfile is not None:
      if hasattr(self._outputfile, methodname):
        method = getattr(self._outputfile, methodname)
        if method and callable(method):
          method(*args)

  def write(self, data):
    """
    Write the given data to stdout and to the log file.
    """
    self._realstdout.write(data)
    self._safe_outputfile_method('write', data)

  def flush(self):
    """
    Flush the internal file buffers.
    """
    self._realstdout.flush()
    self._safe_outputfile_method('flush')


def main():
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

  if (not helpOpt) and (not showDefaultsOpt):
    logFile = file("checkin-test.out", "w")
  else:
    logFile = None

  # There are a lot of print statements in the implementation. It's
  # easier to reset sys.stdout and sys.stderr to our WritableTee object
  # than to replace them.
  teeOutput = WritableTee(logFile)
  originalStdout = sys.stdout
  originalStderr = sys.stderr
  try:
    sys.stdout = teeOutput
    sys.stderr = teeOutput
    success = runProjectTestsWithCommandLineArgs(sys.argv[1:])
  except Exception:
    success = False
    traceback.print_exc(file=teeOutput)
  finally:
    # Reset stdout and stderr
    sys.stdout = originalStdout
    sys.stderr = originalStderr

  if success:
    return 0
  else:
    return 1

if __name__ == '__main__':
  sys.exit(main())
