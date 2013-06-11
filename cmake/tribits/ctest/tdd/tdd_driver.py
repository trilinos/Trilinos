#!/usr/bin/env python

# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
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
# ************************************************************************
# @HEADER

# Prerequisites:
# - python is version 2.4 or later

import glob
import os
import shutil
import sys

# See if we should be verbose and print everything!
if "TDD_DEBUG_VERBOSE" in os.environ:
  verboseEnvValue = os.environ["TDD_DEBUG_VERBOSE"]
  if verboseEnvValue == "1":
    verbose = True
  else:
    verbose = False
else:
  verbose = False


# tribitsTddDriverDir is the directory where *this* script is:
#
this_path = os.path.abspath(os.path.realpath(__file__))
tribitsTddDriverDir = os.path.dirname(this_path)
print "tribitsTddDriverDir = '"+tribitsTddDriverDir+"'"


# Load the general script support python code
sys.path.insert(0, tribitsTddDriverDir+"/../../python")
from GeneralScriptSupport import *


def invoke_ctest(ctestExe, script, tddDashboardRootDir, environment = {}):
  """
  Invokes CTest using the executable given by the ctestExe argument,
  the script file specified by the script argument, in the working
  directory specified by the tddDashboardRootDir argument and the set of
  environment variables specified in the environment map.
  """
  # We have to pass parameters to CTest through environment
  # variables. It would be nice to have another way to do this, but
  # until ctest supports something like CMake's -D argument, this is
  # how it has to be done.
  ctest_environ = None
  if environment:
    print "environment =", environment
    # Use a dictionary that is distinct from the os environment dictionary.
    ctest_environ = os.environ.copy()
    # Append additional environment variables to the
    ctest_environ.update(environment)

  cmd = ctestExe
  if verbose:
    cmd = ctestExe + " -VV"
  ctestRtn = echoRunSysCmnd(
    cmd + " -S" + " " + script,
    throwExcept=False,
    timeCmnd=True,
    workingDir=tddDashboardRootDir,
    environment = ctest_environ
    )

  print "ctestRtn: '" + str(ctestRtn) + "'"
  
  if ctestRtn != 0:
    print "error: ctest returned non-zero error value, script will exit with " + str(ctestRtn)
    
  # Propagate ctest return value
  #
  return ctestRtn


def run_driver(ctestSourceDirectory, projectRepoBaseDir):
  """
  Run the dashboard driver. The ctestSourceDirectory argument specifies
  where the directory that CTest will run over. There should be a
  CMakeLists.txt file in this location. The projectRepoBaseDir argument
  specifies root of the source code repository for the project.
  """
  origDir = os.getcwd()
  try:
    print "\n******************************************************************"
    print "***        Tribits Driver Dashboard tdd_driver.py              ***" 
    print "******************************************************************\n"

    print "\nPWD=\""+os.getcwd()+"\"...\n"
    print "projectRepoBaseDir = '" + projectRepoBaseDir + "'"
    print "tribitsTddDriverDir = '" + tribitsTddDriverDir + "'"

    # tribitsDir is the root directory of the TriBITS system:
    #
    tribitsDir = os.path.dirname(os.path.dirname(tribitsTddDriverDir))
    print "tribitsDir = '"+tribitsDir+"'"

    # dashboardBaseDir is the parent directory of our containing source tree,
    # which we compute relative to tribitsDir:
    #
    tddDashboardRootDir = os.path.dirname(projectRepoBaseDir)
    if "TDD_DASHBOARD_ROOT" in os.environ:
      tddDashboardRootDir = os.environ["TDD_DASHBOARD_ROOT"]
    print "tddDashboardRootDir = '"+tddDashboardRootDir+"'"

    # dashboardToolsDir is the directory to which any needed tools will be downloaded.
    #
    dashboardToolsDir = tddDashboardRootDir + "/tools"
    print "dashboardToolsDir = '"+dashboardToolsDir+"'"

    # Make sure tools directory exists:
    #
    if not os.path.exists(dashboardToolsDir):
      os.makedirs(dashboardToolsDir)
      if not os.path.exists(dashboardToolsDir):
        print "error: could not create directory \"" + dashboardToolsDir + "\""
        sys.exit(1)

    os.chdir(tddDashboardRootDir)
    if verbose: "\nNew PWD = '"+os.getcwd()+"'"

    # Download and install CMake/CTest to use for the outer driver
    #
    cmakeTddDownloadBaseDir = dashboardToolsDir + "/cmake-TDD"

    TDD_CMAKE_INSTALLER_TYPE = "release"
    if "TDD_CMAKE_INSTALLER_TYPE" in os.environ:
      TDD_CMAKE_INSTALLER_TYPE = os.environ["TDD_CMAKE_INSTALLER_TYPE"]

    TDD_FORCE_CMAKE_INSTALL = "1"
    if "TDD_FORCE_CMAKE_INSTALL" in os.environ:
      TDD_FORCE_CMAKE_INSTALL = os.environ["TDD_FORCE_CMAKE_INSTALL"]

    TDD_HTTP_PROXY = ""
    if "TDD_HTTP_PROXY" in os.environ:
      TDD_HTTP_PROXY = os.environ["TDD_HTTP_PROXY"]

    # Only install cmake-TDD if necessary or if requested.
    # (Requires network connectivity; avoid when possible.)
    #

    print "\n***"
    print "*** Downloading and installing CMake to \"" + cmakeTddDownloadBaseDir + "\"..."
    print "***\n"

    installMasterCMake = False
    if not os.path.exists(cmakeTddDownloadBaseDir):
      print "Forcing install of master CMake because '"+cmakeTddDownloadBaseDir+"' does not exist!"
      installMasterCMake = True
    elif TDD_FORCE_CMAKE_INSTALL == "1":
      print "Forcing install of master CMake because" \
        + " TDD_FORCE_CMAKE_INSTALL == 1!"
      installMasterCMake = True
    else:
      print "Leaving current CMake in place ..." \

    if installMasterCMake:

      cmnd =  sys.executable + " " \
        + tribitsDir + "/python/download-cmake.py" \
        + " --skip-detect" \
        + " --install-dir="+cmakeTddDownloadBaseDir \
        + " --installer-type="+TDD_CMAKE_INSTALLER_TYPE

      if TDD_HTTP_PROXY:
        cmnd += " --http-proxy="+TDD_HTTP_PROXY

      try:
        echoRunSysCmnd( cmnd,
          timeCmnd = True,
          workingDir = dashboardToolsDir \
          )
      except Exception, e:
        print "WARNING! The following command failed!\n"+cmnd
        print "However, not updating CMake is not the end of the world!"


    # Find ctest under cmakeTddDownloadBaseDir:
    #
    ctestGlobStr = glob.glob(cmakeTddDownloadBaseDir + "/bin/ctest*")
    if 0 == len(ctestGlobStr):
      ctestGlobStr = glob.glob(cmakeTddDownloadBaseDir + "/*/bin/ctest*")
    if 0 == len(ctestGlobStr):
      ctestGlobStr = glob.glob(cmakeTddDownloadBaseDir + "/*/*/bin/ctest*")
    if 1 != len(ctestGlobStr):
      print "error: could not find ctest executable after download..."
      sys.exit(2)

    ctestExe = ctestGlobStr[0]
    print "\nctestExe = '" + ctestExe + "'"

    if not os.path.exists(ctestExe):
      print "error: ctest does not exist after installation..."
      sys.exit(3)

    # Escape any spaces in the path of the ctest exe. This has to be done
    # here instead of where we set the ctestExe the first time because
    # the check for existence cannot handle the "\"
    #
    ctestExe = ctestExe.replace(" ",  "\ ")

    # Verify ctest works with a simple --version call first:
    #

    ctestVersion = getCmndOutput(ctestExe+" --version", True, False)
    print "ctestVersion = '"+ctestVersion+"'"

    # Run one driver dashboard for this source tree:
    #

    print "\n***"
    print "*** Running the main dashboards as CTest tests .."
    print "***\n"
    sys.exit(
      invoke_ctest(ctestExe,
        os.path.join(tribitsTddDriverDir, "TribitsDriverDashboard.cmake"),
        tddDashboardRootDir,
        {
         "TDD_DASHBOARD_ROOT" : tddDashboardRootDir,
         "CTEST_SOURCE_DIRECTORY" : ctestSourceDirectory,
         "CTEST_UPDATE_DIRECTORY" : projectRepoBaseDir,
         "CTEST_BINARY_DIRECTORY" : tddDashboardRootDir+"/TDD_BUILD",
         }
        )
      )
    
  finally:
    os.chdir(origDir)
