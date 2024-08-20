#!/usr/bin/env python

# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
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


# tribitsDDDir is the directory where *this* script is:
#
this_path = os.path.abspath(os.path.realpath(__file__))
tribitsDDDir = os.path.dirname(this_path)
print "tribitsDDDir = '"+tribitsDDDir+"'"

# Load the general script support python code
pythonUtilsDir = os.path.join(tribitsDDDir, "../python_utils")
sys.path = [pythonUtilsDir] + sys.path
from GeneralScriptSupport import *


def install_ctest(tddDashboardRootDir, tribitsDir):

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
      + tribitsDir + "/python_utils/download-cmake.py" \
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

  return ctestExe


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
  if environment:
    print "environment =", environment

  cmd = ctestExe
  if verbose:
    cmd = ctestExe + " -VV"
  ctestRtn = echoRunSysCmnd(
    cmd + " -S" + " " + script,
    throwExcept=False,
    timeCmnd=True,
    workingDir=tddDashboardRootDir,
    extraEnv = environment
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
    print "tribitsDDDir = '" + tribitsDDDir + "'"

    # tribitsDir is the root directory of the TriBITS system:
    #
    tribitsDir = os.path.abspath(os.path.join(tribitsDDDir, ".."))
    print "tribitsDir = '"+tribitsDir+"'"

    # dashboardBaseDir is the parent directory of our containing source tree,
    # which we compute relative to tribitsDir:
    #
    tddDashboardRootDir = os.path.dirname(projectRepoBaseDir)
    if "TDD_DASHBOARD_ROOT" in os.environ:
      tddDashboardRootDir = os.environ["TDD_DASHBOARD_ROOT"]
    print "tddDashboardRootDir = '"+tddDashboardRootDir+"'"

    os.chdir(tddDashboardRootDir)
    if verbose: "\nNew PWD = '"+os.getcwd()+"'"

    tddUseSystemCTest = False
    if "TRIBITS_TDD_USE_SYSTEM_CTEST" in os.environ \
      and os.environ["TRIBITS_TDD_USE_SYSTEM_CTEST"] == "1" \
      :
      tddUseSystemCTest = True
    print "tddUseSystemCTest =", tddUseSystemCTest

    if tddUseSystemCTest:
      ctestExe = getCmndOutput("which ctest", True, False)
    else:
      ctestExe = install_ctest(tddDashboardRootDir, tribitsDir)

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
        os.path.join(tribitsDDDir, "TribitsDriverDashboard.cmake"),
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
