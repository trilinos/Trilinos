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

# Prerequisites:
# - python is version 2.4 or later

import glob
import os
import shutil
import sys

# SCRIPT_DIR is the directory where *this* script is:
#
this_path = os.path.abspath(os.path.realpath(__file__))
SCRIPT_DIR = os.path.dirname(this_path)
sys.path.insert(0, SCRIPT_DIR+"/../../python")
from GeneralScriptSupport import *

def invoke_ctest(ctest_exe, script, working_dir, environment = {}, verbose=False):
  """
  Invokes CTest using the executable given by the ctest_exe argument,
  the script file specified by the script argument, in the working
  directory specified by the working_dir argument and the set of
  environment variables specified in the environment map.
  """
  # We have to pass parameters to CTest through environment
  # variables. It would be nice to have another way to do this, but
  # until ctest supports something like CMake's -D argument, this is
  # how it has to be done.
  ctest_environ = None
  if environment:
    # Use a dictionary that is distinct from the os environment dictionary.
    ctest_environ = os.environ.copy()
    # Append additional environment variables to the
    ctest_environ.update(environment)

  cmd = ctest_exe
  if verbose:
    cmd = ctest_exe + " -VV"
  CTEST_RESULT = echoRunSysCmnd(
    cmd + " -S" + " " + script,
    throwExcept=False,
    timeCmnd=True,
    workingDir=working_dir,
    environment = ctest_environ
    )

  print "CTEST_RESULT: +" + str(CTEST_RESULT) + "+"
  
  if CTEST_RESULT != 0:
    print "error: ctest returned non-zero error value, script will exit with " + str(CTEST_RESULT)
    
  # Propagate ctest return value
  #
  return CTEST_RESULT


def run_driver(ctest_directory, repo_directory, ctest_is_verbose=False):
  """
  Run the dashboard driver. The ctest_directory argument specifies
  where the directory that CTest will run over. There should be a
  CMakeLists.txt file in this location. The repo_directory argument
  specifies root of the source code repository for the project.
  """
  origDir = os.getcwd()
  try:
    print "\n*********************************************************"
    print "***        TrilinosDriver cron_driver.py              ***" 
    print "*********************************************************\n"

    print "\nPWD=\""+os.getcwd()+"\"...\n"

    # TRIBITS_DIR is the root directory of the TriBITS system:
    #
    TRIBITS_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))

    # REPO_DIR is the root directory of the source code repository.
    #
    REPO_DIR = repo_directory

    # BASE_DIR is the parent directory of our containing "Trilinos" source tree,
    # which we compute relative to TRIBITS_DIR:
    #
    BASE_DIR = os.path.dirname(REPO_DIR)

    # TOOLS_DIR is the directory to which any needed tools will be downloaded.
    #
    TOOLS_DIR = BASE_DIR + "/tools"

    # Make sure tools directory exists:
    #
    if not os.path.exists(TOOLS_DIR):
      os.makedirs(TOOLS_DIR)
      if not os.path.exists(TOOLS_DIR):
        print "error: could not create directory \"" + TOOLS_DIR + "\""
        sys.exit(1)

    os.chdir(BASE_DIR)

    print "SCRIPT_DIR: +" + SCRIPT_DIR + "+"
    print "TRIBITS_DIR: +" + TRIBITS_DIR + "+"
    print "BASE_DIR: +" + BASE_DIR + "+"
    print "TOOLS_DIR: +" + TOOLS_DIR + "+"
    print "REPO_DIR: +" + REPO_DIR + "+"

    # Download and install CMake/CTest 'release' build
    #
    CMAKE_DIR = TOOLS_DIR + "/cmake-TDD"

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
    print "*** Downloading and installing CMake to \"" + CMAKE_DIR + "\"..."
    print "***\n"

    installMasterCMake = False
    if not os.path.exists(CMAKE_DIR):
      print "Forcing install of master CMake because '"+CMAKE_DIR+"' does not exist!"
      installMasterCMake = True
    elif TDD_FORCE_CMAKE_INSTALL != "0":
      print "Forcing install of master CMake because" \
        + " TDD_FORCE_CMAKE_INSTALL="+TDD_FORCE_CMAKE_INSTALL+" != 0!"
      installMasterCMake = True
    else:
      print "Leaving current CMake in place ..." \

    if installMasterCMake:

      cmnd =  sys.executable + " " \
        + TRIBITS_DIR + "/python/download-cmake.py" \
        + " --skip-detect" \
        + " --install-dir="+CMAKE_DIR \
        + " --installer-type="+TDD_CMAKE_INSTALLER_TYPE

      if TDD_HTTP_PROXY:
        cmnd += " --http-proxy="+TDD_HTTP_PROXY

      try:
        echoRunSysCmnd( cmnd,
          timeCmnd = True,
          workingDir = TOOLS_DIR \
          )
      except Exception, e:
        print "WARNING! The following command failed!\n"+cmnd
        print "However, not updating CMake is not the end of the world!"


    # Find ctest under CMAKE_DIR:
    #
    gr = glob.glob(CMAKE_DIR + "/bin/ctest*")
    if 0 == len(gr):
      gr = glob.glob(CMAKE_DIR + "/*/bin/ctest*")
    if 0 == len(gr):
      gr = glob.glob(CMAKE_DIR + "/*/*/bin/ctest*")
    if 1 != len(gr):
      print "error: could not find ctest executable after download..."
      sys.exit(2)

    CTEST_EXE = gr[0]
    print "\nCTEST_EXE: +" + CTEST_EXE + "+"

    if not os.path.exists(CTEST_EXE):
      print "error: ctest does not exist after installation..."
      sys.exit(3)

    # Escape any spaces in the path of the ctest exe. This has to be done
    # here instead of where we set the CTEST_EXE the first time because
    # the check for existence cannot handle the "\"
    #
    CTEST_EXE = CTEST_EXE.replace(" ",  "\ ")

    # Verify ctest works with a simple --version call first:
    #

    CTEST_VERSION = getCmndOutput(CTEST_EXE+" --version", True, False)
    print "CTEST_VERSION: +" + CTEST_VERSION + "+"

    # Run one TrilinosDriver dashboard for this Trilinos source tree:
    #

    print "\n***"
    print "*** Running the main dashboards as CTest tests .."
    print "***\n"
    sys.exit(invoke_ctest(CTEST_EXE,
                           os.path.join(SCRIPT_DIR, "TribitsDriverDashboard.cmake"),
                           BASE_DIR,
                           {"CTEST_SOURCE_DIRECTORY": ctest_directory}))
    
  finally:
    os.chdir(origDir)

