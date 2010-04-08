#!/usr/bin/env python

# Prerequisites:
# - python is version 2.4 or later

import glob
import os
import shutil
import subprocess
import sys

# SCRIPT_DIR is the directory where *this* script is:
#
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# BASE_DIR is the parent directory of our containing "Trilinos" source tree,
# which we compute relative to SCRIPT_DIR:
#
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
  SCRIPT_DIR))))

TOOLS_DIR = BASE_DIR + "/tools"

# Make sure tools directory exists:
#
if not os.path.exists(TOOLS_DIR):
  os.makedirs(TOOLS_DIR)
  if not os.path.exists(TOOLS_DIR):
    print "error: could not create directory \"" + TOOLS_DIR + "\""
    sys.exit(1)

origDir = os.getcwd()
os.chdir(BASE_DIR)

print "SCRIPT_DIR: +" + SCRIPT_DIR + "+"
print "BASE_DIR: +" + BASE_DIR + "+"
print "TOOLS_DIR: +" + TOOLS_DIR + "+"

# Download and install CMake/CTest 'release' build
#
CMAKE_DIR = TOOLS_DIR + "/cmake-TDD"

TDD_CMAKE_INSTALLER_TYPE = "release"
if "TDD_CMAKE_INSTALLER_TYPE" in os.environ:
  if (os.environ["TDD_CMAKE_INSTALLER_TYPE"] != ""):
    TDD_CMAKE_INSTALLER_TYPE = os.environ["TDD_CMAKE_INSTALLER_TYPE"]

TDD_FORCE_CMAKE_INSTALL = 0
if "TDD_FORCE_CMAKE_INSTALL" in os.environ:
  if (os.environ["TDD_FORCE_CMAKE_INSTALL"] != ""):
    TDD_FORCE_CMAKE_INSTALL = os.environ["TDD_FORCE_CMAKE_INSTALL"]

# Only install cmake-TDD if necessary or if requested.
# (Requires network connectivity; avoid when possible.)
#
if not os.path.exists(CMAKE_DIR) or TDD_FORCE_CMAKE_INSTALL != 0:
  print "Downloading and installing CMake to \"" + CMAKE_DIR + "\"..."

  if os.path.exists(CMAKE_DIR):
    shutil.rmtree(CMAKE_DIR)

  cmake_installer_output = (subprocess.Popen([
    sys.executable
    , BASE_DIR + "/Trilinos/cmake/python/download-cmake.py"
    , "--skip-detect"
    , "--install-dir="+CMAKE_DIR
    , "--installer-type="+TDD_CMAKE_INSTALLER_TYPE
    ], stdout=subprocess.PIPE).communicate()[0])

  print cmake_installer_output

else:
  print "Skipped download and install of CMake"

# Find ctest under CMAKE_DIR:
#
gr = glob.glob(CMAKE_DIR + "/bin/ctest*")
if 0 == len(gr):
  gr = glob.glob(CMAKE_DIR + "/*/bin/ctest*")
if 0 == len(gr):
  gr = glob.glob(CMAKE_DIR + "/*/*/bin/ctest*")
if 1 != len(gr):
  print "error: could not find ctest executable after download..."
  os.chdir(origDir)
  sys.exit(2)

CTEST_EXE = gr[0]
print "CTEST_EXE: +" + CTEST_EXE + "+"

if not os.path.exists(CTEST_EXE):
  print "error: ctest does not exist after installation..."
  os.chdir(origDir)
  sys.exit(3)

# Verify ctest works with a simple --version call first:
#
CTEST_VERSION = (subprocess.Popen([
  CTEST_EXE
  , "--version"
  ], stdout=subprocess.PIPE).communicate()[0])
CTEST_VERSION = CTEST_VERSION.strip()
print "CTEST_VERSION: +" + CTEST_VERSION + "+"

# Run one TrilinosDriver dashboard for this Trilinos source tree:
#
print "Running ctest -S TrilinosDriverDashboard.cmake..."

pr = (subprocess.Popen([
  CTEST_EXE
  , "-S"
  , SCRIPT_DIR+"/TrilinosDriverDashboard.cmake"
  ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT))

CTEST_OUTPUT = pr.communicate()[0]
print CTEST_OUTPUT

CTEST_RESULT = pr.returncode
print "CTEST_RESULT: +" + str(CTEST_RESULT) + "+"

if CTEST_RESULT != 0:
  print "error: ctest returned non-zero error value, script will exit with " + str(CTEST_RESULT)

# Propagate ctest return value
#
os.chdir(origDir)
sys.exit(CTEST_RESULT)
