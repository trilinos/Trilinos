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

# This is a modification of the standard PythonInterp module that adds
# version checking.

# - Find python interpreter
# This module finds if Python interpreter is installed and determines where the
# executables are. This code sets the following variables:
#
#  PYTHONINTERP_FOUND - Was the Python executable found
#  PYTHON_EXECUTABLE  - path to the Python interpreter
#
# This module also responds to the following CMake cache varibles:
#
#  PythonInterp_FIND_VERSION - Requires that a given version be found
#    if some Python executable is found (e.g. '2.4').  Note that this
#    results in the most recent Python version that that is supported
#    on the system, not the minimum version!
#  PythonInterp_MUST_BE_FOUND - If set to TRUE, then the python executable
#    must be found or the CMake configure will stop immediately.
#

#MESSAGE("Calling TriBITS version of FindPythonInterp.cmake")

FIND_PROGRAM(PYTHON_EXECUTABLE
  NAMES python2.7 python2.6 python2.5 python2.4 python2.3 python2.2 python2.1 python2.0 python1.6 python1.5 python
  PATHS
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.7\\InstallPath]
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.6\\InstallPath]
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.5\\InstallPath]
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.4\\InstallPath]
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.3\\InstallPath]
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.2\\InstallPath]
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.1\\InstallPath]
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.0\\InstallPath]
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\1.6\\InstallPath]
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\1.5\\InstallPath]
  )

# handle the QUIETLY and REQUIRED arguments and set PYTHONINTERP_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PythonInterp DEFAULT_MSG PYTHON_EXECUTABLE)

MARK_AS_ADVANCED(PYTHON_EXECUTABLE)

#
# Version checking: If a version check is requested, set
# PythonInterp_VERSION, convert it to a list, extract the last element
# of the list, and compare it to the requested version
#MESSAGE("Initial: PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
#MESSAGE("PythonInterp_FIND_VERSION = ${PythonInterp_FIND_VERSION}")
#MESSAGE("PythonInterp_MUST_BE_FOUND = ${PythonInterp_MUST_BE_FOUND}")
IF (PYTHON_EXECUTABLE OR PythonInterp_MUST_BE_FOUND)
  IF (NOT PYTHON_EXECUTABLE)
    MESSAGE(FATAL_ERROR "Error, Python must be found!")
  ENDIF()
  IF (PythonInterp_FIND_VERSION)
    EXECUTE_PROCESS(COMMAND
      ${PYTHON_EXECUTABLE} -c "import sys; print sys.version.split()[0]"
      OUTPUT_VARIABLE PythonInterp_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    MESSAGE(STATUS "Python version ${PythonInterp_VERSION}")
    #SEPARATE_ARGUMENTS(PythonInterp_VERSION)
    #LIST(GET PythonInterp_VERSION 1 PythonInterp_VERSION)
    IF(${PythonInterp_VERSION} VERSION_LESS ${PythonInterp_FIND_VERSION})
      MESSAGE(WARNING
        "Python version ${PythonInterp_VERSION}"
        " is less than required version ${PythonInterp_FIND_VERSION}!"
        "  Disabling Python!"
        )
      SET(PYTHONINTERP_FOUND FALSE)
      UNSET(PYTHON_EXECUTABLE)
      UNSET(PYTHON_EXECUTABLE CACHE)
    ENDIF()
  ENDIF()
ENDIF()
#MESSAGE("Final: PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
