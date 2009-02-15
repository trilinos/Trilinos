# This is a modification of the standard PythonInterp module that adds
# version checking.

# - Find python interpreter
# This module finds if Python interpreter is installed and determines where the
# executables are. This code sets the following variables:
#
#  PYTHONINTERP_FOUND - Was the Python executable found
#  PYTHON_EXECUTABLE  - path to the Python interpreter
#

FIND_PROGRAM(PYTHON_EXECUTABLE
  NAMES python2.6 python2.5 python2.4 python2.3 python2.2 python2.1 python2.0 python1.6 python1.5 python
  PATHS
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
IF(PythonInterp_FIND_VERSION)
  EXECUTE_PROCESS(COMMAND
    ${PYTHON_EXECUTABLE} "--version"
    ERROR_VARIABLE PythonInterp_VERSION
    ERROR_STRIP_TRAILING_WHITESPACE
    )
  SEPARATE_ARGUMENTS(PythonInterp_VERSION)
  LIST(GET PythonInterp_VERSION -1 PythonInterp_VERSION)
  IF(${PythonInterp_VERSION} VERSION_LESS ${PythonInterp_FIND_VERSION})
    MESSAGE(FATAL_ERROR
      "Python version " ${PythonInterp_VERSION}
      " is less than required version " ${PythonInterp_FIND_VERSION}
      )
  ENDIF(${PythonInterp_VERSION} VERSION_LESS ${PythonInterp_FIND_VERSION})
ENDIF(PythonInterp_FIND_VERSION)
