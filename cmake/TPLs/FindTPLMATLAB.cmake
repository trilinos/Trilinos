# TPL finding routine for MATLAB.  A bit more advanced and general than the FindMatlab.cmake included w/ cmake.
#
# Any questions as to how this works? Ask Chris Siefert <csiefer@sandia.gov>


# Check to make sure MATLAB_ARCH is set
IF(NOT MATLAB_ARCH)
  MESSAGE(FATAL_ERROR "You need to set MATLAB_ARCH to use MATLAB with Trilinos (hint: It's probably something like glnx86, glnxa64, mac, maci, maci64, sol2, or sol64).")
ENDIF()

# Check to make sure MATLAB_ROOT is set
IF(NOT MATLAB_ROOT)
  MESSAGE(FATAL_ERROR "You need to set MATLAB_ROOT to use MATLAB with Trilinos.  This should be set to the root of your MATLAB install.")
ENDIF()

# Check to make sure MPI is OFF
IF(TPL_ENABLE_MPI)
  MESSAGE(FATAL_ERROR "MATLAB TPL is incompatible with MPI.  Please disable MPI if you want to use MATLAB.")
ENDIF()

# Add Include/Library directories
# Note #1: I have to add this to the cache otherwise TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES will clobber them completely.
# Note #2: ADVANCED_SET does not handle extra parameters (like CACHE) correctly, hence why I can't use it here
SET(MATLAB_INCLUDE_DIRS ${MATLAB_ROOT}/extern/include CACHE PATH "Include directories for Matlab")
MARK_AS_ADVANCED(MATLAB_INCLUDE_DIRS)
SET(MATLAB_LIBRARY_DIRS ${MATLAB_ROOT}/bin/${MATLAB_ARCH} CACHE PATH "Lib directories for Matlab")
MARK_AS_ADVANCED(MATLAB_LIBRARY_DIRS)

# Make sure we can find the matlab libs
TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( MATLAB
  REQUIRED_HEADERS engine.h mex.h
  REQUIRED_LIBS_NAMES mex mx eng
  )

# Find mex & mexext
SET(MATLAB_MEX_DIR ${MATLAB_ROOT}/bin CACHE PATH "Directory of MATLAB mex compiler wrapper")
FIND_PROGRAM(MEX_COMPILER mex ${MATLAB_MEX_DIR} NO_DEFAULT_PATH)
IF(NOT MEX_COMPILER)
  MESSAGE(FATAL_ERROR " Could not find mex.")
ENDIF()
IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
  MESSAGE("Found mex: " ${MEX_COMPILER})
ENDIF()

FIND_PROGRAM(MEX_MEXEXT mexext ${MATLAB_MEX_DIR})
IF(NOT MEX_MEXEXT)
  MESSAGE(FATAL_ERROR " Could not find mexext.")
ENDIF()
IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
  MESSAGE("Found mexext: " ${MEX_MEXEXT})
ENDIF()

# Find the matlab file exension from mexext
EXECUTE_PROCESS(COMMAND ${MEX_MEXEXT} OUTPUT_VARIABLE MEX_EXTENSION OUTPUT_STRIP_TRAILING_WHITESPACE)
SET(MEX_EXTENSION ${MEX_EXTENSION} CACHE STRING "MATLAB mex file extension")
IF(NOT MEX_EXTENSION)
  MESSAGE(FATAL_ERROR " Platform-specific mex extension could not be found (hint: check to be sure mexext runs correctly).")
ENDIF()
IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
  MESSAGE("Platform-specific mex extension: " ${MEX_EXTENSION})
ENDIF()
 
