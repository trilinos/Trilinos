# - Find the python include path

# This module determines where the python include file Python.h is,
# based on the current value of PYTHON_EXECUTABLE. This code sets the
# following variable:
#
#  PYTHON_INCLUDE_PATH  = path to where Python.h is found
#

IF(PYTHON_EXECUTABLE)

  EXECUTE_PROCESS(
    COMMAND ${PYTHON_EXECUTABLE} -c 
    "import sys, os.path; print os.path.join(sys.prefix, 'include', 'python' + sys.version[:3])"
    OUTPUT_VARIABLE PYTHON_INCLUDE_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

ELSE()

  SET(PYTHON_INCLUDE_PATH-NOTFOUND TRUE)

ENDIF()