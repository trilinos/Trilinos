# - Find the python include path

# This module determines where the python include file Python.h is,
# based on the current value of Python3_EXECUTABLE. This code sets the
# following variable:
#
#  PYTHON_INCLUDE_PATH  = path to where Python.h is found
#

IF(Python3_EXECUTABLE)
  # Obtain the candidate path for python include
  EXECUTE_PROCESS(COMMAND
    ${Python3_EXECUTABLE} -c 
    "import sys; print(sys.prefix + '/include/python' + sys.version[:3])"
    OUTPUT_VARIABLE CANDIDATE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  # Verify theat Python.h exists within the candidate path
  IF(EXISTS "${CANDIDATE}")
    IF(EXISTS ${CANDIDATE}/Python.h)
      SET(PYTHON_INCLUDE_FOUND TRUE)
      SET(PYTHON_INCLUDE_PATH ${CANDIDATE})
    ENDIF(EXISTS ${CANDIDATE}/Python.h)
  ENDIF(EXISTS "${CANDIDATE}")
  # Obtain the candidate path for python library
  EXECUTE_PROCESS(COMMAND
    ${Python3_EXECUTABLE} -c 
    "import sys; print('python' + sys.version[:3])"
    OUTPUT_VARIABLE PYVERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  EXECUTE_PROCESS(COMMAND
    ${Python3_EXECUTABLE} -c 
    "import sys; print(sys.prefix + '/lib/python' + sys.version[:3] + '/config')"
    OUTPUT_VARIABLE CANDIDATE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  # Verify theat Python.h exists within the candidate path
  IF(EXISTS "${CANDIDATE}")
    FOREACH(SUFFIX dll dylib so a)
      SET(LIB_CANDIDATE ${CANDIDATE}/lib${PYVERSION}.${SUFFIX})
      IF(EXISTS ${LIB_CANDIDATE})
	SET(PYTHON_LIB_FOUND TRUE)
	SET(PYTHON_LIBRARIES ${LIB_CANDIDATE})
	BREAK()
      ENDIF(EXISTS ${LIB_CANDIDATE})
    ENDFOREACH(SUFFIX)
  ENDIF(EXISTS "${CANDIDATE}")
ENDIF(Python3_EXECUTABLE)

# Set any python variables that did not get set by the above logic
IF(PYTHON_INCLUDE_FOUND AND PYTHON_LIB_FOUND)
  # We're good
ELSE(PYTHON_INCLUDE_FOUND AND PYTHON_LIB_FOUND)
  FIND_PACKAGE(PythonLibs)
ENDIF(PYTHON_INCLUDE_FOUND AND PYTHON_LIB_FOUND)
