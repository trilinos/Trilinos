# - Find the python numpy module

# This cmake module determines whether the python module numpy exists
# for the current PYTHON_EXECUTABLE. This code sets the following
# variable:
#
#  NUMPY_FOUND       = Set to aTRUE if numpy is found
#  NUMPY_VERSION     = NumPy version number
#  NUMPY_INCLUDE_DIR = Path to numpy include files
#

SET(NUMPY_FOUND FALSE)

IF(PYTHON_EXECUTABLE)

  EXECUTE_PROCESS(
    COMMAND ${PYTHON_EXECUTABLE} -c 
    "import numpy; print numpy.__version__"
    OUTPUT_VARIABLE NUMPY_VERSION
    ERROR_VARIABLE  NUMPY_VERSION_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  IF(NOT NUMPY_VERSION_ERROR)

    SET(NUMPY_FOUND TRUE)

    EXECUTE_PROCESS(
      COMMAND ${PYTHON_EXECUTABLE} -c 
      "import numpy; print numpy.get_include()"
      OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
      ERROR_VARIABLE  NUMPY_INCLUDE_ERROR
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    IF(NUMPY_INCLUDE_ERROR)

      EXECUTE_PROCESS(
        COMMAND ${PYTHON_EXECUTABLE} -c 
        "import numpy; print numpy.get_numpy_include()"
        OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )

    ENDIF(NUMPY_INCLUDE_ERROR)

  ENDIF(NOT NUMPY_VERSION_ERROR)

ENDIF(PYTHON_EXECUTABLE)
