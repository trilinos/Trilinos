# - Find the python numpy module

# This cmake module determines whether the python module numpy exists
# for the current PYTHON_EXECUTABLE. This code sets the following
# variable:
#
#  NumPy_FOUND       = Set to TRUE if numpy is found
#  NumPy_VERSION     = NumPy version number
#  NumPy_INCLUDE_DIR = Path to numpy include files
#

#
# If NumPy is required and python executable does not exist, then send
# an error
IF(NOT PYTHON_EXECUTABLE)
  IF(NumPy_FIND_REQUIRED)
    MESSAGE(SEND_ERROR
      "Python executable not found, so required NumPy module not found"
      )
  ENDIF(NumPy_FIND_REQUIRED)
#
# Continue processing if python executable is known
ELSE(NOT PYTHON_EXECUTABLE)

  # Retrieve the NumPy version
  EXECUTE_PROCESS(COMMAND
    ${PYTHON_EXECUTABLE} -c "import numpy; print numpy.__version__"
    OUTPUT_VARIABLE NumPy_VERSION
    ERROR_VARIABLE  NumPy_VERSION_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  # If there is no error retrieving the error, then we know NumPy
  # exists; now look for the NumPy include directory
  IF(NOT NumPy_VERSION_ERROR)
    EXECUTE_PROCESS(COMMAND
      ${PYTHON_EXECUTABLE} -c "import numpy; print numpy.get_include()"
      OUTPUT_VARIABLE NumPy_INCLUDE_DIR
      ERROR_VARIABLE  NumPy_INCLUDE_ERROR
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )

    # If there is a NumPy include error, it is because NumPy is older
    # and the wrong function name was called
    IF(NumPy_INCLUDE_ERROR)
      EXECUTE_PROCESS(COMMAND
	${PYTHON_EXECUTABLE} -c "import numpy; print numpy.get_numpy_include()"
        OUTPUT_VARIABLE NumPy_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
	)
    ENDIF(NumPy_INCLUDE_ERROR)

    # Handle the QUIETLY and REQUIRED arguments and set NumPy_FOUND to
    # TRUE if all listed variables are TRUE
    INCLUDE(FindPackageHandleStandardArgs)
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(NumPy DEFAULT_MSG
      NumPy_VERSION NumPy_INCLUDE_DIR)

    #
    # Version checking: If a version check is requested, compare
    # NumPy_VERSION to the requested version
    IF(NumPy_FIND_VERSION)
      IF(${NumPy_VERSION} VERSION_LESS ${NumPy_FIND_VERSION})
	MESSAGE(FATAL_ERROR
	  "NumPy version " ${NumPy_VERSION}
	  " is less than required version " ${NumPy_FIND_VERSION}
	  )
      ENDIF(${NumPy_VERSION} VERSION_LESS ${NumPy_FIND_VERSION})
    ENDIF(NumPy_FIND_VERSION)

  #
  # A NumPy version error means that NumPy was not found
  ELSE(NOT NumPy_VERSION_ERROR)
    IF(NumPy_FIND_REQUIRED)
      MESSAGE(SEND_ERROR
	"Required NumPy python module not found"
	)
    ENDIF(NumPy_FIND_REQUIRED)

  ENDIF(NOT NumPy_VERSION_ERROR)

ENDIF(NOT PYTHON_EXECUTABLE)
