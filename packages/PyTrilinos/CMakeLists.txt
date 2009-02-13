
INCLUDE(PackageMacros)

PACKAGE(PyTrilinos)

SET(${PACKAGE_NAME}_VERSION "4.2d")

# Python interpreter and header files are required
FIND_PACKAGE(PythonInterp  REQUIRED)
FIND_PACKAGE(PythonInclude REQUIRED)
MESSAGE("   Python executable   ... " ${PYTHON_EXECUTABLE})
MESSAGE("   Python include path ... " ${PYTHON_INCLUDE_PATH})

# Python module numpy is required
FIND_PACKAGE(NumPy REQUIRED)
MESSAGE("   NumPy version       ... " ${NUMPY_VERSION})
MESSAGE("   NumPy include dir   ... " ${NUMPY_INCLUDE_DIR})

# Simple wrapper interface generator (SWIG) is required
FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})
MESSAGE("   SWIG executable     ... " ${SWIG_EXECUTABLE})
MESSAGE("   SWIG version        ... " ${SWIG_VERSION})
MESSAGE("   SWIG dir            ... " ${SWIG_DIR})

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

#PACKAGE_ADD_TEST_DIRECTORIES(test)

#PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

PACKAGE_POSTPROCESS()
