# -*- cmake -*-

#
# Include the TriBITS library macros
INCLUDE(TribitsLibraryMacros)

#
# Package-specific configuration options
IF(Domi_ENABLE_Epetra)
  SET(HAVE_EPETRA ON)
ENDIF()

IF(Domi_ENABLE_Tpetra)
  SET(HAVE_TPETRA ON)
ENDIF()

#
# Here we set the Domi_ORDINAL_TYPE to be int.  This should be
# sufficient for all conceivable situations, even for global ordinals
# for distributed MDVectors.  We do not run into the same upper bound
# on problem size that, say, Epetra does, because Domi is a
# multi-dimensional package, and each dimension can be the maximum
# value supported by int.
SET(${PACKAGE_NAME}_ORDINAL_TYPE int
  CACHE STRING
  "Select the data type for Domi_Ordinal (Domi::Ordinal)."
  )

#
# Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# Set the include directories
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF (Domi_ENABLE_PyTrilinos)
  INCLUDE_DIRECTORIES(PYTHON_INCLUDE_PATH)
  INCLUDE_DIRECTORIES(NumPy_INCLUDE_DIR)
ENDIF (Domi_ENABLE_PyTrilinos)

ADD_SUBDIRECTORY(PyTrilinos)

#
# Define the header and source files (and directories)
SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

APPEND_SET(HEADERS
  Domi_ConfigDefs.hpp
  Domi_Version.hpp
  Domi_Utils.hpp
  Domi_Exceptions.hpp
  Domi_Slice.hpp
  Domi_MDIterator.hpp
  Domi_MDRevIterator.hpp
  Domi_MDArray.hpp
  Domi_MDArrayView.hpp
  Domi_MDArrayRCP.hpp
  Domi_MDComm.hpp
  Domi_MDMap.hpp
  Domi_MDVector.hpp
  Domi_getValidParameters.hpp
  )

APPEND_SET(SOURCES
  Domi_Version.cpp
  Domi_Utils.cpp
  Domi_Slice.cpp
  Domi_MDComm.cpp
  Domi_getValidParameters.cpp
  )

#
# Define the targets for package's library(s)
TRIBITS_ADD_LIBRARY(
  domi
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

#
# Python interface
IF (Domi_ENABLE_PyTrilinos)

  SET(CMAKE_SWIG_FLAGS "")
  GET_FILENAME_COMPONENT(PYTRILINOS_DOXDIR
    ${CMAKE_CURRENT_SOURCE_DIR}/../../../packages/PyTrilinos/doc/Doxygen ABSOLUTE
    )
  LIST(APPEND CMAKE_SWIG_FLAGS -I${PYTRILINOS_DOXDIR})
  LIST(APPEND CMAKE_SWIG_FLAGS -noproxydel)

  SET_SOURCE_FILES_PROPERTIES(
    Domi.i
    PROPERTIES CPLUSPLUS ON
    )

  SWIG_MODULE_GET_OUTDIR_AND_MODULE(Domi.i OUTDIR MODULE)

  #
  # On Mac OS X Gnu compilers, add dynamic lookup for undefined symbols
  # to the pytrilinos library and PyTrilinos extension modules
  SET(EXTRA_LINK_ARGS "${CMAKE_SHARED_LINKER_FLAGS}")
  IF(APPLE)
    IF((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR
       (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
      SET(EXTRA_LINK_ARGS "${EXTRA_LINK_ARGS} -undefined dynamic_lookup")
    ENDIF()
  ENDIF(APPLE)

  # Note that I use a module name PyDomi here because the Mac OS can't
  # handle domi.dir and Domi.dir being different directories
  SWIG_ADD_MODULE(PyDomi
    python
    Domi.i
    ${OUTDIR}
    ${MODULE}
    )

  SET_TARGET_PROPERTIES(PyDomi
    PROPERTIES
    LINK_FLAGS ${EXTRA_LINK_ARGS}
    NO_SONAME  1
    )

  SWIG_LINK_LIBRARIES(PyDomi pytrilinos domi)

  LIST(APPEND ADDITIONAL_CLEAN_FILES "${MODULE}PYTHON_wrap.h")
ENDIF (Domi_ENABLE_PyTrilinos)
