INCLUDE(TrilinosCreateClientTemplateHeaders)

# Function to generate ETI (explicit template instantiation) files
# from a template and list of class names
FUNCTION(IFPACK2_PROCESS_ETI_TEMPLATE ETI_CLASSES TEMPLATE_FILE PROCESSED_FILE SOURCES_LIST)
  SET(SRCS "")
  FOREACH(CLASS ${ETI_CLASSES})
    string(REPLACE "::" "_" CLASS_FILE_NAME "${CLASS}")
    string(TOUPPER "${CLASS_FILE_NAME}" UPPER_CASE_CLASS)
    string(REPLACE "CLASS_FILE_NAME" "${CLASS_FILE_NAME}" FINAL_FILE_NAME "${PROCESSED_FILE}")
    CONFIGURE_FILE(${TEMPLATE_FILE} ${FINAL_FILE_NAME})
    SET(SRCS ${SRCS} ${FINAL_FILE_NAME})
  ENDFOREACH()
  SET(${SOURCES_LIST} ${SRCS} PARENT_SCOPE)
ENDFUNCTION(IFPACK2_PROCESS_ETI_TEMPLATE)

# Function that follows the Tpetra convention for mangling C++ types
# so that they can be used as C preprocessor macro arguments.
FUNCTION(IFPACK2_MANGLE_TEMPLATE_PARAMETER TYPE_MANGLED_OUT TYPE_IN)
  STRING(REPLACE "<" "0" TMP0 "${TYPE_IN}")
  STRING(REPLACE ">" "0" TMP1 "${TMP0}")
  STRING(REPLACE "::" "_" TMP2 "${TMP1}")
  # Spaces (as in "long long") get squished out.
  STRING(REPLACE " " "" TMP3 "${TMP2}")
  SET(${TYPE_MANGLED_OUT} ${TMP3} PARENT_SCOPE)
ENDFUNCTION(IFPACK2_MANGLE_TEMPLATE_PARAMETER)

# Function that turns a valid Scalar template parameter into a macro
# name (all caps, with no white space and no punctuation other than
# underscore).
FUNCTION(IFPACK2_SCALAR_MACRO_NAME NAME_OUT NAME_IN)
  STRING(COMPARE EQUAL "${NAME_IN}" "__float128" IS_FLOAT128)
  IF(IS_FLOAT128)
    # __float128 is a special case; we remove the __ from the macro name.
    SET(${NAME_OUT} "FLOAT128" PARENT_SCOPE)
  ELSE()
    STRING(COMPARE EQUAL "${NAME_IN}" "std::complex<float>" IS_COMPLEX_FLOAT)
    IF(IS_COMPLEX_FLOAT)
      SET(${NAME_OUT} "COMPLEX_FLOAT" PARENT_SCOPE)
    ELSE()
      STRING(COMPARE EQUAL "${NAME_IN}" "std::complex<double>" IS_COMPLEX_DOUBLE)
      IF(IS_COMPLEX_DOUBLE)
        SET(${NAME_OUT} "COMPLEX_DOUBLE" PARENT_SCOPE)
      ELSE()
        # Convert to upper case, convert double colons to underscores,
        # and hope for the best.
        #
        # It would be nice if CMake were consistent about where output
        # arguments go.  Alas, this is not to be.  TOUPPER puts the
        # output argument last; REPLACE puts it after the search and
        # substitute strings, before the input string.
        STRING(TOUPPER "${NAME_IN}" TMP0)
        STRING(REPLACE "::" "_" TMP1 "${TMP0}")
        SET(${NAME_OUT} ${TMP1} PARENT_SCOPE)
      ENDIF()
    ENDIF()
  ENDIF()
ENDFUNCTION(IFPACK2_SCALAR_MACRO_NAME)

# Function that turns a valid LocalOrdinal or GlobalOrdinal template
# parameter (which must be a C++ built-in integer type) into a macro
# name (all caps, with no white space and no punctuation other than
# underscore).
FUNCTION(IFPACK2_ORDINAL_MACRO_NAME NAME_OUT NAME_IN)
  # It would be nice if CMake were consistent about where output
  # arguments go.  Alas, this is not to be.  TOUPPER puts the
  # output argument last; REPLACE puts it after the search and
  # substitute strings, before the input string.
  #
  # "long long" -> "LONG LONG"
  STRING(TOUPPER "${NAME_IN}" TMP0)
  # "LONG LONG" -> "LONG_LONG"
  STRING(REPLACE " " "_" TMP1 "${TMP0}")
  SET(${NAME_OUT} ${TMP1} PARENT_SCOPE)
ENDFUNCTION(IFPACK2_ORDINAL_MACRO_NAME)

# Function that turns a valid Node template parameter into a macro
# name (all caps, with no white space and no punctuation other than
# underscore).
FUNCTION(IFPACK2_NODE_MACRO_NAME NAME_OUT NAME_IN)
  STRING(REGEX MATCH "Kokkos::Compat::Kokkos(.*)WrapperNode" TMP0 "${NAME_IN}")
  STRING(COMPARE EQUAL "${TMP0}" "" DOES_NOT_MATCH)
  IF(DOES_NOT_MATCH)
    MESSAGE(FATAL_ERROR "Ifpack2: Node $NAME_IN is not a supported Node type.")
  ELSE()
    # Extract the Kokkos execution space (KOKKOS_EXEC_SPACE) from the Node name.
    STRING(REGEX REPLACE "Kokkos::Compat::Kokkos(.*)WrapperNode" "\\1" KOKKOS_EXEC_SPACE "${NAME_IN}")

    # Special case: Threads.  The macro name unfortunately differs
    # from the execution space name in a way that doesn't fit the
    # pattern of the other execution spaces.
    STRING(COMPARE EQUAL "${KOKKOS_EXEC_SPACE}" "Threads" IS_THREADS)
    IF(IS_THREADS)
      SET(${NAME_OUT} "PTHREAD" PARENT_SCOPE)
    ELSE()
      # The other cases (Cuda, Serial, OpenMP) are easy.
      STRING(TOUPPER "${KOKKOS_EXEC_SPACE}" NAME_OUT_TMP)
      SET(${NAME_OUT} ${NAME_OUT_TMP} PARENT_SCOPE)
    ENDIF()
  ENDIF()
ENDFUNCTION(IFPACK2_NODE_MACRO_NAME)

# Function to generate one .cpp file for the given (Scalar,
# LocalOrdinal, GlobalOrdinal, Node) template parameter combination,
# for run-time registration of Ifpack2's LinearSolverFactory over
# those template parameters.  This is meant to be called by
# IFPACK2_PROCESS_ALL_FACTORY_TEMPLATES.  This function takes the
# names already mangled, to avoid unnecessary string processing
# overhead.
FUNCTION(IFPACK2_PROCESS_ONE_FACTORY_TEMPLATE OUTPUT_FILE TEMPLATE_FILE SC LO GO NT SC_MANGLED LO_MANGLED GO_MANGLED NT_MANGLED SC_MACRO_NAME LO_MACRO_NAME GO_MACRO_NAME NT_MACRO_NAME)
  STRING(REPLACE ".tmpl" "_${SC_MACRO_NAME}_${LO_MACRO_NAME}_${GO_MACRO_NAME}_${NT_MACRO_NAME}.cpp" OUT_FILE "${TEMPLATE_FILE}")
  CONFIGURE_FILE("${TEMPLATE_FILE}" "${OUT_FILE}")

  SET(${OUTPUT_FILE} ${OUT_FILE} PARENT_SCOPE)
ENDFUNCTION(IFPACK2_PROCESS_ONE_FACTORY_TEMPLATE)

# Function to generate .cpp files for run-time registration of
# Ifpack2's LinearSolverFactory, from a template.  We generate one
# .cpp file for each (Scalar, LocalOrdinal, GlobalOrdinal, Node) type
# combination over which Ifpack2 does tests and/or ETI.
FUNCTION(IFPACK2_PROCESS_ALL_FACTORY_TEMPLATES OUTPUT_FILES TEMPLATE_FILE SCALAR_TYPES LOCALORDINAL_TYPES GLOBALORDINAL_TYPES NODE_TYPES)
  SET(OUT_FILES "")
  FOREACH(NT ${NODE_TYPES})
    IFPACK2_MANGLE_TEMPLATE_PARAMETER(NT_MANGLED "${NT}")
    IFPACK2_NODE_MACRO_NAME(NT_MACRO_NAME "${NT}")
    FOREACH(GO ${GLOBALORDINAL_TYPES})
      IFPACK2_MANGLE_TEMPLATE_PARAMETER(GO_MANGLED "${GO}")
      IFPACK2_ORDINAL_MACRO_NAME(GO_MACRO_NAME "${GO}")
      FOREACH(LO ${LOCALORDINAL_TYPES})
        IFPACK2_MANGLE_TEMPLATE_PARAMETER(LO_MANGLED "${LO}")
        IFPACK2_ORDINAL_MACRO_NAME(LO_MACRO_NAME "${LO}")
        FOREACH(SC ${SCALAR_TYPES})
          IFPACK2_MANGLE_TEMPLATE_PARAMETER(SC_MANGLED "${SC}")
          IFPACK2_SCALAR_MACRO_NAME(SC_MACRO_NAME "${SC}")
          IFPACK2_PROCESS_ONE_FACTORY_TEMPLATE(OUT_FILE "${TEMPLATE_FILE}" "${SC}" "${LO}" "${GO}" "${NT}" "${SC_MANGLED}" "${LO_MANGLED}" "${GO_MANGLED}" "${NT_MANGLED}" "${SC_MACRO_NAME}" "${LO_MACRO_NAME}" "${GO_MACRO_NAME}" "${NT_MACRO_NAME}")
          LIST(APPEND OUT_FILES ${OUT_FILE})
        ENDFOREACH() # SC
      ENDFOREACH() # LO
    ENDFOREACH() # GO
  ENDFOREACH() # NT

  # This is the standard CMake idiom for setting an output variable so
  # that the caller can see the result.
  SET(${OUTPUT_FILES} ${OUT_FILES} PARENT_SCOPE)
ENDFUNCTION(IFPACK2_PROCESS_ALL_FACTORY_TEMPLATES)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# SupportGraph requires some TPLs and options to be enabled.
IF (${PACKAGE_NAME}_ENABLE_Experimental AND ${PACKAGE_NAME}_ENABLE_Lemon AND ${PACKAGE_NAME}_ENABLE_Amesos2 AND ${PACKAGE_NAME}_ENABLE_Cholmod)
  TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/supportgraph)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/supportgraph)
  SET(HEADERS ${HEADERS}
    supportgraph/Ifpack2_SupportGraph_decl.hpp
    supportgraph/Ifpack2_SupportGraph_def.hpp
    )
  SET(SOURCES ${SOURCES}
    supportgraph/Ifpack2_SupportGraph.cpp
    )
ENDIF()

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_config.h )
APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_ETIHelperMacros.h )

# Automatically generate ETI (explicit template instanatiation) files
# for Node types that use the (new) Kokkos Devices.
SET(IFPACK2_CPP_SOURCES "")
IF(Ifpack2_ENABLE_EXPLICIT_INSTANTIATION)
  # Set the list of Ifpack2 classes templated on <Scalar, LO, GO,
  # Node> for which we want to do ETI using this system.  These
  # classes usually operate on sparse matrices (instances of
  # Tpetra::CrsMatrix or Tpetra::RowMatrix, which also take these
  # template parameters).
  GLOBAL_SET(IFPACK2_ETI_CLASSES
    AdditiveSchwarz
    BlockRelaxation
    BorderedOperator
    Chebyshev
    BandedContainer
    DenseContainer
    Diagonal
    DiagonalFilter
    DropFilter
    Factory
    Hiptmair
    IdentitySolver
    ILUT
    LinePartitioner
    LocalFilter
    LocalSparseTriangularSolver
    OverlappingRowMatrix
    Relaxation
    ReorderFilter
    RILUK
    SingletonFilter
    SparseContainer
    SparsityFilter
    ContainerFactory
    TriDiContainer
    Details::Chebyshev
    Details::DenseSolver
    Details::Factory
    Details::LinearSolver
    Details::OneLevelFactory
    Details::TriDiSolver
    Experimental::RBILUK
    )
  IF(Ifpack2_ENABLE_Amesos2)
    APPEND_GLOBAL_SET(IFPACK2_ETI_CLASSES
      Details::Amesos2Wrapper
      )
  ENDIF()
  IF(Ifpack2_ENABLE_DEPRECATED_CODE)
    APPEND_GLOBAL_SET(IFPACK2_ETI_CLASSES Krylov)
  ENDIF()
  IF(${PACKAGE_NAME}_ENABLE_Experimental_KokkosKernels_Features)
    APPEND_GLOBAL_SET(IFPACK2_ETI_CLASSES BlockTriDiContainer)
  ENDIF()

  IF(Ifpack2_ENABLE_ShyLU_NodeFastILU)
    APPEND_GLOBAL_SET(IFPACK2_ETI_CLASSES
      Details::FastILU_Base
      Details::Filu
      Details::Fic
      Details::Fildl
    )
  ENDIF()

  # Set the list of Ifpack2 classes templated on <LO, GO, Node> for
  # which we want to do ETI using this system.  These classes usually
  # operate on sparse graphs (instances of Tpetra::CrsGraph or
  # Tpetra::RowGraph, which also take these template parameters).
  GLOBAL_SET(IFPACK2_ETI_LO_GO_CLASSES
    LinearPartitioner
    OverlappingPartitioner
    Details::OverlappingRowGraph
    )

  # ETI for the Kokkos::Serial device.
  IF (Tpetra_INST_SERIAL)
    IFPACK2_PROCESS_ETI_TEMPLATE(
      "${IFPACK2_ETI_CLASSES}"
      Ifpack2_ETI_SC_LO_GO_Serial.tmpl
      "Ifpack2_CLASS_FILE_NAME_Serial.cpp"
      IFPACK2_SRCS)
    LIST(APPEND IFPACK2_CPP_SOURCES ${IFPACK2_SRCS})
    IFPACK2_PROCESS_ETI_TEMPLATE(
      "${IFPACK2_ETI_LO_GO_CLASSES}"
      Ifpack2_ETI_LO_GO_Serial.tmpl
      "Ifpack2_CLASS_FILE_NAME_Serial.cpp"
      IFPACK2_SRCS)
    LIST(APPEND IFPACK2_CPP_SOURCES ${IFPACK2_SRCS})
  ENDIF()

  # ETI for the Kokkos::Classic::Serialnode device.
  IF (Tpetra_INST_SERIALCLASSIC)
    IFPACK2_PROCESS_ETI_TEMPLATE(
      "${IFPACK2_ETI_CLASSES}"
      Ifpack2_ETI_SC_LO_GO_Classic.tmpl
      "Ifpack2_CLASS_FILE_NAME_Classic.cpp"
      IFPACK2_SRCS)
    LIST(APPEND IFPACK2_CPP_SOURCES ${IFPACK2_SRCS})
    IFPACK2_PROCESS_ETI_TEMPLATE(
      "${IFPACK2_ETI_LO_GO_CLASSES}"
      Ifpack2_ETI_LO_GO_Classic.tmpl
      "Ifpack2_CLASS_FILE_NAME_Classic.cpp"
      IFPACK2_SRCS)
    LIST(APPEND IFPACK2_CPP_SOURCES ${IFPACK2_SRCS})
  ENDIF()

  # ETI for the Kokkos::Threads device.
  IF (Tpetra_INST_PTHREAD)
    IFPACK2_PROCESS_ETI_TEMPLATE(
      "${IFPACK2_ETI_CLASSES}"
      Ifpack2_ETI_SC_LO_GO_Threads.tmpl
      "Ifpack2_CLASS_FILE_NAME_Threads.cpp"
      IFPACK2_SRCS)
    LIST(APPEND IFPACK2_CPP_SOURCES ${IFPACK2_SRCS})
    IFPACK2_PROCESS_ETI_TEMPLATE(
      "${IFPACK2_ETI_LO_GO_CLASSES}"
      Ifpack2_ETI_LO_GO_Threads.tmpl
      "Ifpack2_CLASS_FILE_NAME_Threads.cpp"
      IFPACK2_SRCS)
    LIST(APPEND IFPACK2_CPP_SOURCES ${IFPACK2_SRCS})
  ENDIF()

  # ETI for the Kokkos::OpenMP device.
  IF (Tpetra_INST_OPENMP)
    IFPACK2_PROCESS_ETI_TEMPLATE(
      "${IFPACK2_ETI_CLASSES}"
      Ifpack2_ETI_SC_LO_GO_OpenMP.tmpl
      "Ifpack2_CLASS_FILE_NAME_OpenMP.cpp"
      IFPACK2_SRCS)
    LIST(APPEND IFPACK2_CPP_SOURCES ${IFPACK2_SRCS})
    IFPACK2_PROCESS_ETI_TEMPLATE(
      "${IFPACK2_ETI_LO_GO_CLASSES}"
      Ifpack2_ETI_LO_GO_OpenMP.tmpl
      "Ifpack2_CLASS_FILE_NAME_OpenMP.cpp"
      IFPACK2_SRCS)
    LIST(APPEND IFPACK2_CPP_SOURCES ${IFPACK2_SRCS})
  ENDIF()

  # ETI for the Kokkos::Cuda device.
  IF (Tpetra_INST_CUDA)
    IFPACK2_PROCESS_ETI_TEMPLATE(
      "${IFPACK2_ETI_CLASSES}"
      Ifpack2_ETI_SC_LO_GO_Cuda.tmpl
      "Ifpack2_CLASS_FILE_NAME_Cuda.cpp"
      IFPACK2_SRCS)
    LIST(APPEND IFPACK2_CPP_SOURCES ${IFPACK2_SRCS})
    IFPACK2_PROCESS_ETI_TEMPLATE(
      "${IFPACK2_ETI_LO_GO_CLASSES}"
      Ifpack2_ETI_LO_GO_Cuda.tmpl
      "Ifpack2_CLASS_FILE_NAME_Cuda.cpp"
      IFPACK2_SRCS)
    LIST(APPEND IFPACK2_CPP_SOURCES ${IFPACK2_SRCS})
  ENDIF()
ENDIF()

# Generate .cpp files for LinearSolverFactory registration,
# from the template (.tmpl) file in this directory.
#
# mfh 03 Aug 2015: Bug 6380 is fixed now.
SET(IFPACK2_FACTORY_NODE_TYPES ${Ifpack2_ETI_NODES})
SET(IFPACK2_FACTORY_LOCALORDINAL_TYPES ${Ifpack2_ETI_LORDS})
SET(IFPACK2_FACTORY_GLOBALORDINAL_TYPES ${Ifpack2_ETI_GORDS})
SET(IFPACK2_FACTORY_SCALAR_TYPES ${Ifpack2_ETI_SCALARS})

SET(TEMPLATE_FILE "Ifpack2_Details_LinearSolverFactory.tmpl")
IFPACK2_PROCESS_ALL_FACTORY_TEMPLATES(IFPACK2_FACTORY_OUTPUT_FILES "${TEMPLATE_FILE}" "${IFPACK2_FACTORY_SCALAR_TYPES}" "${IFPACK2_FACTORY_LOCALORDINAL_TYPES}" "${IFPACK2_FACTORY_GLOBALORDINAL_TYPES}" "${IFPACK2_FACTORY_NODE_TYPES}")
APPEND_SET(IFPACK2_CPP_SOURCES ${IFPACK2_FACTORY_OUTPUT_FILES})

#MESSAGE(DEBUG " *** IFPACK2_CPP_SOURCES = ${IFPACK2_CPP_SOURCES}")

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  ifpack2
  HEADERS ${HEADERS}
  SOURCES ${SOURCES} ${IFPACK2_CPP_SOURCES}
  )

#
# Make a trivial change here if you want CMake to run, due to changes
# you make to files in Ifpack2.  Here is another such change.
#
