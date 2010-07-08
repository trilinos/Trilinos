
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

# Set the base Trilinos directory in order to find Trilinos_version.h

ASSERT_DEFINED(PROJECT_SOURCE_DIR)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

APPEND_SET(HEADERS
  ${PROJECT_SOURCE_DIR}/Trilinos_version.h
  )

# Set headers (for IDEs and installing)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

# NOTE: Above, we add the XXX_config.h file so tht MS Visual Studio
# will know to include this file in the dependency analysis

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Amesos2_Control.hpp
  Amesos2_Factory.hpp
  Amesos2_Factory_decl.hpp
  Amesos2_Factory_def.hpp
  Amesos2_FunctionMap.hpp
  Amesos2_MatrixAdapter.hpp
  Amesos2_MatrixAdapter_decl.hpp
  Amesos2_MatrixHelper.hpp
  Amesos2_MultiVecAdapter.hpp
  Amesos2_MultiVecAdapter_decl.hpp
  Amesos2_Solver.hpp
  Amesos2_SolverBase.hpp
  Amesos2_SolverBase_decl.hpp
  Amesos2_Solver_decl.hpp
  Amesos2_Solver_def.hpp
  Amesos2_Status.hpp
  Amesos2_Superlu.hpp
  Amesos2_Superlu_FunctionMap.hpp
  Amesos2_Superlu_MatrixHelper.hpp
  Amesos2_Superlu_TypeMap.hpp
  Amesos2_Superlu_decl.hpp
  Amesos2_Superlu_def.hpp
  Amesos2_Timers.hpp
  Amesos2_TpetraCrsMatrixAdapter.hpp
  Amesos2_TpetraCrsMatrixAdapter_decl.hpp
  Amesos2_TpetraCrsMatrixAdapter_def.hpp
  Amesos2_TpetraMultiVecAdapter.hpp
  Amesos2_TpetraMultiVecAdapter_decl.hpp
  Amesos2_TpetraMultiVecAdapter_def.hpp
  Amesos2_TypeMap.hpp
  Amesos2_Util.hpp
  Amesos2_Util_decl.hpp
  Amesos2_Util_def.hpp
  Amesos2_Version.hpp
  )

# Set sources

APPEND_SET(SOURCES
  Amesos2_Control.cpp
  Amesos2_Factory.cpp
  Amesos2_MatrixAdapter.cpp
  Amesos2_MultiVecAdapter.cpp
  Amesos2_Solver.cpp
  Amesos2_SolverBase.cpp
  Amesos2_Status.cpp
  Amesos2_Superlu.cpp
  Amesos2_TpetraCrsMatrixAdapter.cpp
  Amesos2_TpetraMultiVecAdapter.cpp
  Amesos2_Util.cpp
  Amesos2_Version.cpp
  )

# Set addition headers and source based on package-specific
# configuration options


#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  amesos2
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
