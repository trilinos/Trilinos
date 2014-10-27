INCLUDE(TrilinosCreateClientTemplateHeaders)

# Note about the use of wildcard in CMakeLists.txt: CMake dont't know
# when new files is added. You need to re-run CMake manually to
# regenerate a build system incorporating the new file.
# YOU MUST ALSO TOUCH A CMAKE CONFIGURATION FILE WHEN YOU PUSH THE NEW
# FILE TO FORCE THE RECONFIGURE ON OTHER PEOPLE'S BUILDS.

#
# Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# Define the header and source files (and directories)
#

SET(DIR  ${CMAKE_CURRENT_SOURCE_DIR})
SET(BDIR ${CMAKE_CURRENT_BINARY_DIR})


INCLUDE_DIRECTORIES(${DIR})
INCLUDE_DIRECTORIES(${DIR}/Graph)
INCLUDE_DIRECTORIES(${DIR}/Graph/BrickAggregation)
INCLUDE_DIRECTORIES(${DIR}/Graph/UncoupledAggregation)
INCLUDE_DIRECTORIES(${DIR}/Graph/UserAggregation)
INCLUDE_DIRECTORIES(${DIR}/Headers)
INCLUDE_DIRECTORIES(${DIR}/Interface)
INCLUDE_DIRECTORIES(${DIR}/Misc)
INCLUDE_DIRECTORIES(${DIR}/MueCentral)
INCLUDE_DIRECTORIES(${DIR}/Rebalancing)
INCLUDE_DIRECTORIES(${DIR}/Smoothers)
INCLUDE_DIRECTORIES(${DIR}/Transfers/BaseClass)
IF (${PACKAGE_NAME}_ENABLE_Experimental)
  INCLUDE_DIRECTORIES(${DIR}/Smoothers/BlockedSmoothers)
  INCLUDE_DIRECTORIES(${DIR}/Transfers/BlockedTransfers)
ENDIF()
INCLUDE_DIRECTORIES(${DIR}/Transfers/Generic)
INCLUDE_DIRECTORIES(${DIR}/Transfers/Energy-Minimization)
INCLUDE_DIRECTORIES(${DIR}/Transfers/Energy-Minimization/Solvers)
INCLUDE_DIRECTORIES(${DIR}/Transfers/Geometric-Interpolation)
INCLUDE_DIRECTORIES(${DIR}/Transfers/Smoothed-Aggregation)
INCLUDE_DIRECTORIES(${DIR}/Transfers/Petrov-Galerkin-SA)
INCLUDE_DIRECTORIES(${DIR}/Transfers/User)
INCLUDE_DIRECTORIES(${DIR}/Transfers/SemiCoarsen)
INCLUDE_DIRECTORIES(${DIR}/Utils)
INCLUDE_DIRECTORIES(${DIR}/Utils/ForwardDeclaration)
if (${PACKAGE_NAME}_ENABLE_Tpetra)
  INCLUDE_DIRECTORIES(${DIR}/../adapters/tpetra)
ENDIF()

#The HEADERS variable is the list of headers that will be copied to the install directory during "make install".
SET(HEADERS "")
APPEND_GLOB(HEADERS ${BDIR}/${PACKAGE_NAME}_config.hpp)
APPEND_GLOB(HEADERS */*.hpp)
APPEND_GLOB(HEADERS */*/*.hpp)
APPEND_GLOB(HEADERS */*/*/*.hpp)
# The Tpetra adapter headers are already pulled in by muelu-adapters.
#if (${PACKAGE_NAME}_ENABLE_Tpetra)
#  APPEND_GLOB(HEADERS ${DIR}/../adapters/tpetra/*.hpp)
#ENDIF()

SET(SOURCES "")
APPEND_GLOB(SOURCES */*.cpp)

#APPEND_GLOB(SOURCES */*/*.cpp) -> ExplicitInstantiation files must be included only if explicit instantiation is enabled
APPEND_GLOB(SOURCES Interface/*.cpp)
APPEND_GLOB(SOURCES Smoothers/BlockedSmoothers/*.cpp)
APPEND_GLOB(SOURCES Transfers/BaseClass/*.cpp)
APPEND_GLOB(SOURCES Transfers/BlockedTransfers/*.cpp)
APPEND_GLOB(SOURCES Transfers/Energy-Minimization/*.cpp)
APPEND_GLOB(SOURCES Transfers/Energy-Minimization/Solvers/*.cpp)
APPEND_GLOB(SOURCES Transfers/Generic/*.cpp)
APPEND_GLOB(SOURCES Transfers/Geometric-Interpolation/*.cpp)
APPEND_GLOB(SOURCES Transfers/Petrov-Galerkin-SA/*.cpp)
APPEND_GLOB(SOURCES Transfers/Smoothed-Aggregation/*.cpp)
APPEND_GLOB(SOURCES Utils/ForwardDeclaration/*.cpp)

#
# Explicit instantiation
#

# Header file
# TODO: make a loop
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Graph NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Graph/BrickAggregation NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Graph/UncoupledAggregation NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Graph/UserAggregation NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Interface NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Misc NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/MueCentral NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Rebalancing NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Smoothers NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Smoothers/BlockedSmoothers NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/BaseClass NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/BlockedTransfers NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/Generic NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/Energy-Minimization NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/Energy-Minimization/Solvers NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/Geometric-Interpolation NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/Petrov-Galerkin-SA NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/Smoothed-Aggregation NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/User NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/SemiCoarsen NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Utils NOSIERRABJAM)
if (${PACKAGE_NAME}_ENABLE_Tpetra)
  TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/../adapters/tpetra NOSIERRABJAM)
ENDIF()

INCLUDE_DIRECTORIES(${BDIR}/)
APPEND_GLOB(HEADERS ${BDIR}/*.hpp)

# Cpp file
SET(CU_SOURCES "")
IF(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
  CONFIGURE_FILE(${DIR}/Utils/ExplicitInstantiation/${PACKAGE_NAME}_ExplicitInstantiation.hpp.in ${BDIR}/MueLu_ExplicitInstantiation.hpp)
  APPEND_GLOB(SOURCES ${DIR}/Utils/ExplicitInstantiation/*.cpp)
  APPEND_GLOB(CU_SOURCES ${DIR}/Utils/ExplicitInstantiation/*.cu)
ENDIF()

if (NOT ${PACKAGE_NAME}_ENABLE_Ifpack2)
  LIST(REMOVE_ITEM HEADERS ${DIR}/Smoothers/MueLu_Ifpack2Smoother.hpp)
  LIST(REMOVE_ITEM SOURCES
    ${DIR}/Utils/ExplicitInstantiation/MueLu_Ifpack2Smoother.cpp
    ${DIR}/Utils/ExplicitInstantiation/MueLu_Ifpack2Smoother_Threads.cpp
    ${DIR}/Utils/ExplicitInstantiation/MueLu_Ifpack2Smoother_OpenMP.cpp)
  LIST(REMOVE_ITEM CU_SOURCES
    ${DIR}/Utils/ExplicitInstantiation/MueLu_Ifpack2Smoother_Cuda.cu)
ENDIF()

if (NOT ${PACKAGE_NAME}_ENABLE_Tpetra)
  LIST(REMOVE_ITEM SOURCES
    ${DIR}/Utils/ExplicitInstantiation/MueLu_TpetraOperator.cpp
    ${DIR}/Utils/ExplicitInstantiation/MueLu_TpetraOperator_Threads.cpp
    ${DIR}/Utils/ExplicitInstantiation/MueLu_TpetraOperator_OpenMP.cpp)
  LIST(REMOVE_ITEM CU_SOURCES
    ${DIR}/Utils/ExplicitInstantiation/MueLu_TpetraOperator_Cuda.cu)
ENDIF()


#
# Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  muelu
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

IF(KokkosClassic_ENABLE_Thrust OR
   (KokkosClassic_ENABLE_KokkosCompat AND Kokkos_ENABLE_Cuda))
  IF( TPL_ENABLE_MPI )
    CUDA_INCLUDE_DIRECTORIES( ${MPI_BASE_DIR}/include )
  ENDIF()
  TRIBITS_ADD_LIBRARY(
    muelu_cuda
    SOURCES ${CU_SOURCES}
    CUDALIBRARY
    )
ENDIF()

# for debugging
#MESSAGE("")
#MESSAGE("**************** in muelu src ***********")
#MESSAGE("")
#MESSAGE("HEADERS = ${HEADERS}")
#MESSAGE("")
#MESSAGE("INCLUDE_DIRECTORIES = ")
#get_directory_property(output INCLUDE_DIRECTORIES)
#message(STATUS ${output})
# touch CMakeLists.txt because a new file was created in Utils/ExplicitInstantiation of Utils/ForwardDeclaration
