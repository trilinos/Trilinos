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
INCLUDE_DIRECTORIES(${DIR}/../Graph)
INCLUDE_DIRECTORIES(${DIR}/../Graph/BrickAggregation)
INCLUDE_DIRECTORIES(${DIR}/../Graph/UncoupledAggregation)
INCLUDE_DIRECTORIES(${DIR}/../Graph/UserAggregation)
INCLUDE_DIRECTORIES(${DIR}/../Headers)
INCLUDE_DIRECTORIES(${DIR}/../Misc)
INCLUDE_DIRECTORIES(${DIR}/../MueCentral)
INCLUDE_DIRECTORIES(${DIR}/../Rebalancing)
INCLUDE_DIRECTORIES(${DIR}/../Smoothers)
INCLUDE_DIRECTORIES(${DIR}/../Transfers/BaseClass)
IF (${PACKAGE_NAME}_ENABLE_Experimental)
  INCLUDE_DIRECTORIES(${DIR}/../Smoothers/BlockedSmoothers)
  INCLUDE_DIRECTORIES(${DIR}..//Transfers/BlockedTransfers)
ENDIF()
INCLUDE_DIRECTORIES(${DIR}/../Transfers/Generic)
INCLUDE_DIRECTORIES(${DIR}/../Transfers/Energy-Minimization)
INCLUDE_DIRECTORIES(${DIR}/../Transfers/Energy-Minimization/Solvers)
INCLUDE_DIRECTORIES(${DIR}/../Transfers/Geometric-Interpolation)
INCLUDE_DIRECTORIES(${DIR}/../Transfers/Smoothed-Aggregation)
INCLUDE_DIRECTORIES(${DIR}/../Transfers/Petrov-Galerkin-SA)
INCLUDE_DIRECTORIES(${DIR}/../Transfers/User)
INCLUDE_DIRECTORIES(${DIR}/../Transfers/SemiCoarsen)
INCLUDE_DIRECTORIES(${DIR}/../Utils)
INCLUDE_DIRECTORIES(${DIR}/../Utils/ForwardDeclaration)
if (${PACKAGE_NAME}_ENABLE_Tpetra)
  INCLUDE_DIRECTORIES(${DIR}/../../adapters/tpetra)
ENDIF()

###############################################
#The HEADERS variable is the list of headers that will be copied to the install directory during "make install".
SET(HEADERS "")
APPEND_GLOB(HEADERS ${BDIR}/${PACKAGE_NAME}_config.hpp)
APPEND_GLOB(HEADERS *.hpp)

###############################################
SET(SOURCES "")
APPEND_GLOB(SOURCES *.cpp)

#
# Explicit instantiation
#

# Header file
# TODO: make a loop
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/ NOSIERRABJAM)

INCLUDE_DIRECTORIES(${BDIR}/)
APPEND_GLOB(HEADERS ${BDIR}/*.hpp)

# Cpp file
IF(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
  APPEND_SET(SOURCES ${DIR}/../Utils/ExplicitInstantiation/MueLu_ParameterListInterpreter.cpp)
  APPEND_SET(SOURCES ${DIR}/../Utils/ExplicitInstantiation/MueLu_MLParameterListInterpreter.cpp)
  APPEND_SET(SOURCES ${DIR}/../Utils/ExplicitInstantiation/MueLu_FactoryFactory.cpp)
  APPEND_SET(SOURCES ${DIR}/../Utils/ExplicitInstantiation/MueLu_AdaptiveSaMLParameterListInterpreter.cpp)
ENDIF()

#
# Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  muelu-interface
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

# for debugging
#MESSAGE("")
#MESSAGE("**************** in muelu src ***********")
#MESSAGE("")
#MESSAGE("HEADERS = ${HEADERS}")
#MESSAGE("")
MESSAGE("SOURCES_INTERFACE = ${SOURCES}")

