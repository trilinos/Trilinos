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


TRIBITS_INCLUDE_DIRECTORIES(${DIR})
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Graph)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Graph/BrickAggregation)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Graph/UncoupledAggregation)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Graph/UserAggregation)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Headers)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Misc)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../MueCentral)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Rebalancing)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Smoothers)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Transfers/BaseClass)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Smoothers/BlockedSmoothers)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Transfers/BlockedTransfers)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Transfers/Generic)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Transfers/Energy-Minimization)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Transfers/Energy-Minimization/Solvers)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Transfers/GeneralGeometric)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Transfers/Geometric-Interpolation)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Transfers/Matrix-Free)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Transfers/Smoothed-Aggregation)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Transfers/Petrov-Galerkin-SA)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Transfers/User)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Transfers/SemiCoarsen)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Utils)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../Utils/ForwardDeclaration)
if (${PACKAGE_NAME}_ENABLE_Tpetra)
  TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../../adapters/tpetra)
ENDIF()

TRIBITS_INCLUDE_DIRECTORIES(${DIR}/FacadeClasses)


###############################################
#The HEADERS variable is the list of headers that will be copied to the install directory during "make install".
SET(HEADERS "")
APPEND_GLOB(HEADERS ${BDIR}/${PACKAGE_NAME}_config.hpp)
APPEND_GLOB(HEADERS *.hpp)
APPEND_GLOB(HEADERS ${DIR}/FacadeClasses/*.hpp)

###############################################
SET(SOURCES "")
APPEND_GLOB(SOURCES *.cpp)

#
# Explicit instantiation
#

# Header file
# TODO: make a loop
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/ NOSIERRABJAM)

TRIBITS_INCLUDE_DIRECTORIES(${BDIR}/)
APPEND_GLOB(HEADERS ${BDIR}/*.hpp)

IF (${PACKAGE_NAME}_ENABLE_SPLIT_ETI_CPP_FILES)

  # new ETI system which autogenerates splitted cpp files
  # not used by default since it might slow down the compilation process
  SET(MUELUINTERFACE_SC_LO_GO_NO_ETI_CLASSES )
  APPEND_SET(MUELUINTERFACE_SC_LO_GO_NO_ETI_CLASSES MueLu::FacadeClassFactory)
  APPEND_SET(MUELUINTERFACE_SC_LO_GO_NO_ETI_CLASSES MueLu::FacadeClassBase)
  APPEND_SET(MUELUINTERFACE_SC_LO_GO_NO_ETI_CLASSES MueLu::FacadeSimple)
  APPEND_SET(MUELUINTERFACE_SC_LO_GO_NO_ETI_CLASSES MueLu::BGS2x2)
  APPEND_SET(MUELUINTERFACE_SC_LO_GO_NO_ETI_CLASSES MueLu::ParameterListInterpreter)
  APPEND_SET(MUELUINTERFACE_SC_LO_GO_NO_ETI_CLASSES MueLu::MLParameterListInterpreter)
  APPEND_SET(MUELUINTERFACE_SC_LO_GO_NO_ETI_CLASSES MueLu::FactoryFactory)
  APPEND_SET(MUELUINTERFACE_SC_LO_GO_NO_ETI_CLASSES MueLu::AdaptiveSaMLParameterListInterpreter)

  MUELU_PROCESS_ETI_TEMPLATE(
    "../Utils/ExplicitInstantiation"
    "${MUELUINTERFACE_SC_LO_GO_NO_ETI_CLASSES}"
    "ExplicitInstantiation/CLASS_FILE_NAME"
    MueLuInterface_ETI_SC_LO_GO_NO_SRC
    0)
  LIST(APPEND SOURCES ${MueLuInterface_ETI_SC_LO_GO_NO_SRC})

ELSE()

# Old pseudo ETI system for muelu-interface: just use the cpp files that are contained in the respository

# Cpp file
IF(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
  APPEND_SET(SOURCES ${DIR}/../Utils/ExplicitInstantiation/ETI_MueLu_FacadeClassFactory.cpp)
  APPEND_SET(SOURCES ${DIR}/../Utils/ExplicitInstantiation/ETI_MueLu_FacadeClassBase.cpp)
  APPEND_SET(SOURCES ${DIR}/../Utils/ExplicitInstantiation/ETI_MueLu_FacadeSimple.cpp)
  APPEND_SET(SOURCES ${DIR}/../Utils/ExplicitInstantiation/ETI_MueLu_FacadeBGS2x2.cpp)
  APPEND_SET(SOURCES ${DIR}/../Utils/ExplicitInstantiation/ETI_MueLu_ParameterListInterpreter.cpp)
  APPEND_SET(SOURCES ${DIR}/../Utils/ExplicitInstantiation/ETI_MueLu_MLParameterListInterpreter.cpp)
  APPEND_SET(SOURCES ${DIR}/../Utils/ExplicitInstantiation/ETI_MueLu_FactoryFactory.cpp)
  APPEND_SET(SOURCES ${DIR}/../Utils/ExplicitInstantiation/ETI_MueLu_AdaptiveSaMLParameterListInterpreter.cpp)
ENDIF()

ENDIF() # split cpp files

#
# Define the targets for package's library(s)
#

# TRIBITS_ADD_LIBRARY(
#   muelu-interface
#   HEADERS ${HEADERS}
#   SOURCES ${SOURCES}
#   DEPLIBS muelu
#   )

# for debugging
#MESSAGE("")
#MESSAGE("**************** in muelu src ***********")
#MESSAGE("")
#MESSAGE("HEADERS = ${HEADERS}")
#MESSAGE("")
#MESSAGE("SOURCES_INTERFACE = ${SOURCES}")
