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
if (${PACKAGE_NAME}_ENABLE_Epetra)
  INCLUDE_DIRECTORIES(${DIR}/../adapters/epetra)
ENDIF()
if (${PACKAGE_NAME}_ENABLE_Tpetra)
  INCLUDE_DIRECTORIES(${DIR}/../adapters/tpetra)
ENDIF()
if (${PACKAGE_NAME}_ENABLE_AmgX)
  INCLUDE_DIRECTORIES(${DIR}/../adapters/amgx)
ENDIF()

#########################################################
## TAW: 9/22/15
## CMake-based ETI system inspired by Ifpack2 ETI
## - auto generate *.cpp files for (all) MueLu classes
#########################################################

# Function to generate ETI (explicit template instantiation) files
# from a template and list of class names
## TODO introduce REGEX expressions here
## https://cmake.org/cmake/help/v3.0/command/string.html
FUNCTION(MUELU_PROCESS_ETI_TEMPLATE ETI_CLASSES TEMPLATE_FILE PROCESSED_FILE SOURCES_LIST)
  SET(SRCS "")
  FOREACH(CLASS ${ETI_CLASSES})
    # find spearator of class name and (optional) conditions
    string(FIND ${CLASS} "-" DASHPOS)
    IF ('${DASHPOS}' STREQUAL '-1')
      string(SUBSTRING ${CLASS} 0 ${DASHPOS} CLASS_NAME)
      #MESSAGE("New class name = ${CLASS_NAME}")
      set (CONDITION_NAME "")
      set (CONDITION_NAME_END "")      
    ELSE()
      string(SUBSTRING ${CLASS} 0 ${DASHPOS} CLASS_NAME)
      #MESSAGE("New class name = ${CLASS_NAME}")
      string(SUBSTRING ${CLASS} ${DASHPOS} -1 CONDITION_NAME)
      string(REPLACE "[" "(" CONDITION_NAME ${CONDITION_NAME})
      string(REPLACE "]" ")" CONDITION_NAME ${CONDITION_NAME})
      string(REPLACE "." " " CONDITION_NAME ${CONDITION_NAME})
      string(REPLACE "-" "" CONDITION_NAME ${CONDITION_NAME})
      string(REPLACE "?" "#" CONDITION_NAME ${CONDITION_NAME})
      string(STRIP CONDITION_NAME ${CONDITION_NAME})
      #MESSAGE("New condition name = ${CONDITION_NAME}")
      set (CONDITION_NAME_END "#endif")
    ENDIF()
    
    string(REPLACE "::" "_" CLASS_FILE_NAME "${CLASS_NAME}")
    string(TOUPPER "${CLASS_FILE_NAME}" UPPER_CASE_CLASS)
    string(REPLACE "CLASS_FILE_NAME" "${CLASS_FILE_NAME}" FINAL_FILE_NAME "${PROCESSED_FILE}")
    
    # the following generates one cpp file for all instantiations and 
    # enabled configurations
    CONFIGURE_FILE(${TEMPLATE_FILE} ${FINAL_FILE_NAME})
    #MESSAGE("Configure file ${FINAL_FILE_NAME} using template ${TEMPLATE_FILE}")
    SET(SRCS ${SRCS} ${FINAL_FILE_NAME})
    
  ENDFOREACH()
  SET(${SOURCES_LIST} ${SRCS} PARENT_SCOPE)
ENDFUNCTION(MUELU_PROCESS_ETI_TEMPLATE)

##
## MUELU_ETI_CPP_SOURCES collects all filenames of the MueLu ETI *.cpp files
## These are later added to the gloabl list of source files (SOURCES)
##
SET(MUELU_ETI_CPP_SOURCES "")

##
## only build ETI *.cpp files if ETI is enabled.
##
IF(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
  ##
  ## MUELU_SC_LO_GO_NO_ETI_CLASSES is a list of all class names which ETI 
  ## *.cpp files have to be generated for. Contains basically all MueLu 
  ## classes (with some exceptions, though).
  ##
  INCLUDE(Utils/ExplicitInstantiation/ETI_SC_LO_GO_NO_classes.cmake)

  ##
  ## MUELU_LO_GO_NO_ETI_CLASSES is a list of all class names which ETI 
  ## *.cpp files have to be generated for. Contains basically all MueLu 
  ## classes (with some exceptions, though).
  ##
  INCLUDE(Utils/ExplicitInstantiation/ETI_LO_GO_NO_classes.cmake)

  MUELU_PROCESS_ETI_TEMPLATE(
    "${MUELU_SC_LO_GO_NO_ETI_CLASSES}"
    Utils/ExplicitInstantiation/MueLu_ETI_SC_LO_GO_NO.tmpl
    "Utils/ExplicitInstantiation/CLASS_FILE_NAME.cpp"
    MueLu_ETI_SC_LO_GO_NO_SRC)
  LIST(APPEND MUELU_ETI_CPP_SOURCES ${MueLu_ETI_SC_LO_GO_NO_SRC})
  
  MUELU_PROCESS_ETI_TEMPLATE(
    "${MUELU_LO_GO_NO_ETI_CLASSES}"
    Utils/ExplicitInstantiation/MueLu_ETI_LO_GO_NO.tmpl
    "Utils/ExplicitInstantiation/CLASS_FILE_NAME.cpp"
    MueLu_ETI_LO_GO_NO_SRC)
  LIST(APPEND MUELU_ETI_CPP_SOURCES ${MueLu_ETI_LO_GO_NO_SRC})
ENDIF()

###############################################
#The HEADERS variable is the list of headers that will be copied to the install directory during "make install".
## TODO check this
SET(HEADERS "")
APPEND_GLOB(HEADERS ${BDIR}/${PACKAGE_NAME}_config.hpp)
APPEND_GLOB(HEADERS */*.hpp)
APPEND_GLOB(HEADERS */*/*.hpp)
APPEND_GLOB(HEADERS */*/*/*.hpp)

# Remove everything in the interface directory and add it to a list for muelu-interface instead
# NOTE: Because of the funny way the glob mactros work, the use of two slashses for the APPEND_GLOB is critical
SET(HEADERS_INTERFACE "")
APPEND_GLOB(HEADERS_INTERFACE Interface//*.hpp)
LIST(REMOVE_ITEM HEADERS ${HEADERS_INTERFACE})


###############################################
SET(SOURCES ${MUELU_ETI_CPP_SOURCES})
APPEND_GLOB(SOURCES */*.cpp)

#APPEND_GLOB(SOURCES ${MUELU_ETI_CPP_SOURCES})  # append ETI cpp sources (if any)

## TODO are the following lines necessary???
#APPEND_GLOB(SOURCES */*/*.cpp) -> ExplicitInstantiation files must be included only if explicit instantiation is enabled
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

#MESSAGE("ETI_SOURCES = ${MUELU_ETI_CPP_SOURCES}")
#MESSAGE("SOURCES = ${SOURCES}")

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
IF(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
  CONFIGURE_FILE(${DIR}/Utils/ExplicitInstantiation/${PACKAGE_NAME}_ExplicitInstantiation.hpp.in ${BDIR}/MueLu_ExplicitInstantiation.hpp)
  APPEND_GLOB(SOURCES ${DIR}/Utils/ExplicitInstantiation/*.cpp)
ENDIF()

# Remove everything in the interface directory and add it to a list for muelu-interface instead
# NOTE: Because of the funny way the glob mactros work, the use of two slashses for the APPEND_GLOB is critical
SET(SOURCES_INTERFACE "")
APPEND_GLOB(SOURCES_INTERFACE Interface//*.cpp)
IF(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
  APPEND_SET(SOURCES_INTERFACE ${DIR}/Utils/ExplicitInstantiation/MueLu_ParameterListInterpreter.cpp)
  APPEND_SET(SOURCES_INTERFACE ${DIR}/Utils/ExplicitInstantiation/MueLu_MLParameterListInterpreter.cpp)
  APPEND_SET(SOURCES_INTERFACE ${DIR}/Utils/ExplicitInstantiation/MueLu_FactoryFactory.cpp)
  APPEND_SET(SOURCES_INTERFACE ${DIR}/Utils/ExplicitInstantiation/MueLu_AdaptiveSaMLParameterListInterpreter.cpp)
ENDIF()
LIST(REMOVE_ITEM SOURCES ${SOURCES_INTERFACE})


if (NOT ${PACKAGE_NAME}_ENABLE_Ifpack2)
  LIST(REMOVE_ITEM HEADERS ${DIR}/Smoothers/MueLu_Ifpack2Smoother.hpp)
  LIST(REMOVE_ITEM SOURCES
    ${DIR}/Utils/ExplicitInstantiation/MueLu_Ifpack2Smoother.cpp)
ENDIF()

if (NOT ${PACKAGE_NAME}_ENABLE_Tpetra)
  LIST(REMOVE_ITEM SOURCES
    ${DIR}/Utils/ExplicitInstantiation/MueLu_TpetraOperator.cpp)
ENDIF()


#
# Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  muelu
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

# for debugging
#MESSAGE("")
#MESSAGE("**************** in muelu src ***********")
#MESSAGE("")
#MESSAGE("HEADERS = ${HEADERS}")
#MESSAGE("")
#MESSAGE("SOURCES = ${SOURCES}")

#MESSAGE("INCLUDE_DIRECTORIES = ")
#get_directory_property(output INCLUDE_DIRECTORIES)
#message(STATUS ${output})

#MESSAGE("")
#MESSAGE("**************** in muelu src ***********")
#MESSAGE("")
#MESSAGE("HEADERS_INTERFACE = ${HEADERS_INTERFACE}")
#MESSAGE("")
#MESSAGE("SOURCES_INTERFACE = ${SOURCES_INTERFACE}")
# touch CMakeLists.txt because a new file was created in Utils/ExplicitInstantiation of Utils/ForwardDeclaration
# touch CMakeLists.txt because a new file was created in Utils/ExplicitInstantiation of Utils/ForwardDeclaration
# touch CMakeLists.txt because a new file was created in Utils/ExplicitInstantiation of Utils/ForwardDeclaration
