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
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Graph)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Graph/BrickAggregation)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Graph/Containers)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Graph/HybridAggregation)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Graph/MatrixTransformation)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Graph/PairwiseAggregation)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Graph/StructuredAggregation)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Graph/StructuredAggregation/coupled)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Graph/StructuredAggregation/uncoupled)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Graph/UncoupledAggregation)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Graph/UserAggregation)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Headers)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Interface)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Interface/FacadeClasses)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Misc)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/MueCentral)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Operators)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Rebalancing)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Smoothers)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Transfers/BaseClass)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Transfers/BlackBox)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Smoothers/BlockedSmoothers)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Transfers/BlockedTransfers)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Transfers/Classical)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Transfers/Energy-Minimization)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Transfers/Energy-Minimization/Solvers)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Transfers/GeneralGeometric)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Transfers/Generic)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Transfers/Geometric-Interpolation)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Transfers/Matrix-Free)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Transfers/Petrov-Galerkin-SA)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Transfers/SemiCoarsen)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Transfers/Smoothed-Aggregation)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Transfers/User)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Utils)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Utils/ForwardDeclaration)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../adapters/xpetra)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../adapters/tpetra)
IF (${PACKAGE_NAME}_ENABLE_Epetra)
  TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../adapters/epetra)
ENDIF()
IF (${PACKAGE_NAME}_ENABLE_AmgX)
  TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../adapters/amgx)
ENDIF()
IF (${PACKAGE_NAME}_ENABLE_Intrepid2)
  TRIBITS_INCLUDE_DIRECTORIES(${DIR}/Transfers/PCoarsen)
ENDIF()
IF (${PACKAGE_NAME}_ENABLE_Experimental)
  TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../research/regionMG/src)
ENDIF()
IF (${PACKAGE_NAME}_ENABLE_Stratimikos)
   TRIBITS_INCLUDE_DIRECTORIES(${DIR}/../adapters/stratimikos)
ENDIF()

# Function to generate ETI (explicit template instantiation) files
# from a template and list of class names

IF (${PACKAGE_NAME}_ENABLE_SPLIT_ETI_CPP_FILES)

#########################################################
## TAW: 2/1/16
## CMake-based ETI system inspired by Tpetra ETI
## - auto generate *.cpp files for (all) MueLu classes
## - split *.cpp files into extra cpp files for
##   each configuration. This may slow down the compilation
##   process.
#########################################################
FUNCTION(MUELU_PROCESS_ETI_TEMPLATE TEMPLATE_FILE_DIR ETI_CLASSES PROCESSED_FILE SOURCES_LIST SKIP_SCALAR_INSTANTIATION)

  # Make sure that Tpetra actually defined these variables, even if they
  # are empty.
  # TODO: workaround for Tpetra is disabled but Epetra is enabled?
  ASSERT_DEFINED(TpetraCore_ETI_SCALARS)
  ASSERT_DEFINED(TpetraCore_ETI_LORDS)
  ASSERT_DEFINED(TpetraCore_ETI_GORDS)
  ASSERT_DEFINED(TpetraCore_ETI_NODES)

  # Exclude GlobalOrdinal types from the list of Scalar types over which
  # MueLu does ETI.
  #
  # Look at each entry in Tpetra's list of ETI Scalar types.  If it is
  # also a Tpetra ETI GlobalOrdinal type, don't append it to MueLu's
  # list of ETI Scalar types.  Otherwise, do.
  SET (MueLu_ETI_SCALARS_TEMP ) # start with the empty list
  FOREACH (scalarName IN LISTS TpetraCore_ETI_SCALARS)
    LIST (FIND TpetraCore_ETI_GORDS ${scalarName} scalarNameIndex)
    IF (-1 EQUAL scalarNameIndex) # not in the list of integer types
      LIST (APPEND MueLu_ETI_SCALARS_TEMP ${scalarName})
    ENDIF ()
  ENDFOREACH ()

  SET(MueLu_ETI_SCALARS ${MueLu_ETI_SCALARS_TEMP})
  SET(MueLu_ETI_LORDS   ${TpetraCore_ETI_LORDS})
  SET(MueLu_ETI_GORDS   ${TpetraCore_ETI_GORDS})
  SET(MueLu_ETI_NODES   ${TpetraCore_ETI_NODES})

  #MESSAGE("MueLu_ETI_SCALARS = ${MueLu_ETI_SCALARS}")
  #MESSAGE("MueLu_ETI_LORDS   = ${MueLu_ETI_LORDS}")
  #MESSAGE("MueLu_ETI_GORDS   = ${MueLu_ETI_GORDS}")
  #MESSAGE("MueLu_ETI_NODES   = ${MueLu_ETI_NODES}")

  # check whether both Epetra and Tpetra are enabled
  SET(MueLu_ETI_EpetraGuards "")
  SET(MueLu_ETI_EpetraGuards_END "")
  IF (${PACKAGE_NAME}_ENABLE_Epetra)
  # both Epetra and Tpetra are enabled. The Epetra instantiation needs to be guarded
  SET(MueLu_ETI_EpetraGuards "#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))")
    SET(MueLu_ETI_EpetraGuards_END "#endif")
  ENDIF()

  # loop over all ETI_CLASSES
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

    # If Epetra is enabled, instantiate it on <double,int,int,EpetraNode>
    # Note: no Epetra64 support in MueLu
    # EpetraNode is defined in the template file depending on the enabled
    # Kokkos nodes (either Serial or OpenMP)
    IF (${PACKAGE_NAME}_ENABLE_Epetra)
      IF(${SKIP_SCALAR_INSTANTIATION})
        CONFIGURE_FILE("${TEMPLATE_FILE_DIR}/MueLu_ETI_LO_GO_NO_Epetra.tmpl" ${FINAL_FILE_NAME}_Epetra.cpp)
        SET(SRCS ${SRCS} ${FINAL_FILE_NAME}_Epetra.cpp)
      ELSE()
        CONFIGURE_FILE("${TEMPLATE_FILE_DIR}/MueLu_ETI_SC_LO_GO_NO_Epetra.tmpl" ${FINAL_FILE_NAME}_Epetra.cpp)
        SET(SRCS ${SRCS} ${FINAL_FILE_NAME}_Epetra.cpp)
      ENDIF()
    ENDIF()

    # Create the cpp file for all combinations of template parameters
    FOREACH(NT ${MueLu_ETI_NODES})
      TPETRA_MANGLE_TEMPLATE_PARAMETER(NT_MANGLED "${NT}")
      TPETRA_NODE_MACRO_NAME(NT_MACRO_NAME "${NT}")
      FOREACH(GO ${MueLu_ETI_GORDS})
        TPETRA_MANGLE_TEMPLATE_PARAMETER(GO_MANGLED "${GO}")
        TPETRA_SLG_MACRO_NAME(GO_MACRO_NAME "${GO}")
        FOREACH(LO ${MueLu_ETI_LORDS})
          TPETRA_MANGLE_TEMPLATE_PARAMETER(LO_MANGLED "${LO}")
          TPETRA_SLG_MACRO_NAME(LO_MACRO_NAME "${LO}")

          IF(${SKIP_SCALAR_INSTANTIATION})
            CONFIGURE_FILE("${TEMPLATE_FILE_DIR}/MueLu_ETI_LO_GO_NO_Tpetra.tmpl" ${FINAL_FILE_NAME}_${LO_MACRO_NAME}_${GO_MACRO_NAME}_${NT_MACRO_NAME}.cpp)
            SET(SRCS ${SRCS} ${FINAL_FILE_NAME}_${LO_MACRO_NAME}_${GO_MACRO_NAME}_${NT_MACRO_NAME}.cpp)
          ELSE()
            FOREACH(SC ${MueLu_ETI_SCALARS})
            TPETRA_MANGLE_TEMPLATE_PARAMETER(SC_MANGLED "${SC}")
            TPETRA_SLG_MACRO_NAME(SC_MACRO_NAME "${SC}")
            TPETRA_SC_MACRO_EXPR(SC_MACRO_EXPR "${SC}" "${GO}" "${SC_MACRO_NAME}")

            CONFIGURE_FILE("${TEMPLATE_FILE_DIR}/MueLu_ETI_SC_LO_GO_NO_Tpetra.tmpl" ${FINAL_FILE_NAME}_${SC_MACRO_NAME}_${LO_MACRO_NAME}_${GO_MACRO_NAME}_${NT_MACRO_NAME}.cpp)
          SET(SRCS ${SRCS} ${FINAL_FILE_NAME}_${SC_MACRO_NAME}_${LO_MACRO_NAME}_${GO_MACRO_NAME}_${NT_MACRO_NAME}.cpp)
            ENDFOREACH() # SC
          ENDIF() # skip loop over SC
        ENDFOREACH() # LO
      ENDFOREACH() # GO
    ENDFOREACH() # NT
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
  GLOBAL_SET(MUELU_SC_LO_GO_NO_ETI_CLASSES  )
  INCLUDE(Utils/ExplicitInstantiation/ETI_SC_LO_GO_NO_classes.cmake)

  ##
  ## MUELU_LO_GO_NO_ETI_CLASSES is a list of all class names which ETI
  ## *.cpp files have to be generated for. Contains basically all MueLu
  ## classes (with some exceptions, though).
  ##
  GLOBAL_SET ( MUELU_LO_GO_NO_ETI_CLASSES  )
  INCLUDE(Utils/ExplicitInstantiation/ETI_LO_GO_NO_classes.cmake)

  MUELU_PROCESS_ETI_TEMPLATE(
    "Utils/ExplicitInstantiation"
    "${MUELU_SC_LO_GO_NO_ETI_CLASSES}"
    "Utils/ExplicitInstantiation/CLASS_FILE_NAME"
    MueLu_ETI_SC_LO_GO_NO_SRC
    0)
  LIST(APPEND MUELU_ETI_CPP_SOURCES ${MueLu_ETI_SC_LO_GO_NO_SRC})

  MUELU_PROCESS_ETI_TEMPLATE(
    "Utils/ExplicitInstantiation"
    "${MUELU_LO_GO_NO_ETI_CLASSES}"
    "Utils/ExplicitInstantiation/CLASS_FILE_NAME"
    MueLu_ETI_LO_GO_NO_SRC
    1)
  LIST(APPEND MUELU_ETI_CPP_SOURCES ${MueLu_ETI_LO_GO_NO_SRC})
ENDIF()

ELSE() # do not split cpp files

#########################################################
## TAW: 9/22/15
## CMake-based ETI system inspired by Ifpack2 ETI
## - auto generate *.cpp files for (all) MueLu classes
## - do not split *.cpp files into extra cpp files for
##   each configuration
## - This is the default ETI generation system
#########################################################

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
    string(REPLACE "CLASS_FILE_NAME" "ETI_${CLASS_FILE_NAME}" FINAL_FILE_NAME "${PROCESSED_FILE}")

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
  GLOBAL_SET(MUELU_SC_LO_GO_NO_ETI_CLASSES  )
  INCLUDE(Utils/ExplicitInstantiation/ETI_SC_LO_GO_NO_classes.cmake)

  ##
  ## MUELU_LO_GO_NO_ETI_CLASSES is a list of all class names which ETI
  ## *.cpp files have to be generated for. Contains basically all MueLu
  ## classes (with some exceptions, though).
  ##
  GLOBAL_SET ( MUELU_LO_GO_NO_ETI_CLASSES  )
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

ENDIF() # split ETI files

###############################################
#The HEADERS variable is the list of headers that will be copied to the install directory during "make install".
SET(HEADERS "")
APPEND_GLOB(HEADERS ${BDIR}/${PACKAGE_NAME}_config.hpp)
APPEND_GLOB(HEADERS */*.hpp)
APPEND_GLOB(HEADERS */*/*.hpp)
APPEND_GLOB(HEADERS */*/*/*.hpp)

###############################################
SET(SOURCES ${MUELU_ETI_CPP_SOURCES})
APPEND_GLOB(SOURCES */*.cpp)
APPEND_GLOB(SOURCES ${DIR}/Graph/Containers/MueLu_LinkedList.cpp)

#
# Explicit instantiation
#

# Header file
# TODO: make a loop
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Graph NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Graph/BrickAggregation NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Graph/Containers NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Graph/HybridAggregation NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Graph/MatrixTransformation NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Graph/PairwiseAggregation NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Graph/StructuredAggregation NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Graph/StructuredAggregation/coupled NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Graph/StructuredAggregation/uncoupled NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Graph/UncoupledAggregation NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Graph/UserAggregation NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Interface NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Interface/FacadeClasses NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Misc NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/MueCentral NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Operators NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Rebalancing NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Smoothers NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Smoothers/BlockedSmoothers NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/BaseClass NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/BlackBox NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/BlockedTransfers NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/Classical NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/Generic NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/Energy-Minimization NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/Energy-Minimization/Solvers NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/GeneralGeometric NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/Geometric-Interpolation NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/Matrix-Free NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/Petrov-Galerkin-SA NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/Smoothed-Aggregation NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/User NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/SemiCoarsen NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Utils NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/../adapters/xpetra NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/../adapters/tpetra NOSIERRABJAM)
IF (${PACKAGE_NAME}_ENABLE_Intrepid2)
  TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Transfers/PCoarsen NOSIERRABJAM)
ENDIF()

# Cpp file
IF(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
  CONFIGURE_FILE(${DIR}/Utils/ExplicitInstantiation/${PACKAGE_NAME}_ExplicitInstantiation.hpp.in ${BDIR}/MueLu_ExplicitInstantiation.hpp)
  APPEND_GLOB(SOURCES ${DIR}/Utils/ExplicitInstantiation/*.cpp)
ENDIF()

TRIBITS_INCLUDE_DIRECTORIES(${BDIR}/)
APPEND_GLOB(HEADERS ${BDIR}/*.hpp)


if (NOT ${PACKAGE_NAME}_ENABLE_Ifpack2)
  LIST(REMOVE_ITEM HEADERS ${DIR}/Smoothers/MueLu_Ifpack2Smoother.hpp)
  LIST(REMOVE_ITEM SOURCES ${DIR}/Utils/ExplicitInstantiation/MueLu_Ifpack2Smoother.cpp)
ENDIF()


#
# Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  muelu
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
# touch CMakeLists.txt because a new file was created in Utils/ExplicitInstantiation of Utils/ForwardDeclaration
# touch CMakeLists.txt because a new file was created in Utils/ExplicitInstantiation of Utils/ForwardDeclaration
# touch CMakeLists.txt because a new file was created in Utils/ExplicitInstantiation of Utils/ForwardDeclaration
# touch CMakeLists.txt because a new file was created in Utils/ExplicitInstantiation of Utils/ForwardDeclaration
# touch CMakeLists.txt because a new file was created in Utils/ExplicitInstantiation of Utils/ForwardDeclaration
# touch CMakeLists.txt because a new file was created in Utils/ExplicitInstantiation of Utils/ForwardDeclaration
# touch CMakeLists.txt because a new file was created in Utils/ExplicitInstantiation of Utils/ForwardDeclaration
# touch CMakeLists.txt because a new file was created in Utils/ExplicitInstantiation of Utils/ForwardDeclaration
# touch CMakeLists.txt because a new file was created in Utils/ExplicitInstantiation of Utils/ForwardDeclaration
# touch CMakeLists.txt because a new file was created in Utils/ExplicitInstantiation of Utils/ForwardDeclaration
# touch CMakeLists.txt because a new file was created in Utils/ExplicitInstantiation of Utils/ForwardDeclaration
