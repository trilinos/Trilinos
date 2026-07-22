INCLUDE(TrilinosCreateClientTemplateHeaders)

#
# A) Package-specific configuration options
#

#setting up optional package defines since the default
#names are not used

IF(HAVE_TRILINOSCOUPLINGS_EPETRAEXT)
 SET(HAVE_EPETRAEXT ON)
ENDIF()

IF(HAVE_TRILINOSCOUPLINGS_ML)
 SET(HAVE_ML ON)
ENDIF()

IF(HAVE_TRILINOSCOUPLINGS_NOX)
 SET(HAVE_NOX ON)
ENDIF()

IF(HAVE_TRILINOSCOUPLINGS_AMESOS)
 SET(HAVE_AMESOS ON)
ENDIF()

IF(HAVE_TRILINOSCOUPLINGS_IFPACK)
 SET(HAVE_IFPACK ON)
ENDIF()

IF(HAVE_TRILINOSCOUPLINGS_ISORROPIA)
 SET(HAVE_ISORROPIA ON)
ENDIF()

ASSERT_DEFINED(
  HAVE_TRILINOSCOUPLINGS_TPETRA
  HAVE_TRILINOSCOUPLINGS_ZOLTAN2CORE
)

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

#
# Core files
#

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  TrilinosCouplings_Version.h
  )

APPEND_SET(SOURCES
  dummy.cpp
  )

IF(HAVE_TRILINOSCOUPLINGS_ML)
  #Have Nox and ML
  IF(HAVE_TRILINOSCOUPLINGS_NOX)
    TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ml/NonlinML)
    APPEND_SET(HEADERS
      ml/NonlinML/nlnml_ConstrainedMultiLevelOperator.H
      ml/NonlinML/nlnml_finelevelnoxinterface.H
      ml/NonlinML/nlnml_coarselevelnoxinterface.H
      ml/NonlinML/nlnml_preconditioner.H
      ml/NonlinML/nlnml_preconditioner_utils.H
      ml/NonlinML/nlnml_linearsystem.H
      ml/NonlinML/nlnml_nonlinearlevel.H
      ml/NonlinML/nlnml_prepostoperator.H
      )

    APPEND_SET(SOURCES
      ml/NonlinML/nlnml_prepostoperator.cpp
      ml/NonlinML/nlnml_ConstrainedMultiLevelOperator.cpp
      ml/NonlinML/nlnml_linearsystem.cpp
      ml/NonlinML/nlnml_nonlinearlevel.cpp
      ml/NonlinML/nlnml_coarselevelnoxinterface.cpp
      ml/NonlinML/nlnml_preconditioner1.cpp
      ml/NonlinML/nlnml_preconditioner_utils.cpp
      )
  ENDIF()

  #Have Ifpack and ML
  IF(HAVE_TRILINOSCOUPLINGS_IFPACK)
    TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ifpack)
    APPEND_SET(HEADERS
      ifpack/Ifpack_ML.h
      )

    APPEND_SET(SOURCES
      ifpack/Ifpack_ML.cpp
      )
  ENDIF()
ENDIF()

IF(HAVE_TRILINOSCOUPLINGS_EPETRAEXT)
  #Have Isorropia and Epetra ext
  IF(HAVE_TRILINOSCOUPLINGS_ISORROPIA)
    TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/epetraext)
    APPEND_SET(HEADERS
      epetraext/EpetraExt_Isorropia_CrsGraph.h
      )

    APPEND_SET(SOURCES
      epetraext/EpetraExt_Isorropia_CrsGraph.cpp
      )
  ENDIF()

  #Have Amesos and Epetra ext
  IF(HAVE_TRILINOSCOUPLINGS_AMESOS)
    TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/epetraext)
    APPEND_SET(HEADERS
      epetraext/EpetraExt_AmesosBTF_CrsGraph.h
      epetraext/EpetraExt_AmesosBTF_CrsMatrix.h
      epetraext/EpetraExt_AmesosAMDGlobal_CrsGraph.h
      )

    APPEND_SET(SOURCES
      epetraext/EpetraExt_AmesosBTF_CrsGraph.cpp
      epetraext/EpetraExt_AmesosBTF_CrsMatrix.cpp
      epetraext/EpetraExt_AmesosAMDGlobal_CrsGraph.cpp
      )
  ENDIF()

  #Have all three Epetra ext, Amesos and Isorropia
  IF(HAVE_TRILINOSCOUPLINGS_AMESOS AND HAVE_TRILINOSCOUPLINGS_ISORROPIA)
    TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/epetraext)
    APPEND_SET(HEADERS
      epetraext/EpetraExt_AmesosBTFGlobal_LinearProblem.h
      )

    APPEND_SET(SOURCES
      epetraext/EpetraExt_AmesosBTFGlobal_LinearProblem.cpp
      )
  ENDIF()
ENDIF()

IF(HAVE_TRILINOSCOUPLINGS_TPETRA)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/tpetra)
  APPEND_GLOB(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/tpetra/*.hpp)
  APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tpetra/*.cpp)
  TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/tpetra)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

  IF(TrilinosCouplings_ENABLE_EXPLICIT_INSTANTIATION)

    IF(HAVE_TRILINOSCOUPLINGS_ZOLTAN2CORE)
      SET(TRILINOSCOUPLINGS_ETI_CLASSES
        Rebalance_LinearProblem)
    ENDIF()

    # Make sure that Tpetra actually defined these variables, even if they
    # are empty.
    ASSERT_DEFINED(Tpetra_ETI_SCALARS_NO_ORDS)
    ASSERT_DEFINED(Tpetra_ETI_LORDS)
    ASSERT_DEFINED(Tpetra_ETI_GORDS)
    ASSERT_DEFINED(Tpetra_ETI_NODES)

    TPETRA_PROCESS_ETI_TEMPLATES_SLGN(SRCS
      "TrilinosCouplings_ETI_SC_LO_GO_NT.tmpl"
      "${TRILINOSCOUPLINGS_ETI_CLASSES}"
      "${Tpetra_ETI_SCALARS_NO_ORDS}" "${Tpetra_ETI_LORDS}" "${Tpetra_ETI_GORDS}" "${Tpetra_ETI_NODES}" FALSE)
    APPEND_SET(SOURCES ${SRCS})

  ENDIF()
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  trilinoscouplings
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
