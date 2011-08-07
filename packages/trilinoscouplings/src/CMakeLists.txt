
INCLUDE(PackageLibraryMacros)

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

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

#
# Core files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  TrilinosCouplings_ConfigDefs.h
  TrilinosCouplings_Version.h
  )

APPEND_SET(SOURCES
  dummy.cpp
  )

IF(HAVE_TRILINOSCOUPLINGS_ML)
  #Have Nox and ML
  IF(HAVE_TRILINOSCOUPLINGS_NOX)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ml/NonlinML)
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
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ifpack)
    APPEND_SET(HEADERS
      ifpack/Ifpack_ML.h
      )

    APPEND_SET(SOURCES
      ifpack/Ifpack_ML.cpp
      )
  ENDIF()
ENDIF()

IF(HAVE_TRILINOSCOUPLINGS_EPETRAEXT)
  #Hae Isorropia and Epetra ext
  IF(HAVE_TRILINOSCOUPLINGS_ISORROPIA)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/epetraext)
    APPEND_SET(HEADERS
      epetraext/EpetraExt_Isorropia_CrsGraph.h
      )

    APPEND_SET(SOURCES
      epetraext/EpetraExt_Isorropia_CrsGraph.cpp
      )
  ENDIF()

  #Have Amesos and Epetra ext
  IF(HAVE_TRILINOSCOUPLINGS_AMESOS)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/epetraext)
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
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/epetraext)
    APPEND_SET(HEADERS
      epetraext/EpetraExt_AmesosBTFGlobal_LinearProblem.h
      )

    APPEND_SET(SOURCES
      epetraext/EpetraExt_AmesosBTFGlobal_LinearProblem.cpp
      )
  ENDIF()
ENDIF()

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  trilinoscouplings
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
