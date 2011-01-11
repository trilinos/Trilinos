INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.hpp
  )

#
# B.1) utilities directory
#
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/utilities)
APPEND_SET(HEADERS utilities/Piro_ValidPiroParameters.hpp
            utilities/Piro_Version.hpp
            utilities/Piro_ConfigDefs.hpp)
APPEND_SET(SOURCES utilities/Piro_ValidPiroParameters.cpp
            utilities/Piro_Version.cpp)

#
# B.2) epetra directory
#
# Optional NOX Dependency (which includes LOCA)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/epetra)
IF (Piro_ENABLE_NOX) 
  APPEND_SET(HEADERS epetra/Piro_Epetra_NOXSolver.hpp
                     epetra/Piro_Epetra_LOCASolver.hpp
                     epetra/Piro_Epetra_VelocityVerletSolver.hpp
                     epetra/Piro_Epetra_MatrixFreeDecorator.hpp)
  APPEND_SET(SOURCES epetra/Piro_Epetra_NOXSolver.cpp
                     epetra/Piro_Epetra_LOCASolver.cpp
                     epetra/Piro_Epetra_VelocityVerletSolver.cpp
                     epetra/Piro_Epetra_MatrixFreeDecorator.cpp)
ENDIF()
# Optional Rythmos Dependency
IF (Piro_ENABLE_Rythmos) 
  APPEND_SET(HEADERS epetra/Piro_Epetra_RythmosSolver.hpp
                     thyra/Piro_Thyra_RythmosSolver.hpp
                     thyra/Piro_Thyra_RythmosSolver_Def.hpp
                     epetra/Piro_Epetra_InvertMassMatrixDecorator.hpp) 
  APPEND_SET(SOURCES epetra/Piro_Epetra_RythmosSolver.cpp
                     epetra/Piro_Epetra_InvertMassMatrixDecorator.cpp)
ENDIF()
# Optional StochasticGalerkin capability, dependion upon Stokhos and NOX
# NEED TO CREATE  epetra/Piro_Epetra_StokhosNOXSolver from ENAT
IF (Piro_ENABLE_Stokhos AND Piro_ENABLE_NOX) 
  APPEND_SET(HEADERS epetra/Piro_Epetra_StokhosNOXObserver.hpp) 
  APPEND_SET(SOURCES epetra/Piro_Epetra_StokhosNOXObserver.cpp)
ENDIF()
#

#
# B.2) thyra directory
#

IF (Piro_ENABLE_TriKota)
  ADD_DEFINITIONS(-DHAVE_CONFIG_H)
ENDIF() 

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/thyra)
  APPEND_SET(HEADERS thyra/Piro_Thyra_PerformAnalysis.hpp) 
  APPEND_SET(SOURCES thyra/Piro_Thyra_PerformAnalysis.cpp)

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  piro
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

