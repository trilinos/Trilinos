INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.hpp
  )

#
# B.1) utilities 
#
APPEND_SET(HEADERS Piro_ValidPiroParameters.hpp
                   Piro_Version.hpp
                   Piro_Epetra_VelocityVerletSolver.hpp
                   Piro_Epetra_TrapezoidRuleSolver.hpp
                   Piro_Epetra_InvertMassMatrixDecorator.hpp
		   Piro_Epetra_Factory.hpp
                   Piro_ConfigDefs.hpp)
APPEND_SET(SOURCES Piro_ValidPiroParameters.cpp
                   Piro_Epetra_VelocityVerletSolver.cpp
                   Piro_Epetra_TrapezoidRuleSolver.cpp
                   Piro_Epetra_InvertMassMatrixDecorator.cpp
		   Piro_Epetra_Factory.cpp
                   Piro_Version.cpp)

#
# B.2) epetra files
#
# Optional NOX Dependency (which includes LOCA)
IF (Piro_ENABLE_NOX) 
  APPEND_SET(HEADERS Piro_Epetra_NOXSolver.hpp
                     Piro_NOXSolver.hpp
                     Piro_NOXSolver_Def.hpp
                     Piro_Epetra_LOCASolver.hpp
                     Piro_Epetra_MatrixFreeDecorator.hpp
		     Piro_Epetra_SensitivityOperator.hpp)
  APPEND_SET(SOURCES Piro_Epetra_NOXSolver.cpp
                     Piro_Epetra_LOCASolver.cpp
                     Piro_Epetra_MatrixFreeDecorator.cpp
		     Piro_Epetra_SensitivityOperator.cpp)
ENDIF()
# Optional Rythmos Dependency
IF (Piro_ENABLE_Rythmos) 
  APPEND_SET(HEADERS Piro_Epetra_RythmosSolver.hpp
                     Piro_RythmosSolver.hpp
                     Piro_RythmosSolver_Def.hpp) 
  APPEND_SET(SOURCES Piro_Epetra_RythmosSolver.cpp)
ENDIF()

# Optional StochasticGalerkin capability, depending upon Stokhos and NOX
IF (Piro_ENABLE_Stokhos AND Piro_ENABLE_NOX) 
  APPEND_SET(HEADERS Piro_Epetra_StokhosNOXObserver.hpp
  		     Piro_Epetra_StokhosSolverFactory.hpp
  		     Piro_Epetra_StokhosSolver.hpp) 
  APPEND_SET(SOURCES Piro_Epetra_StokhosNOXObserver.cpp
  		     Piro_Epetra_StokhosSolverFactory.cpp
  		     Piro_Epetra_StokhosSolver.cpp)
ENDIF()
#

#
# B.2) Analysis Tools
#

IF (Piro_ENABLE_TriKota)
  ADD_DEFINITIONS(${TRIKOTA_DAKOTA_DEFINITIONS})
ENDIF() 

APPEND_SET(HEADERS Piro_PerformAnalysis.hpp) 
APPEND_SET(SOURCES Piro_PerformAnalysis.cpp)

## DEPRECATED Thyra Namespace in Piro: REMOVE THESE FOR TRILINOS 11.0)
APPEND_SET(HEADERS Piro_Thyra_PerformAnalysis.hpp) 
APPEND_SET(SOURCES Piro_Thyra_PerformAnalysis.cpp)

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  piro
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
