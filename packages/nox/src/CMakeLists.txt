
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_Config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(NOINSTALLHEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_Config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Base headers
APPEND_SET(HEADERS
  NOX.H
  NOX_Common.H
  ) 

# Abstract
APPEND_SET(HEADERS
  NOX_Abstract_Vector.H
  NOX_Abstract_MultiVector.H
  NOX_Abstract_Group.H
  NOX_Abstract_PrePostOperator.H
  )

APPEND_SET(SOURCES
  NOX_Abstract_Vector.C
  NOX_Abstract_Group.C
  )

# Linesearch
APPEND_SET(HEADERS
  NOX_LineSearch_Generic.H
  NOX_LineSearch_Factory.H
  NOX_LineSearch_UserDefinedFactory.H
  NOX_LineSearch_UserDefinedFactoryT.H
  NOX_LineSearch_Utils_Slope.H
  NOX_LineSearch_Utils_Printing.H
  NOX_LineSearch_Utils_Counters.H
  NOX_LineSearch_FullStep.H
  NOX_LineSearch_Backtrack.H
  NOX_LineSearch_MoreThuente.H
  NOX_LineSearch_NonlinearCG.H
  NOX_LineSearch_Polynomial.H 
)
APPEND_SET(SOURCES
  NOX_LineSearch_Factory.C
  NOX_LineSearch_Utils_Slope.C
  NOX_LineSearch_Utils_Printing.C
  NOX_LineSearch_Utils_Counters.C
  NOX_LineSearch_FullStep.C
  NOX_LineSearch_Backtrack.C
  NOX_LineSearch_MoreThuente.C
  NOX_LineSearch_NonlinearCG.C
  NOX_LineSearch_Polynomial.C 
)

# Direction
APPEND_SET(HEADERS
  NOX_Direction_Generic.H
  NOX_Direction_Factory.H
  NOX_Direction_UserDefinedFactory.H
  NOX_Direction_UserDefinedFactoryT.H
  NOX_Direction_Newton.H
  NOX_Direction_Broyden.H
  NOX_Direction_NonlinearCG.H
  NOX_Direction_SteepestDescent.H
  NOX_Direction_Utils_InexactNewton.H
  )
APPEND_SET(SOURCES
  NOX_Direction_Generic.C
  NOX_Direction_Factory.C
  NOX_Direction_Newton.C
  NOX_Direction_Broyden.C
  NOX_Direction_NonlinearCG.C
  NOX_Direction_SteepestDescent.C
  NOX_Direction_Utils_InexactNewton.C
  )

# Solvers
APPEND_SET(HEADERS
  NOX_Solver_Generic.H
  NOX_Solver_Factory.H
  NOX_Solver_SolverUtils.H
  NOX_Solver_LineSearchBased.H
  NOX_Solver_TrustRegionBased.H
  NOX_Solver_InexactTrustRegionBased.H
  NOX_Solver_TensorBased.H
  NOX_Solver_PrePostOperator.H
  )
APPEND_SET(SOURCES
  NOX_Solver_Factory.C
  NOX_Solver_SolverUtils.C
  NOX_Solver_LineSearchBased.C
  NOX_Solver_TrustRegionBased.C
  NOX_Solver_InexactTrustRegionBased.C
  NOX_Solver_TensorBased.C
  NOX_Solver_PrePostOperator.C
  )

# Merit Function
APPEND_SET(HEADERS
  NOX_MeritFunction_Generic.H
  NOX_MeritFunction_SumOfSquares.H
  )
APPEND_SET(SOURCES
  NOX_MeritFunction_SumOfSquares.C
  )

# Status Tests
APPEND_SET(HEADERS
  NOX_StatusTest_Generic.H
  NOX_StatusTest_Factory.H
  NOX_StatusTest_NormF.H
  NOX_StatusTest_FiniteValue.H
  NOX_StatusTest_NormUpdate.H
  NOX_StatusTest_Combo.H
  NOX_StatusTest_MaxIters.H
  NOX_StatusTest_NormWRMS.H
  NOX_StatusTest_Stagnation.H
  NOX_StatusTest_Divergence.H 
  )
APPEND_SET(SOURCES
  NOX_StatusTest_Generic.C
  NOX_StatusTest_Factory.C
  NOX_StatusTest_NormF.C
  NOX_StatusTest_FiniteValue.C
  NOX_StatusTest_NormUpdate.C
  NOX_StatusTest_Combo.C
  NOX_StatusTest_MaxIters.C
  NOX_StatusTest_NormWRMS.C
  NOX_StatusTest_Stagnation.C
  NOX_StatusTest_Divergence.C 
  )

# Utils
APPEND_SET(HEADERS
  NOX_Utils.H
  NOX_SharedObjectTemplate.H
  NOX_Random.H
  NOX_Version.H
  NOX_GlobalData.H
  )
APPEND_SET(SOURCES
  NOX_Utils.C
  NOX_Random.C
  NOX_Version.C
  NOX_GlobalData.C
  )

# Multivector
APPEND_SET(HEADERS
  NOX_MultiVector.H
  )
APPEND_SET(SOURCES
  NOX_MultiVector.C
  )

# Multi-physics
APPEND_SET(HEADERS
  NOX_Multiphysics_Group.H
  NOX_Multiphysics_Solver_Manager.H
  NOX_Multiphysics_Solver_Generic.H
  NOX_Multiphysics_Solver_FixedPointBased.H
  NOX_Multiphysics_DataExchange_Interface.H
  )
APPEND_SET(SOURCES
  NOX_Multiphysics_Group.C
  NOX_Multiphysics_Solver_Manager.C
  NOX_Multiphysics_Solver_FixedPointBased.C
  )

# Prerelease Code
IF (${PACKAGE_NAME}_BUILD_PRERELEASE)
  
  APPEND_SET(HEADERS
    NOX_Solver_TensorBasedTest.H
    NOX_LineSearch_Tensor.H
    NOX_Direction_Tensor.H
    NOX_Direction_QuasiNewton.H
    NOX_Direction_ModifiedNewton.H
    )

  APPEND_SET(SOURCES
    NOX_Solver_TensorBasedTest.C
    NOX_LineSearch_Tensor.C
    NOX_Direction_Tensor.C
    NOX_Direction_QuasiNewton.C
    NOX_Direction_ModifiedNewton.C
    )

ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  nox
  HEADERS ${HEADERS}
  NOINSTALLHEADERS ${NOINSTALLHEADERS}
  SOURCES ${SOURCES}
  )
