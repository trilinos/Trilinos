
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_Config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

#
# src
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_Config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  ${PACKAGE_NAME}_ConfigDefs.hpp
   RTOpPack_version.hpp
  )

SET(SOURCES ${SOURCES}
  RTOpPack_version.cpp
  )

#
# src/interfaces
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/interfaces)

SET(HEADERS ${HEADERS}
  interfaces/RTOpPack_RTOpT.hpp
  interfaces/RTOpPack_RTOpT_decl.hpp
  interfaces/RTOpPack_RTOpT_def.hpp
  interfaces/RTOpPack_Types.hpp
  )

SET(SOURCES ${SOURCES}
  interfaces/RTOpPack_RTOpT.cpp
  )

#
# src/support
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/support)

SET(HEADERS ${HEADERS}
  support/RTOpPack_RTOpServer.hpp
  support/RTOpPack_RTOpServerDecl.hpp
  support/RTOpPack_RTOpTHelpers.hpp
  support/RTOpPack_RTOpTHelpers_decl.hpp
  support/RTOpPack_RTOpTHelpers_def.hpp
  support/RTOpPack_SPMD_apply_op.hpp
  support/RTOpPack_SPMD_apply_op_decl.hpp
  support/RTOpPack_SPMD_apply_op_def.hpp
  support/RTOpPack_SparseSubVectorT.hpp
  support/RTOpPack_TOpUnaryFuncPtr.hpp
  support/RTOpPack_TOpUnaryFuncPtrDecl.hpp
  )

SET(SOURCES ${SOURCES}
  support/RTOpPack_RTOpTHelpers.cpp
  support/RTOpPack_SPMD_apply_op.cpp
  )

SET(HEADERS ${HEADERS}
  support/RTOp_parallel_helpers.h
  )

SET(SOURCES ${SOURCES}
  support/RTOp_parallel_helpers.c
  )

#
# src/ops_lib
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ops_lib)

SET(HEADERS ${HEADERS}
  ops_lib/RTOpPack_ROpCountNanInf.hpp
  ops_lib/RTOpPack_ROpDotProd.hpp
  ops_lib/RTOpPack_ROpGetSubVector.hpp
  ops_lib/RTOpPack_ROpGetSubVector_decl.hpp
  ops_lib/RTOpPack_ROpGetSubVector_def.hpp
  ops_lib/RTOpPack_ROpMax.hpp
  ops_lib/RTOpPack_ROpMaxIndex.hpp
  ops_lib/RTOpPack_ROpMaxIndexLessThanBound.hpp
  ops_lib/RTOpPack_ROpMin.hpp
  ops_lib/RTOpPack_ROpMinIndex.hpp
  ops_lib/RTOpPack_ROpMinIndexGreaterThanBound.hpp
  ops_lib/RTOpPack_ROpNorm1.hpp
  ops_lib/RTOpPack_ROpNorm2.hpp
  ops_lib/RTOpPack_ROpNormInf.hpp
  ops_lib/RTOpPack_ROpSum.hpp
  ops_lib/RTOpPack_ROpWeightedNorm2.hpp
  ops_lib/RTOpPack_TOpAXPY.hpp
  ops_lib/RTOpPack_TOpAbs.hpp
  ops_lib/RTOpPack_TOpAddScalar.hpp
  ops_lib/RTOpPack_TOpAssignScalar.hpp
  ops_lib/RTOpPack_TOpAssignVectors.hpp
  ops_lib/RTOpPack_TOpEleWiseConjProd.hpp
  ops_lib/RTOpPack_TOpEleWiseDivide.hpp
  ops_lib/RTOpPack_TOpEleWiseProd.hpp
  ops_lib/RTOpPack_TOpEleWiseProdUpdate.hpp
  ops_lib/RTOpPack_TOpLinearCombination.hpp
  ops_lib/RTOpPack_TOpLinearCombination_def.hpp
  ops_lib/RTOpPack_TOpRandomize.hpp
  ops_lib/RTOpPack_TOpReciprocal.hpp
  ops_lib/RTOpPack_TOpScaleVector.hpp
  ops_lib/RTOpPack_TOpSetSubVector.hpp
  ops_lib/RTOpPack_TOpSetSubVector_def.hpp
  )

SET(SOURCES ${SOURCES}
  ops_lib/RTOpPack_ROpGetSubVector.cpp
  ops_lib/RTOpPack_TOpLinearCombination.cpp
  ops_lib/RTOpPack_TOpSetSubVector.cpp
  )

#
# src/lapack
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/lapack)

SET(HEADERS ${HEADERS}
  lapack/RTOpPack_LapackWrappers.hpp
  )

SET(SOURCES ${SOURCES}
  lapack/RTOpPack_LapackWrappers.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  rtop
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
