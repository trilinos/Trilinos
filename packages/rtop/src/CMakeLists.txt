
INCLUDE(TrilinosPackageLibraryHelpers)

#
# A) Package-specific configuration options
#

TRILINOS_PACKAGE_CONFIGURE_FILE(${PROJECT_NAME}_Config.h)

#
# B) Define the header and source files (and directories)
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
LINK_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

#
# src/interfaces
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/interfaces)

SET(RTOP_INTERFACES_H
  interfaces/RTOpPack_RTOpT.hpp
  interfaces/RTOpPack_RTOpT_decl.hpp
  interfaces/RTOpPack_RTOpT_def.hpp
  interfaces/RTOpPack_Types.hpp
  )

SET(RTOP_INTERFACES
  interfaces/RTOpPack_RTOpT.cpp
  )

#
# src/support
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/support)

SET(RTOP_SUPPORT_H
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

SET(RTOP_SUPPORT
  support/RTOpPack_RTOpTHelpers.cpp
  support/RTOpPack_SPMD_apply_op.cpp
  )

SET(RTOP_SUPPORT_C_H
  support/RTOp_parallel_helpers.h
  )

SET(RTOP_SUPPORT_C
  support/RTOp_parallel_helpers.c
  )

#
# src/ops_lib
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ops_lib)

SET(RTOP_OPS_LIB_H
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

SET(RTOP_OPS_LIB
  ops_lib/RTOpPack_ROpGetSubVector.cpp
  ops_lib/RTOpPack_TOpLinearCombination.cpp
  ops_lib/RTOpPack_TOpSetSubVector.cpp
  )

#
# src/lapack
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/lapack)

SET(RTOP_LAPACK_H
  lapack/RTOpPack_LapackWrappers.hpp
  )

SET(RTOP_LAPACK
  lapack/RTOpPack_LapackWrappers.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRILINOS_PACKAGE_ADD_LIBRARY(
  rtop
  HEADERS
   RTOp_ConfigDefs.hpp RTOpPack_version.hpp
    ${RTOP_INTERFACES_H}
    ${RTOP_SUPPORT_H}
    ${RTOP_SUPPORT_C_H}
    ${RTOP_OPS_LIB_H}
    ${RTOP_LAPACK_H}
  SOURCES
    RTOpPack_version.cpp
    ${RTOP_INTERFACES}
    ${RTOP_SUPPORT_C}
    ${RTOP_SUPPORT}
    ${RTOP_OPS_LIB}
    ${RTOP_LAPACK}
  )

#
# D) Export the dependancy variables of this package for other packages to use
#

TRILINOS_PACKAGE_EXPORT_DEPENDENCY_VARIABLES()
