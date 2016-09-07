

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(Intrepid2_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/Intrepid2_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# All headers
#

APPEND_SET(HEADERS
  Intrepid2_ConfigDefs.hpp
  )

#
# Shared headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Shared)
APPEND_SET(HEADERS
  Shared/Intrepid2_ArrayToolsDefCloneScale.hpp
  Shared/Intrepid2_ArrayToolsDefContractions.hpp
  Shared/Intrepid2_ArrayToolsDefDot.hpp
  Shared/Intrepid2_ArrayToolsDefScalar.hpp
  Shared/Intrepid2_ArrayToolsDefTensor.hpp
  Shared/Intrepid2_ArrayTools.hpp
  Shared/Intrepid2_BurkardtRules.hpp
  Shared/Intrepid2_BurkardtRulesDef.hpp
  Shared/Intrepid2_FieldContainerDef.hpp
  Shared/Intrepid2_FieldContainer.hpp
  Shared/Intrepid2_PointToolsDef.hpp
  Shared/Intrepid2_PointTools.hpp
  Shared/Intrepid2_PolylibDef.hpp
  Shared/Intrepid2_Polylib.hpp
  Shared/Intrepid2_RealSpaceToolsDef.hpp
  Shared/Intrepid2_RealSpaceTools.hpp
  Shared/Intrepid2_Types.hpp
  Shared/Intrepid2_Utils.hpp
  )

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Kokkos)
  APPEND_SET(HEADERS ${HEADERS}
  Kokkos/Intrepid2_KokkosRank.hpp
  Kokkos/Intrepid2_FieldContainer_Kokkos.hpp
  Kokkos/Intrepid2_FieldContainer_Kokkos_CUDA_Left.hpp
  Kokkos/Intrepid2_FieldContainer_Kokkos_CUDA_Right.hpp
  Kokkos/Intrepid2_FieldContainer_Kokkos_OpenMP_Left.hpp
  Kokkos/Intrepid2_FieldContainer_Kokkos_OpenMP_Right.hpp 
  Kokkos/Intrepid2_FieldContainer_Kokkos_PThreads_Left.hpp
  Kokkos/Intrepid2_FieldContainer_Kokkos_PThreads_Right.hpp
  Kokkos/Intrepid2_FieldContainer_Kokkos_Serial_Left.hpp
  Kokkos/Intrepid2_FieldContainer_Kokkos_Serial_Right.hpp
 )



#
# MiniTensor headers
#
IF (Intrepid2_ENABLE_Boost)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Shared/MiniTensor)
  APPEND_SET(HEADERS ${HEADERS}
    Shared/MiniTensor/Intrepid2_MiniTensor.h
    Shared/MiniTensor/Intrepid2_MiniTensor_TensorBase.h
    Shared/MiniTensor/Intrepid2_MiniTensor_TensorBase.i.h
    Shared/MiniTensor/Intrepid2_MiniTensor_TensorBase.t.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Definitions.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Geometry.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Geometry.i.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Geometry.t.h
    Shared/MiniTensor/Intrepid2_MiniTensor_LinearAlgebra.h
    Shared/MiniTensor/Intrepid2_MiniTensor_LinearAlgebra.i.h
    Shared/MiniTensor/Intrepid2_MiniTensor_LinearAlgebra.t.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Mechanics.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Mechanics.i.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Mechanics.t.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Matrix.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Matrix.i.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Matrix.t.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Solvers.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Solvers.t.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Storage.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Storage.i.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Tensor.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Tensor.i.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Tensor.t.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Tensor3.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Tensor3.i.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Tensor3.t.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Tensor4.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Tensor4.i.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Tensor4.t.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Utilities.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Utilities.i.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Vector.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Vector.i.h
    Shared/MiniTensor/Intrepid2_MiniTensor_Vector.t.h
  )
ENDIF()

#
# Discretization headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Discretization/Integration)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Discretization/FunctionSpaceTools)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Discretization/Basis)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Discretization/TensorProductSpaceTools)

APPEND_SET(HEADERS  
  Discretization/Integration/Intrepid2_AdaptiveSparseGridInterfaceDef.hpp
  Discretization/Integration/Intrepid2_AdaptiveSparseGridInterface.hpp
  Discretization/Integration/Intrepid2_AdaptiveSparseGridDef.hpp
  Discretization/Integration/Intrepid2_AdaptiveSparseGrid.hpp
  Discretization/Integration/Intrepid2_Cubature.hpp
  Discretization/Integration/Intrepid2_CubatureControlVolumeDef.hpp
  Discretization/Integration/Intrepid2_CubatureControlVolume.hpp
  Discretization/Integration/Intrepid2_CubatureControlVolumeSideDef.hpp
  Discretization/Integration/Intrepid2_CubatureControlVolumeSide.hpp
  Discretization/Integration/Intrepid2_CubatureControlVolumeBoundaryDef.hpp
  Discretization/Integration/Intrepid2_CubatureControlVolumeBoundary.hpp
  Discretization/Integration/Intrepid2_CubatureDirectDef.hpp
  Discretization/Integration/Intrepid2_CubatureDirect.hpp
  Discretization/Integration/Intrepid2_CubatureDirectLineGaussDef.hpp
  Discretization/Integration/Intrepid2_CubatureDirectLineGauss.hpp
  Discretization/Integration/Intrepid2_CubatureDirectLineGaussJacobi20.hpp
  Discretization/Integration/Intrepid2_CubatureDirectLineGaussJacobi20Def.hpp
  Discretization/Integration/Intrepid2_CubatureDirectTetDefaultDef.hpp
  Discretization/Integration/Intrepid2_CubatureDirectTetDefault.hpp
  Discretization/Integration/Intrepid2_CubatureDirectTriDefaultDef.hpp
  Discretization/Integration/Intrepid2_CubatureDirectTriDefault.hpp
  Discretization/Integration/Intrepid2_CubatureCompositeTetDef.hpp
  Discretization/Integration/Intrepid2_CubatureCompositeTet.hpp
  Discretization/Integration/Intrepid2_CubatureGenSparseDef.hpp
  Discretization/Integration/Intrepid2_CubatureGenSparse.hpp
  Discretization/Integration/Intrepid2_CubatureLineSortedDef.hpp
  Discretization/Integration/Intrepid2_CubatureLineSorted.hpp
  Discretization/Integration/Intrepid2_CubatureSparseDef.hpp
  Discretization/Integration/Intrepid2_CubatureSparse.hpp
  Discretization/Integration/Intrepid2_CubatureSparseHelper.hpp
  Discretization/Integration/Intrepid2_CubatureTensorDef.hpp
  Discretization/Integration/Intrepid2_CubatureTensor.hpp
  Discretization/Integration/Intrepid2_CubatureTensorPyr.hpp
  Discretization/Integration/Intrepid2_CubatureTensorPyrDef.hpp
  Discretization/Integration/Intrepid2_CubatureTensorSortedDef.hpp
  Discretization/Integration/Intrepid2_CubatureTensorSorted.hpp
  Discretization/Integration/Intrepid2_CubaturePolygonDef.hpp
  Discretization/Integration/Intrepid2_CubaturePolygon.hpp
  Discretization/Integration/Intrepid2_CubaturePolylib.hpp
  Discretization/Integration/Intrepid2_CubaturePolylibDef.hpp
  Discretization/Integration/Intrepid2_DefaultCubatureFactoryDef.hpp
  Discretization/Integration/Intrepid2_DefaultCubatureFactory.hpp
  Discretization/FunctionSpaceTools/Intrepid2_FunctionSpaceToolsDef.hpp
  Discretization/FunctionSpaceTools/Intrepid2_FunctionSpaceTools.hpp
  Discretization/FunctionSpaceTools/Intrepid2_FunctionSpaceToolsInPlaceDef.hpp
  Discretization/FunctionSpaceTools/Intrepid2_FunctionSpaceToolsInPlace.hpp
  Discretization/Basis/Intrepid2_BasisDef.hpp
  Discretization/Basis/Intrepid2_Basis.hpp
  Discretization/Basis/Intrepid2_HGRAD_HEX_C1_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_HEX_C1_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_HEX_C2_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_HEX_C2_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_HEX_I2_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_HEX_I2_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_HEX_Cn_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_HEX_Cn_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_LINE_C1_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_LINE_C1_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_LINE_Cn_FEM_JACOBI.hpp
  Discretization/Basis/Intrepid2_HGRAD_LINE_Cn_FEM_JACOBIDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_LINE_Cn_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_LINE_Cn_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_POLY_C1_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_POLY_C1_FEMDef.hpp  
  Discretization/Basis/Intrepid2_HGRAD_QUAD_C1_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_QUAD_C1_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_QUAD_C2_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_QUAD_C2_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_QUAD_Cn_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_QUAD_Cn_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_TET_C1_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_TET_C1_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_TET_C2_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_TET_C2_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_TET_COMP12_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_TET_COMP12_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_TRI_C1_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_TRI_C1_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_TRI_C2_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_TRI_C2_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_WEDGE_C1_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_WEDGE_C1_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_WEDGE_C2_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_WEDGE_C2_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_WEDGE_I2_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_WEDGE_I2_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_PYR_C1_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_PYR_C1_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_PYR_I2_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_PYR_I2_FEMDef.hpp
  Discretization/Basis/Intrepid2_HCURL_HEX_I1_FEM.hpp
  Discretization/Basis/Intrepid2_HCURL_HEX_I1_FEMDef.hpp
  Discretization/Basis/Intrepid2_HCURL_HEX_In_FEM.hpp
  Discretization/Basis/Intrepid2_HCURL_HEX_In_FEMDef.hpp
  Discretization/Basis/Intrepid2_HCURL_TET_I1_FEM.hpp
  Discretization/Basis/Intrepid2_HCURL_TET_I1_FEMDef.hpp
  Discretization/Basis/Intrepid2_HCURL_WEDGE_I1_FEM.hpp
  Discretization/Basis/Intrepid2_HCURL_WEDGE_I1_FEMDef.hpp
  Discretization/Basis/Intrepid2_HCURL_TRI_I1_FEM.hpp
  Discretization/Basis/Intrepid2_HCURL_TRI_I1_FEMDef.hpp
  Discretization/Basis/Intrepid2_HCURL_QUAD_I1_FEM.hpp
  Discretization/Basis/Intrepid2_HCURL_QUAD_I1_FEMDef.hpp
  Discretization/Basis/Intrepid2_HCURL_QUAD_In_FEM.hpp
  Discretization/Basis/Intrepid2_HCURL_QUAD_In_FEMDef.hpp
  Discretization/Basis/Intrepid2_HDIV_HEX_I1_FEM.hpp
  Discretization/Basis/Intrepid2_HDIV_HEX_I1_FEMDef.hpp
  Discretization/Basis/Intrepid2_HDIV_HEX_In_FEM.hpp
  Discretization/Basis/Intrepid2_HDIV_HEX_In_FEMDef.hpp
  Discretization/Basis/Intrepid2_HDIV_TET_I1_FEM.hpp
  Discretization/Basis/Intrepid2_HDIV_TET_I1_FEMDef.hpp
  Discretization/Basis/Intrepid2_HDIV_TRI_I1_FEM.hpp
  Discretization/Basis/Intrepid2_HDIV_TRI_I1_FEMDef.hpp
  Discretization/Basis/Intrepid2_HDIV_QUAD_I1_FEM.hpp
  Discretization/Basis/Intrepid2_HDIV_QUAD_I1_FEMDef.hpp
  Discretization/Basis/Intrepid2_HDIV_QUAD_In_FEM.hpp
  Discretization/Basis/Intrepid2_HDIV_QUAD_In_FEMDef.hpp
  Discretization/Basis/Intrepid2_HDIV_WEDGE_I1_FEM.hpp
  Discretization/Basis/Intrepid2_HDIV_WEDGE_I1_FEMDef.hpp 
  Discretization/Basis/Intrepid2_HGRAD_TRI_Cn_FEM_ORTH.hpp
  Discretization/Basis/Intrepid2_HGRAD_TRI_Cn_FEM_ORTHDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_TET_Cn_FEM_ORTH.hpp
  Discretization/Basis/Intrepid2_HGRAD_TET_Cn_FEM_ORTHDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_TRI_Cn_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_TRI_Cn_FEMDef.hpp
  Discretization/Basis/Intrepid2_HGRAD_TET_Cn_FEM.hpp
  Discretization/Basis/Intrepid2_HGRAD_TET_Cn_FEMDef.hpp
  Discretization/Basis/Intrepid2_HDIV_TRI_In_FEM.hpp
  Discretization/Basis/Intrepid2_HDIV_TRI_In_FEMDef.hpp
  Discretization/Basis/Intrepid2_HDIV_TET_In_FEM.hpp
  Discretization/Basis/Intrepid2_HDIV_TET_In_FEMDef.hpp
  Discretization/Basis/Intrepid2_HCURL_TRI_In_FEM.hpp
  Discretization/Basis/Intrepid2_HCURL_TRI_In_FEMDef.hpp
  Discretization/Basis/Intrepid2_HCURL_TET_In_FEM.hpp
  Discretization/Basis/Intrepid2_HCURL_TET_In_FEMDef.hpp
  Discretization/Basis/Intrepid2_OrthogonalBases.hpp
  Discretization/Basis/Intrepid2_OrthogonalBasesDef.hpp
  Discretization/Basis/Intrepid2_ProductTopology.hpp
  Discretization/Basis/Intrepid2_ProductTopologyDef.hpp
  Discretization/Basis/Intrepid2_TensorBasis.hpp
  Discretization/TensorProductSpaceTools/Intrepid2_TensorProductSpaceTools.hpp
  Discretization/TensorProductSpaceTools/Intrepid2_TensorProductSpaceToolsDef.hpp
  )

#
# Cell headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Cell)
APPEND_SET(HEADERS
  Cell/Intrepid2_CellTools.hpp
  Cell/Intrepid2_CellToolsDef.hpp
  )

#
# Orientation headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Orientation)
APPEND_SET(HEADERS
  Orientation/Intrepid2_OrientationTools.hpp
  Orientation/Intrepid2_OrientationToolsDef.hpp
  )

#
# All sources
#

APPEND_SET(SOURCES
  Intrepid2_Version.cpp
  Shared/Intrepid2_Utils.cpp
  )


IF (TPL_ENABLE_MPI)
ENDIF()


#
# C) Define the targets for package's library(s)
# 

  TRIBITS_ADD_LIBRARY(
    intrepid2
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )
