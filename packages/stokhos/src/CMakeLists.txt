
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

SET(HAVE_STOKHOS_DAKOTA ${${PACKAGE_NAME}_ENABLE_TriKota})

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# Add path to Dakota headers: TODO: make conditional on TriKota?
IF (${PACKAGE_NAME}_ENABLE_TriKota)
INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/../../TriKota/Dakota/packages/quadrature/rules
    ${CMAKE_CURRENT_SOURCE_DIR}/../../TriKota/Dakota/packages/quadrature/sparse_grid
  )
ENDIF()


SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  Stokhos.hpp
  Stokhos_ConfigDefs.h
  Stokhos_Version.hpp
  Stokhos_gaussq.h
  Stokhos_DynamicArrayTraits.hpp
  Stokhos_OneDOrthogPolyBasis.hpp
  Stokhos_RecurrenceBasis.hpp
  Stokhos_RecurrenceBasisImp.hpp
  Stokhos_HermiteBasis.hpp
  Stokhos_HermiteBasisImp.hpp
  Stokhos_DiscretizedStieltjesBasis.hpp
  Stokhos_DiscretizedStieltjesBasisImp.hpp
  Stokhos_RysBasis.hpp
  Stokhos_LegendreBasis.hpp
  Stokhos_LegendreBasisImp.hpp
  Stokhos_OrthogPolyBasis.hpp
  Stokhos_ProductBasis.hpp
  Stokhos_DerivBasis.hpp
  Stokhos_CompletePolynomialBasis.hpp
  Stokhos_CompletePolynomialBasisImp.hpp
  Stokhos_Sparse3Tensor.hpp
  Stokhos_Sparse3TensorImp.hpp
  Stokhos_Dense3Tensor.hpp
  Stokhos_Dense3TensorImp.hpp
  Stokhos_Quadrature.hpp
  Stokhos_TensorProductQuadrature.hpp
  Stokhos_TensorProductQuadratureImp.hpp
  Stokhos_UserDefinedQuadrature.hpp
  Stokhos_UserDefinedQuadratureImp.hpp
  Stokhos_OrthogPolyApprox.hpp
  Stokhos_OrthogPolyApproxImp.hpp
  Stokhos_OrthogPolyExpansion.hpp
  Stokhos_ConstantOrthogPolyExpansion.hpp
  Stokhos_ConstantOrthogPolyExpansionImp.hpp
  Stokhos_QuadOrthogPolyExpansion.hpp
  Stokhos_QuadOrthogPolyExpansionImp.hpp
  Stokhos_DerivOrthogPolyExpansion.hpp
  Stokhos_DerivOrthogPolyExpansionImp.hpp
  Stokhos_ForUQTKOrthogPolyExpansion.hpp
  Stokhos_ForUQTKOrthogPolyExpansionImp.hpp
  Stokhos_MatrixFreeEpetraOp.hpp
  Stokhos_KLMatrixFreeEpetraOp.hpp
  Stokhos_MeanEpetraOp.hpp
  Stokhos_SGModelEvaluator.hpp
  Stokhos_SGQuadModelEvaluator.hpp
  Stokhos_VectorOrthogPoly.hpp
  Stokhos_VectorOrthogPolyImp.hpp
  Stokhos_VectorOrthogPolyTraitsEpetra.hpp
  )

SET(SOURCES ${SOURCES}
  Stokhos_Version.cpp
  Stokhos_MatrixFreeEpetraOp.cpp
  Stokhos_KLMatrixFreeEpetraOp.cpp
  Stokhos_MeanEpetraOp.cpp
  Stokhos_SGModelEvaluator.cpp
  Stokhos_SGQuadModelEvaluator.cpp
  )

IF(Stokhos_ENABLE_TriKota)
  SET(HEADERS ${HEADERS}
    Stokhos_SparseGridQuadrature.hpp
    Stokhos_SparseGridQuadratureImp.hpp
    )
  SET(SOURCES ${SOURCES}
    Stokhos_SparseGridQuadrature.cpp
    )
ENDIF()

IF(Stokhos_ENABLE_ML)
  SET(HEADERS ${HEADERS}
    Stokhos_MLPrecOp.hpp
    )
  SET(SOURCES ${SOURCES}
    Stokhos_MLPrecOp.cpp
    )
ENDIF()


#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  stokhos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
