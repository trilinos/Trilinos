
INCLUDE(TrilinosPackageLibraryHelpers)

#
# A) Package-specific configuration options
#

TRILINOS_PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  Stokhos_ConfigDefs.h
  Stokhos_Version.hpp
  Stokhos_DynamicArrayTraits.hpp
  Stokhos_Polynomial.hpp
  Stokhos_PolynomialImp.hpp
  Stokhos_OrthogPolyBasis.hpp
  Stokhos_OrthogPolyBasisBase.hpp
  Stokhos_OrthogPolyBasisBaseImp.hpp
  Stokhos_HermiteBasis.hpp
  Stokhos_HermiteBasisImp.hpp
  Stokhos_HermiteEBasis.hpp
  Stokhos_HermiteEBasisImp.hpp
  Stokhos_HermiteEBasis2.hpp
  Stokhos_HermiteEBasis2Imp.hpp
  Stokhos_HermiteUnitBasis.hpp
  Stokhos_HermiteUnitBasisImp.hpp
  Stokhos_LegendreBasis.hpp
  Stokhos_LegendreBasisImp.hpp
  Stokhos_CompletePolynomialBasis.hpp
  Stokhos_CompletePolynomialBasisImp.hpp
  Stokhos_TripleProduct.hpp
  Stokhos_TripleProductImp.hpp
  Stokhos_TayTripleProduct.hpp
  Stokhos_TayTripleProductImp.hpp
  Stokhos_OrthogPolyApprox.hpp
  Stokhos_OrthogPolyApproxImp.hpp
  Stokhos_OrthogPolyExpansion.hpp
  Stokhos_OrthogPolyExpansionImp.hpp
  Stokhos_TayOrthogPolyExpansion.hpp
  Stokhos_TayOrthogPolyExpansionImp.hpp
  )

SET(SOURCES ${SOURCES}
  Stokhos_Version.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRILINOS_PACKAGE_ADD_LIBRARY(
  stokhos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

#
# D) Export the dependancy variables of this package for other packages to use
#

TRILINOS_PACKAGE_EXPORT_DEPENDENCY_VARIABLES()
