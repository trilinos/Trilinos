// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   BasisEquivalenceTests_Quad.cpp
    \brief  Tests to verify that quad bases that span the same space are equivalent.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_TestUtils.hpp"

#include "BasisEquivalenceHelpers.hpp"

using namespace Intrepid2;

namespace
{
  TEUCHOS_UNIT_TEST( BasisEquivalence, QuadrilateralHierarchicalDGVersusHierarchicalCG_HGRAD )
  {
    using CGBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_QUAD;
    using DGBasis = DGHierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_QUAD;
    
    // NOTE: for the moment, OPERATOR_Dn for n > 2 not supported by DerivedBasis.  We can support more by either increasing
    //       Parameters::MaxVectorComponents (which is 7 right now), or by changing VectorData to allow a dynamic number of
    //       components.  (We were doing the latter using Kokkos::vector, but have switched to a Kokkos::Array instead to
    //       avoid using UVM.)
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2}; //, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-12; // 2e-13 is sharp on development setup for polyOrder=3; relaxing for potential architectural differences
    const double absTol=1e-11; // 5e-13 is sharp on development setup for polyOrder=3; relaxing for potential architectural differences
    
    for (int polyOrder=1; polyOrder<4; polyOrder++)
    {
      CGBasis cgBasis(polyOrder);
      DGBasis dgBasis(polyOrder);
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(cgBasis, dgBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, QuadrilateralNodalVersusHierarchicalCG_HGRAD )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_QUAD;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HGRAD_QUAD;
    
    // NOTE: for the moment, OPERATOR_Dn for n > 2 not supported by DerivedBasis.  We can support more by either increasing
    //       Parameters::MaxVectorComponents (which is 7 right now), or by changing VectorData to allow a dynamic number of
    //       components.  (We were doing the latter using Kokkos::vector, but have switched to a Kokkos::Array instead to
    //       avoid using UVM.)
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2}; //, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-13; // 8e-14 is sharp on development setup for polyOrder=3; relaxing for potential architectural differences
    const double absTol=1e-12; // 2e-13 is sharp on development setup for polyOrder=3; relaxing for potential architectural differences
    
    for (int polyOrder=1; polyOrder<4; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }
} // namespace
