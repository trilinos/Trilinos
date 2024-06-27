// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   BasisEquivalenceTests_Tri.cpp
    \brief  Tests to verify that triangle bases that span the same space are equivalent.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_TestUtils.hpp"

#include "BasisEquivalenceHelpers.hpp"

using namespace Intrepid2;

namespace
{
  TEUCHOS_UNIT_TEST( BasisEquivalence, TriangleNodalVersusHierarchicalCG_HGRAD )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_TRI;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HGRAD_TRI;
    
    // OPERATOR_D2 and above are not supported by either the nodal or the hierarchical basis at present...
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_GRAD, OPERATOR_D1};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-11; // 7e-13 is sharp on development setup for polyOrder=4; relaxing for potential architectural differences
    const double absTol=1e-12; // 2e-14 is sharp on development setup for polyOrder=4; relaxing for potential architectural differences
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, TriangleHierarchicalCGVersusHierarchicalDG_HGRAD )
  {
    using CGBasis =   HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_TRI;
    using DGBasis = DGHierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_TRI;
    
    // OPERATOR_D2 and above are not supported by either the nodal or the hierarchical basis at present...
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_GRAD, OPERATOR_D1};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-11; // 7e-13 is sharp on development setup for polyOrder=4; relaxing for potential architectural differences
    const double absTol=1e-12; // 2e-14 is sharp on development setup for polyOrder=4; relaxing for potential architectural differences
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      CGBasis cgBasis(polyOrder);
      DGBasis dgBasis(polyOrder);
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(cgBasis, dgBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, TriangleNodalVersusHierarchical_HCURL )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HCURL_TRI;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HCURL_TRI;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_CURL};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-10;
    const double absTol=1e-10;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(hierarchicalBasis, nodalBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, TriangleNodalVersusHierarchical_HDIV )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HDIV_TRI;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HDIV_TRI;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_DIV};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-10;
    const double absTol=1e-10;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(hierarchicalBasis, nodalBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, TriangleNodalVersusHierarchical_HVOL )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HVOL_TRI;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HVOL_TRI;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-11;
    const double absTol=1e-12;
    
    for (int polyOrder=0; polyOrder<5; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
} // namespace
