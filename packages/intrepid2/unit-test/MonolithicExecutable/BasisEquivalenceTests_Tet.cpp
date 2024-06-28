// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   BasisEquivalenceTests_Tet.cpp
    \brief  Tests to verify that tet bases that span the same space are equivalent.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_TestUtils.hpp"

#include "BasisEquivalenceHelpers.hpp"

using namespace Intrepid2;

namespace
{
  TEUCHOS_UNIT_TEST( BasisEquivalence, TetrahedronNodalCnVersusNodalC2_HGRAD )
  {
    using CnBasis = Intrepid2::Basis_HGRAD_TET_Cn_FEM<DefaultTestDeviceType>;
    using C2Basis = Intrepid2::Basis_HGRAD_TET_C2_FEM<DefaultTestDeviceType>;
    
    // OPERATOR_D2 and above are not supported by either the nodal or the hierarchical basis at present...
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-12; // 2e-14 is sharp on development setup; relaxing for potential architectural differences
    const double absTol=1e-13; // 3e-15 is sharp on development setup; relaxing for potential architectural differences
    
    CnBasis cnBasis(2);
    C2Basis c2Basis;
    BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(cnBasis, c2Basis, opsToTest, relTol, absTol, out, success);
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, TetrahedronHierarchicalDGVersusHierarchicalCG_HGRAD )
  {
    using CGBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_TET;
    using DGBasis = DGHierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_TET;
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-6; // 3e-8 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    const double absTol=1e-9; // 5e-11 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1};
    for (int polyOrder=1; polyOrder<7; polyOrder++)
    {
      CGBasis cgBasis(polyOrder);
      DGBasis dgBasis(polyOrder);
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(cgBasis, dgBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, TetrahedronNodalVersusHierarchicalCG_HGRAD )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_TET;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HGRAD_TET;
    
    // OPERATOR_D2 and above are not supported by either the nodal or the hierarchical basis at present...
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-6;  // 3e-08 is sharp on development setup for polyOrder=6; relaxing for potential architectural differences
    const double absTol=1e-10; // 3e-12 is sharp on development setup for polyOrder=6; relaxing for potential architectural differences
    
    for (int polyOrder=1; polyOrder<7; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      
//      auto cellTopo = nodalBasis.getBaseCellTopology();
//      const int faceDim = 2;
//      for (int intrepid2FaceOrdinal=0; intrepid2FaceOrdinal<4; intrepid2FaceOrdinal++)
//      {
//        int vertex0 = cellTopo->getNodeMap(faceDim, intrepid2FaceOrdinal, 0);
//        int vertex1 = cellTopo->getNodeMap(faceDim, intrepid2FaceOrdinal, 1);
//        int vertex2 = cellTopo->getNodeMap(faceDim, intrepid2FaceOrdinal, 2);
//        std::cout << "face " << intrepid2FaceOrdinal << ": " << vertex0 << vertex1 << vertex2 << std::endl;
//      }
      
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, TetrahedronNodalVersusHierarchical_HCURL )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HCURL_TET;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HCURL_TET;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_CURL};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-6;
    const double absTol=1e-10;
    
    for (int polyOrder=1; polyOrder<7; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, TetrahedronNodalVersusHierarchical_HDIV )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HDIV_TET;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HDIV_TET;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_DIV};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-6;
    const double absTol=1e-9;
    
    for (int polyOrder=1; polyOrder<7; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, TetrahedronNodalVersusHierarchical_HVOL )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HVOL_TET;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HVOL_TET;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-6;
    const double absTol=1e-10;
    
    for (int polyOrder=1; polyOrder<7; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }
} // namespace
