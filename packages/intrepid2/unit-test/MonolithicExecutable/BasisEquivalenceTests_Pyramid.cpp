// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   BasisEquivalenceTests_Pyramid.cpp
    \brief  Tests to verify that bases on the pyramid that span the same space are equivalent.  At present, we only have such redundancy between the legacy first-order HGRAD basis and the first-order hierarchical basis.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "BasisEquivalenceHelpers.hpp"

#include "Intrepid2_HGRAD_PYR_C1_FEM.hpp"

using namespace Intrepid2;

namespace
{
  TEUCHOS_UNIT_TEST( BasisEquivalence, PyramidNodalC1VersusHierarchical_HGRAD )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_PYR;
    using NodalC1Basis      = Basis_HGRAD_PYR_C1_FEM<DefaultTestDeviceType,double,double>; // regular nodal basis family does not define pyramid bases; we don't have these implemented except for the first-order (lowest-order) basis.
    
    std::vector<EOperator> opsToTest {OPERATOR_GRAD};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    const ordinal_type p = 1;
    
    HierarchicalBasis hierarchicalBasis(p);
    NodalC1Basis      nodalBasis;
    BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
  }
} // namespace
