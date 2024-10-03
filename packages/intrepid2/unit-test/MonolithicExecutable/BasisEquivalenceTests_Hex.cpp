// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   BasisEquivalenceTests_Hex.cpp
    \brief  Tests to verify that hexahedron bases that span the same space are equivalent.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_TestUtils.hpp"

#include "Intrepid2_HGRAD_LINE_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"

#include <Kokkos_Core.hpp>
//#include "KokkosBlas3.hpp"

#include "BasisEquivalenceHelpers.hpp"

using namespace Intrepid2;

namespace
{
  TEUCHOS_UNIT_TEST( BasisEquivalence, HexahedronHierarchicalDGVersusHierarchicalCG_HGRAD )
  {
    using CGBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_HEX;
    using DGBasis = DGHierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_HEX;
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-11; // 2e-12 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    const double absTol=1e-11; // 9e-13 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    
    // NOTE: for the moment, OPERATOR_Dn for n > 2 not supported by DerivedBasis.  We can support more by either increasing
    //       Parameters::MaxVectorComponents (which is 7 right now), or by changing VectorData to allow a dynamic number of
    //       components.  (We were doing the latter using Kokkos::vector, but have switched to a Kokkos::Array instead to
    //       avoid using UVM.)
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2}; //, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    for (int polyOrder=1; polyOrder<3; polyOrder++)
    {
      CGBasis cgBasis(polyOrder);
      DGBasis dgBasis(polyOrder);
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(cgBasis, dgBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, HexahedronHierarchicalCGVersusHypercube3D_HGRAD )
  {
    using CGBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_HEX;
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-11;
    const double absTol=1e-11;
    
    // NOTE: for the moment, OPERATOR_Dn for n > 2 not supported by DerivedBasis.  We can support more by either increasing
    //       Parameters::MaxVectorComponents (which is 7 right now), or by changing VectorData to allow a dynamic number of
    //       components.  (We were doing the latter using Kokkos::vector, but have switched to a Kokkos::Array instead to
    //       avoid using UVM.)
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2}; //, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    for (int polyOrder=1; polyOrder<3; polyOrder++)
    {
      CGBasis cgBasis(polyOrder);
      const int spaceDim = 3;
      auto hypercubeBasis = getHypercubeBasis_HGRAD<HierarchicalBasisFamily<DefaultTestDeviceType>>(polyOrder, spaceDim);
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(cgBasis, *hypercubeBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, HexahedronNodalVersusHierarchicalCG_HGRAD )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_HEX;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HGRAD_HEX;
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-12; // 2e-13 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    const double absTol=1e-13; // 2e-14 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    
    // NOTE: for the moment, OPERATOR_Dn for n > 2 not supported by DerivedBasis.  We can support more by either increasing
    //       Parameters::MaxVectorComponents (which is 7 right now), or by changing VectorData to allow a dynamic number of
    //       components.  (We were doing the latter using Kokkos::vector, but have switched to a Kokkos::Array instead to
    //       avoid using UVM.)
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2}; //, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    for (int polyOrder=1; polyOrder<3; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, HexahedronNodalCnVersusNodalC1_HGRAD )
  {
    using CnBasis = Intrepid2::Basis_HGRAD_HEX_Cn_FEM<DefaultTestDeviceType>;
    using C1Basis = Intrepid2::Basis_HGRAD_HEX_C1_FEM<DefaultTestDeviceType>;

    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-13; // ____ is sharp on development setup for polyOrder=1; relaxing for potential architectural differences
    const double absTol=1e-13; // ____ is sharp on development setup for polyOrder=1; relaxing for potential architectural differences

    // NOTE: for the moment, OPERATOR_Dn for n > 2 on Hexahedron not supported by BasisValues.  We can support more by either increasing
    //       Parameters::MaxVectorComponents (which is 7 right now), or by changing VectorData to allow a dynamic number of
    //       components.  (We were doing the latter using Kokkos::vector, but have switched to a Kokkos::Array instead to
    //       avoid using UVM.)
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2}; //, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    const int polyOrder = 1;
    CnBasis cnBasis(polyOrder);
    C1Basis c1Basis;
    BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(cnBasis, c1Basis, opsToTest, relTol, absTol, out, success);
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, HexahedronNodalCnVersusNodalC2_HGRAD )
  {
    using CnBasis = Intrepid2::Basis_HGRAD_HEX_Cn_FEM<DefaultTestDeviceType>;
    using C2Basis = Intrepid2::Basis_HGRAD_HEX_C2_FEM<DefaultTestDeviceType>;

    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-13; // 4e-14 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    const double absTol=1e-14; // 2e-15 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences

    // C2 throws an exception for OPERATOR_D5 and OPERATOR_D6, with a message that these are unsupported.
    // I'm not sure why that is, but for that reason we don't test with OPERATOR_D5 here, as we do in other tests
    // NOTE: for the moment, OPERATOR_Dn for n > 2 on Hexahedron not supported by BasisValues.  We can support more by either increasing
    //       Parameters::MaxVectorComponents (which is 7 right now), or by changing VectorData to allow a dynamic number of
    //       components.  (We were doing the latter using Kokkos::vector, but have switched to a Kokkos::Array instead to
    //       avoid using UVM.)
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2}; //, OPERATOR_D3, OPERATOR_D4};
    const int polyOrder = 2;
    CnBasis cnBasis(polyOrder);
    C2Basis c2Basis;
    BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(cnBasis, c2Basis, opsToTest, relTol, absTol, out, success);
  }
} // namespace
