// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   BasisEquivalenceTests.cpp
    \brief  Tests to verify that bases on the wedge that span the same space are equivalent.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "BasisEquivalenceHelpers.hpp"

#include "Intrepid2_HGRAD_WEDGE_C1_FEM.hpp"
#include "Intrepid2_HCURL_WEDGE_I1_FEM.hpp"
#include  "Intrepid2_HDIV_WEDGE_I1_FEM.hpp"

#include "Intrepid2_HGRAD_WEDGE_C2_FEM.hpp"

#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"
#include "Intrepid2_TestUtils.hpp"
#include "Intrepid2_Types.hpp"

using namespace Intrepid2;

namespace
{
  TEUCHOS_UNIT_TEST( BasisEquivalence, WedgeNodalC1VersusDerivedNodal_HGRAD )
  {
    using DerivedNodalBasis = DerivedNodalBasisFamily<DefaultTestDeviceType>::HGRAD_WEDGE;
    using NodalC1Basis      = Basis_HGRAD_WEDGE_C1_FEM<DefaultTestDeviceType,double,double>; // regular nodal basis family does not define wedge bases beyond second order
    
    std::vector<EOperator> opsToTest {OPERATOR_GRAD};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    const ordinal_type p = 1;
    
    DerivedNodalBasis derivedNodalBasis(p);
    NodalC1Basis      nodalBasis;
    BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, derivedNodalBasis, opsToTest, relTol, absTol, out, success);
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, WedgeNodalC1VersusHierarchical_HGRAD )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_WEDGE;
    using NodalC1Basis      = Basis_HGRAD_WEDGE_C1_FEM<DefaultTestDeviceType,double,double>; // regular nodal basis family does not define wedge bases beyond second order
    
    std::vector<EOperator> opsToTest {OPERATOR_GRAD};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    const ordinal_type p = 1;
    
    HierarchicalBasis hierarchicalBasis(p);
    NodalC1Basis      nodalBasis;
    BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, WedgeNodalI1VersusDerivedNodal_HCURL )
  {
    using DerivedNodalBasis = DerivedNodalBasisFamily<DefaultTestDeviceType>::HCURL_WEDGE;
    using NodalI1Basis      = Basis_HCURL_WEDGE_I1_FEM<DefaultTestDeviceType,double,double>; // regular nodal basis family does not define wedge bases beyond second order
    
    std::vector<EOperator> opsToTest {OPERATOR_CURL};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    const ordinal_type p = 1;
    
    DerivedNodalBasis derivedNodalBasis(p);
    NodalI1Basis      nodalBasis;
    BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, derivedNodalBasis, opsToTest, relTol, absTol, out, success);
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, WedgeNodalI1VersusHierarchical_HCURL )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HCURL_WEDGE;
    using NodalI1Basis      = Basis_HCURL_WEDGE_I1_FEM<DefaultTestDeviceType,double,double>; // regular nodal basis family does not define wedge bases beyond second order (and second order only for HGRAD)
    
    std::vector<EOperator> opsToTest {OPERATOR_CURL};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    const ordinal_type p = 1;
    
    HierarchicalBasis hierarchicalBasis(p);
    NodalI1Basis      nodalBasis;
    BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, WedgeNodalI1VersusDerivedNodal_HDIV )
  {
    using DerivedNodalBasis = DerivedNodalBasisFamily<DefaultTestDeviceType>::HDIV_WEDGE;
    using NodalI1Basis      = Basis_HDIV_WEDGE_I1_FEM<DefaultTestDeviceType,double,double>; // regular nodal basis family does not define wedge bases beyond second order
    
    std::vector<EOperator> opsToTest {OPERATOR_DIV};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    const ordinal_type p = 1;
    
    DerivedNodalBasis derivedNodalBasis(p);
    NodalI1Basis      nodalBasis;
    BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, derivedNodalBasis, opsToTest, relTol, absTol, out, success);
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, WedgeNodalI1VersusHierarchical_HDIV )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HDIV_WEDGE;
    using NodalI1Basis      = Basis_HDIV_WEDGE_I1_FEM<DefaultTestDeviceType,double,double>; // regular nodal basis family does not define wedge bases beyond second order (and second order only for HGRAD)
    
    std::vector<EOperator> opsToTest {OPERATOR_DIV};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    const ordinal_type p = 1;
    
    HierarchicalBasis hierarchicalBasis(p);
    NodalI1Basis      nodalBasis;
    BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, WedgeNodalC2VersusHierarchical_HGRAD )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_WEDGE;
    using NodalC2Basis      = Basis_HGRAD_WEDGE_C2_FEM<DefaultTestDeviceType,double,double>; // regular nodal basis family does not define wedge bases beyond second order
    
    std::vector<EOperator> opsToTest {OPERATOR_GRAD};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    const ordinal_type p = 2;
    
    HierarchicalBasis hierarchicalBasis(p);
    NodalC2Basis      nodalBasis;
    BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, WedgeNodalVersusHierarchical_HGRAD )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_WEDGE;
    using NodalBasis        = DerivedNodalBasisFamily<DefaultTestDeviceType>::HGRAD_WEDGE; // regular nodal basis family does not define wedge bases beyond second order
    
    std::vector<EOperator> opsToTest {OPERATOR_GRAD};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (ordinal_type polyOrder=1; polyOrder<5; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, WedgeNodalVersusHierarchical_HCURL )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HCURL_WEDGE;
    using NodalBasis        = DerivedNodalBasisFamily<DefaultTestDeviceType>::HCURL_WEDGE; // regular nodal basis family does not define wedge bases beyond second order
    
    std::vector<EOperator> opsToTest {OPERATOR_CURL};
    
    const double relTol=1e-11;
    const double absTol=1e-11;
    
    for (ordinal_type polyOrder=1; polyOrder<5; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, WedgeNodalVersusHierarchical_HDIV )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HDIV_WEDGE;
    using NodalBasis        = DerivedNodalBasisFamily<DefaultTestDeviceType>::HDIV_WEDGE; // regular nodal basis family does not define wedge bases beyond second order
    
    std::vector<EOperator> opsToTest {OPERATOR_DIV};
    
    const double relTol=1e-11;
    const double absTol=1e-11;
    
    for (ordinal_type polyOrder=1; polyOrder<5; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, WedgeNodalVersusHierarchical_HVOL )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HVOL_WEDGE;
    using NodalBasis        = DerivedNodalBasisFamily<DefaultTestDeviceType>::HVOL_WEDGE; // no nodal bases defined for HVOL beyond the "derived" basis families
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE};
    
    const double relTol=1e-11;
    const double absTol=1e-11;
    
    for (ordinal_type polyOrder=1; polyOrder<5; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }
} // namespace
