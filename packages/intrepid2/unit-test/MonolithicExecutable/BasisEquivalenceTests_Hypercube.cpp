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

/** \file   BasisEquivalenceTests_Hypercube.cpp
    \brief  Tests to verify that hypercube bases that span the same space are equivalent.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_TestUtils.hpp"

#include <Kokkos_Core.hpp>
//#include "KokkosBlas3.hpp"

#include "BasisEquivalenceHelpers.hpp"

using namespace Intrepid2;

namespace
{
  TEUCHOS_UNIT_TEST( BasisEquivalence, HypercubeNodalVersusHypercubeHierarchical_HGRAD )
  {
    const int maxSpaceDim = 7; // we only test polyDegree = 1 for spaceDim > 5, due to test performance considerations
    const int maxDegree = 2;
    
    using HierarchicalBasisFamily = HierarchicalBasisFamily<DefaultTestDeviceType>;
    using NodalBasisFamily = NodalBasisFamily<DefaultTestDeviceType>;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_GRAD};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-14;
    const double absTol=1e-14;
    
    Intrepid2::EFunctionSpace fs = FUNCTION_SPACE_HGRAD;
    
    const int minDegree = (fs == FUNCTION_SPACE_HVOL) ? 0 : 1;
    for (int polyDegree = minDegree; polyDegree <= maxDegree; polyDegree++)
    {
      out << "** polyDegree " << polyDegree << " **\n";
      for (int spaceDim = 1; spaceDim <= maxSpaceDim; spaceDim++)
      {
        if ((polyDegree > 1) && (spaceDim > 5)) continue; // skip this case in the interest of test performance
        out << "** spaceDim " << spaceDim << " **\n";
        auto hierarchicalBasis = getHypercubeBasis_HGRAD<HierarchicalBasisFamily>(polyDegree, spaceDim);
        auto nodalBasis        = getHypercubeBasis_HGRAD<NodalBasisFamily>(polyDegree, spaceDim);
        BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(*nodalBasis, *hierarchicalBasis, opsToTest, relTol, absTol, out, success);
      }
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, HypercubeLowestOrderVersusSerendipity_HGRAD )
  {
    const int maxSpaceDim = 7; // lowest-order hypercube basis should be identitical to its serendipity basis
    const int minDegree = 1;
    const int maxDegree = 1;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_GRAD};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-14;
    const double absTol=1e-14;
    
    using BasisFamily = HierarchicalBasisFamily<DefaultTestDeviceType>;
    using BasisBase = typename BasisFamily::HGRAD_LINE::BasisBase;
    
    const bool compareWithGetValuesViewPath = false; // Serendipity basis does not support the View path -- particularly for high order in high dimensions, this would imply some large allocations.
    
    for (int polyDegree = minDegree; polyDegree <= maxDegree; polyDegree++)
    {
      out << "** polyDegree " << polyDegree << " **\n";
      for (int spaceDim = 1; spaceDim <= maxSpaceDim; spaceDim++)
      {
        out << "** spaceDim " << spaceDim << " **\n";
        auto hierarchicalBasis = getHypercubeBasis_HGRAD<BasisFamily>(polyDegree, spaceDim);
        auto serendipityBasis = Teuchos::rcp(new SerendipityBasis<BasisBase>(hierarchicalBasis));
        
        BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(*hierarchicalBasis, *serendipityBasis,
                                                                             opsToTest, relTol, absTol, out, success, compareWithGetValuesViewPath);
      }
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, HypercubeHigherOrderVersusSerendipity_HGRAD )
  {
    const int maxSpaceDim = 4;
    const int minDegree = 2;
    const int maxDegree = 4;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    using BasisFamily = HierarchicalBasisFamily<DefaultTestDeviceType>;
    using BasisBase = typename BasisFamily::HGRAD_LINE::BasisBase;
    
    const bool compareWithGetValuesViewPath = false; // Serendipity basis does not support the View path -- particularly for high order in high dimensions, this would imply some large allocations.
    
    for (int polyDegree = minDegree; polyDegree <= maxDegree; polyDegree++)
    {
      out << "** polyDegree " << polyDegree << " **\n";
      for (int spaceDim = 1; spaceDim <= maxSpaceDim; spaceDim++)
      {
        out << "** spaceDim " << spaceDim << " **\n";
        auto hierarchicalBasis = getHypercubeBasis_HVOL<BasisFamily>(polyDegree, spaceDim);
        auto serendipityBasis = Teuchos::rcp(new SerendipityBasis<BasisBase>(hierarchicalBasis));
        
        BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(*serendipityBasis, *hierarchicalBasis, serendipityBasis->ordinalMap(),
                                                                             opsToTest, relTol, absTol, out, success, compareWithGetValuesViewPath);
      }
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, HypercubeHigherOrderVersusSerendipity_HVOL )
  {
    const int maxSpaceDim = 4;
    const int minDegree = 2;
    const int maxDegree = 4;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    using BasisFamily = HierarchicalBasisFamily<DefaultTestDeviceType>;
    using BasisBase = typename BasisFamily::HGRAD_LINE::BasisBase;
    
    const bool compareWithGetValuesViewPath = false; // Serendipity basis does not support the View path -- particularly for high order in high dimensions, this would imply some large allocations.
    
    for (int polyDegree = minDegree; polyDegree <= maxDegree; polyDegree++)
    {
      out << "** polyDegree " << polyDegree << " **\n";
      for (int spaceDim = 1; spaceDim <= maxSpaceDim; spaceDim++)
      {
        out << "** spaceDim " << spaceDim << " **\n";
        auto hierarchicalBasis = getHypercubeBasis_HVOL<BasisFamily>(polyDegree, spaceDim);
        auto serendipityBasis = Teuchos::rcp(new SerendipityBasis<BasisBase>(hierarchicalBasis));
        
        BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(*serendipityBasis, *hierarchicalBasis, serendipityBasis->ordinalMap(),
                                                                             opsToTest, relTol, absTol, out, success, compareWithGetValuesViewPath);
      }
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, HypercubeNodalVersusHypercubeHierarchical_HVOL )
  {
    const int maxSpaceDim = 7; // we only test polyDegree = 0,1 for spaceDim > 5, due to test performance considerations
    const int maxDegreeForCardinalityTests = 2;
    
    using HierarchicalBasisFamily = HierarchicalBasisFamily<DefaultTestDeviceType>;
    using NodalBasisFamily = NodalBasisFamily<DefaultTestDeviceType>;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher-order bases in higher dimensions)
    const double relTol=1e-10;
    const double absTol=1e-10;
    
    Intrepid2::EFunctionSpace fs = FUNCTION_SPACE_HVOL;
    
    const int minDegree = (fs == FUNCTION_SPACE_HVOL) ? 0 : 1;
    for (int polyDegree = minDegree; polyDegree <= maxDegreeForCardinalityTests; polyDegree++)
    {
      out << "** polyDegree " << polyDegree << " **\n";
      for (int spaceDim = 1; spaceDim <= maxSpaceDim; spaceDim++)
      {
        if ((polyDegree > 1) && (spaceDim > 5)) continue; // skip this case in the interest of test performance
        out << "** spaceDim " << spaceDim << " **\n";
        auto hierarchicalBasis = getHypercubeBasis_HVOL<HierarchicalBasisFamily>(polyDegree, spaceDim);
        auto nodalBasis        = getHypercubeBasis_HVOL<NodalBasisFamily>(polyDegree, spaceDim);
        BasisEquivalenceHelpers::testBasisEquivalence<DefaultTestDeviceType>(*nodalBasis, *hierarchicalBasis, opsToTest, relTol, absTol, out, success);
      }
    }
  }
} // namespace
