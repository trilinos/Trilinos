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

/** \file   BasisCardinalityTests.cpp
    \brief  Tests to verify that basis implementations return the expected cardinality.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_Types.hpp"

#include "Intrepid2_TestUtils.hpp"

#include "Kokkos_Core.hpp"

namespace
{
  using namespace Intrepid2;

  void testQuadBasisCardinality(Intrepid2::EFunctionSpace fs, int polyOrder_x, int polyOrder_y, Teuchos::FancyOStream &out, bool &success)
  {
    shards::CellTopology cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<> >() );
    const bool defineVertexFunctions = true;
    int expectedCardinality = -1;
    switch (fs)
    {
      case Intrepid2::FUNCTION_SPACE_HVOL:
      case Intrepid2::FUNCTION_SPACE_HGRAD:
        expectedCardinality = (polyOrder_x + 1) * (polyOrder_y + 1);
        break;
      case Intrepid2::FUNCTION_SPACE_HCURL:
      case Intrepid2::FUNCTION_SPACE_HDIV:
        expectedCardinality = (polyOrder_x    ) * (polyOrder_y + 1)
        + (polyOrder_x + 1) * (polyOrder_y    );
        break;
      default:
        out << "Unsupported function space.\n";
        success = false;
        return;
    }
    auto basis = getHierarchicalBasis<defineVertexFunctions>(cellTopo, fs, polyOrder_x, polyOrder_y, -1);
    if (basis->getCardinality() != expectedCardinality)
    {
      out << "FAILURE: expected cardinality of " << expectedCardinality << " but got " << basis->getCardinality() << std::endl;
      success = false;
    }
  }
  
  void testHexBasisCardinality(Intrepid2::EFunctionSpace fs, int polyOrder_x, int polyOrder_y, int polyOrder_z, Teuchos::FancyOStream &out, bool &success)
  {
    shards::CellTopology cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<> >() );
    const bool defineVertexFunctions = true;
    int expectedCardinality = -1;
    switch (fs)
    {
      case Intrepid2::FUNCTION_SPACE_HVOL:
      case Intrepid2::FUNCTION_SPACE_HGRAD:
        expectedCardinality = (polyOrder_x + 1) * (polyOrder_y + 1) * (polyOrder_z + 1);
        break;
      case Intrepid2::FUNCTION_SPACE_HCURL:
        expectedCardinality = (polyOrder_x    ) * (polyOrder_y + 1) * (polyOrder_z + 1)
        + (polyOrder_x + 1) * (polyOrder_y    ) * (polyOrder_z + 1)
        + (polyOrder_x + 1) * (polyOrder_y + 1) * (polyOrder_z    );
        break;
      case Intrepid2::FUNCTION_SPACE_HDIV:
        expectedCardinality = (polyOrder_x + 1) * (polyOrder_y    ) * (polyOrder_z    )
        + (polyOrder_x    ) * (polyOrder_y + 1) * (polyOrder_z    )
        + (polyOrder_x    ) * (polyOrder_y    ) * (polyOrder_z + 1);
        break;
      default:
        out << "Unsupported function space.\n";
        success = false;
        return;
    }
    auto basis = getHierarchicalBasis<defineVertexFunctions>(cellTopo, fs, polyOrder_x, polyOrder_y, polyOrder_z);
    if (basis->getCardinality() != expectedCardinality)
    {
      out << "FAILURE: expected cardinality of " << expectedCardinality << " but got " << basis->getCardinality() << std::endl;
      success = false;
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisCardinality, Quadrilateral )
  {
    std::vector<Intrepid2::EFunctionSpace> functionSpaces_2D = {FUNCTION_SPACE_HGRAD,FUNCTION_SPACE_HCURL,FUNCTION_SPACE_HDIV,FUNCTION_SPACE_HVOL};
    
    const int maxDegreeForCardinalityTests = 3;

    for (auto fs : functionSpaces_2D)
    {
      const int minDegree = (fs == FUNCTION_SPACE_HVOL) ? 0 : 1;
      int spaceDim = 2;
      auto cardinalityTestCases = getBasisTestCasesUpToDegree(spaceDim, minDegree, maxDegreeForCardinalityTests, maxDegreeForCardinalityTests);
      for (auto testCase : cardinalityTestCases)
      {
        testQuadBasisCardinality(fs, testCase[0], testCase[1], out, success);
      }
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisCardinality, Hexahedron )
  {
    std::vector<Intrepid2::EFunctionSpace> functionSpaces_3D = {FUNCTION_SPACE_HGRAD,FUNCTION_SPACE_HCURL,FUNCTION_SPACE_HDIV,FUNCTION_SPACE_HVOL};

    const int maxDegreeForCardinalityTests = 3;
    
    for (auto fs : functionSpaces_3D)
    {
      const int minDegree = (fs == FUNCTION_SPACE_HVOL) ? 0 : 1;
      const int spaceDim = 3;
      auto cardinalityTestCases = getBasisTestCasesUpToDegree(spaceDim, minDegree, maxDegreeForCardinalityTests, maxDegreeForCardinalityTests, maxDegreeForCardinalityTests);
      for (auto testCase : cardinalityTestCases)
      {
        testHexBasisCardinality(fs, testCase[0], testCase[1], testCase[2], out, success);
      }
    }
  }
  
} // namespace
