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

/** \file   SubBasisInclusionTests.cpp
    \brief  Tests to verify that hierarchical bases are hierarchical -- that is, that lower-order instances are subsets of higher-order instances.
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

  bool testSubBasis(shards::CellTopology cellTopo, Intrepid2::EFunctionSpace fs, const double tol, Teuchos::FancyOStream &out, bool &success,
                    int polyOrder_x, int polyOrder_y=-1, int polyOrder_z = -1, bool defineVertexFunctions = true)
  {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using OutputScalar = double;
    using PointScalar  = double;
    using namespace Intrepid2;
    using Basis = Basis<ExecSpace,OutputScalar,PointScalar>;
    
    auto vectorsMatch = [](std::vector<int> lhs,std::vector<int> rhs)
    {
      if (lhs.size() != rhs.size()) return false;
      for (unsigned i=0; i<lhs.size(); i++)
      {
        if (lhs[i] != rhs[i]) return false;
      }
      return true;
    };
    
    Teuchos::RCP<Basis> basis;
    if (defineVertexFunctions)
    {
      basis = getHierarchicalBasis<true>(cellTopo, fs, polyOrder_x, polyOrder_y, polyOrder_z);
    }
    else
    {
      basis = getHierarchicalBasis<false>(cellTopo, fs, polyOrder_x, polyOrder_y, polyOrder_z);
    }
    int spaceDim = cellTopo.getDimension();
    int minDegree = (fs == Intrepid2::FUNCTION_SPACE_HVOL) ? 0 : 1;
    
    auto subBasisDegreeTestCases = getBasisTestCasesUpToDegree(spaceDim, minDegree, polyOrder_x, polyOrder_y, polyOrder_z);
    
    std::vector<int> degrees(spaceDim);
    degrees[0] = polyOrder_x;
    if (spaceDim > 1) degrees[1] = polyOrder_y;
    if (spaceDim > 2) degrees[2] = polyOrder_z;
    //  // TODO: consider what to do when non-hypercubes are being tested
    
    int numPoints_1D = 5;
    auto inputPoints = getInputPointsView<PointScalar>(cellTopo, numPoints_1D);
    int numPoints = inputPoints.extent_int(0);
    
    auto op = Intrepid2::OPERATOR_VALUE;
    auto outputValues = getOutputView<OutputScalar>(fs, op, basis->getCardinality(), numPoints, spaceDim);
    
    basis->getValues(outputValues, inputPoints, op);
    
    for (auto testCase : subBasisDegreeTestCases)
    {
      // pad test case with -1s to get to 3D (so we can make one call to getBasis below)
      auto paddedTestCase = testCase;
      for (int d=testCase.size(); d<3; d++)
      {
        paddedTestCase.push_back(-1);
      }
      Teuchos::RCP<Basis> subBasis;
      if (defineVertexFunctions)
      {
        subBasis = getHierarchicalBasis<true>(cellTopo, fs, paddedTestCase[0], paddedTestCase[1], paddedTestCase[2]);
      }
      else
      {
        subBasis = getHierarchicalBasis<false>(cellTopo, fs, paddedTestCase[0], paddedTestCase[1], paddedTestCase[2]);
      }
      auto allLowerOrderFieldOrdinals = basis->getFieldOrdinalsForDegree(testCase);
      // for H(curl) and H(div), allLowerOrderFieldOrdinals may have some extra entries beyond what the subBasis has
      // (this has to do with the fact that the different families have different degrees in each dimension)
      // we filter according to polynomial degree of the fields, and place the result in lowerOrderFieldOrdinals
      std::vector<int> lowerOrderFieldOrdinals;
      if ((fs == Intrepid2::FUNCTION_SPACE_HCURL) || (fs == Intrepid2::FUNCTION_SPACE_HDIV))
      {
        unsigned allFieldOrdinalIndex = 0;
        for (int subBasisFieldOrdinal=0; subBasisFieldOrdinal < subBasis->getCardinality(); subBasisFieldOrdinal++)
        {
          while (   (allFieldOrdinalIndex<allLowerOrderFieldOrdinals.size())
                 && (
                     !vectorsMatch(
                                   basis   ->getPolynomialDegreeOfFieldAsVector(allLowerOrderFieldOrdinals[allFieldOrdinalIndex]),
                                   subBasis->getPolynomialDegreeOfFieldAsVector(subBasisFieldOrdinal) // same poly order in each dimension
                                   )
                     )
                 )
          {
            allFieldOrdinalIndex++;
          }
          if (allFieldOrdinalIndex < allLowerOrderFieldOrdinals.size())
          {
            // then the subcell dim and subcell ordinal must match -- add this into our filtered list
            lowerOrderFieldOrdinals.push_back(allLowerOrderFieldOrdinals[allFieldOrdinalIndex]);
          }
          // if all has gone well, we've established a correspondence between subBasisFieldOrdinal and allLowerOrderFieldOrdinals[allFieldOrdinalIndex]
          // we therefore should increment allFieldOrdinalIndex
          allFieldOrdinalIndex++;
        }
      }
      else
      {
        lowerOrderFieldOrdinals = allLowerOrderFieldOrdinals;
      }
      if (lowerOrderFieldOrdinals.size() != unsigned(subBasis->getCardinality()))
      {
        success = false;
        out << "FAILURE: for test case {";
        for (unsigned d=0; d<testCase.size(); d++)
        {
          out << testCase[d];
          if (d<testCase.size()-1) out << ",";
        }
        out << "}, expected fieldOrdinals for degree to have " << subBasis->getCardinality() << " entries, but had " << lowerOrderFieldOrdinals.size() << std::endl;
        out.flush();
        continue; // next test case
      }
      
      auto subBasisOutputValues = getOutputView<OutputScalar>(fs, op, subBasis->getCardinality(), numPoints, spaceDim);
      subBasis->getValues(subBasisOutputValues, inputPoints, op);
      
      bool vectorValued = (outputValues.rank() == 3); // F,P,D -- if scalar-valued, F,P
      
      for (int pointOrdinal=0; pointOrdinal<numPoints; pointOrdinal++)
      {
        // by construction, the sub-basis should have fields in the same order as the original basis
        // (the original basis just has extra fields interspersed)
        int subBasisFieldOrdinal = 0;
        for (int fieldOrdinal : lowerOrderFieldOrdinals)
        {
          if (!vectorValued)
          {
            double originalValue = outputValues(fieldOrdinal,pointOrdinal);
            double subBasisValue = subBasisOutputValues(subBasisFieldOrdinal,pointOrdinal);
            
            bool valuesMatch = essentiallyEqual(originalValue, subBasisValue, tol);
            
            if (!valuesMatch)
            {
              if (fs == Intrepid2::FUNCTION_SPACE_HCURL)
              {
                // scalar values are the curls
                out << "curl ";
              }
              else if (fs == Intrepid2::FUNCTION_SPACE_HDIV)
              {
                // scalar values are the div values
                out << "div ";
              }
              double x = inputPoints(pointOrdinal,0);
              double y = (spaceDim > 1) ? inputPoints(pointOrdinal,1) : -2.0;
              double z = (spaceDim > 2) ? inputPoints(pointOrdinal,2) : -2.0;
              
              if (spaceDim == 1)
                out << "values for "  << x  << " differ for field ordinal " << fieldOrdinal;
              else if (spaceDim == 2)
                out << "values for ("  << x << "," << y << ") differ for field ordinal " << fieldOrdinal;
              else
                out << "values for ("  << x << "," << y << "," << z << ") differ for field ordinal " << fieldOrdinal;
              out << ": expected " << subBasisValue << "; actual " << originalValue;
              out << " (diff: " << subBasisValue-originalValue << ")" << std::endl;
              success = false;
            }
          }
          else // vector-valued
          {
            bool valuesMatch = true;
            for (int d=0; d<spaceDim; d++)
            {
              double originalValue = outputValues(fieldOrdinal,pointOrdinal,d);
              double subBasisValue = subBasisOutputValues(subBasisFieldOrdinal,pointOrdinal,d);
              
              if (!essentiallyEqual(originalValue, subBasisValue, tol))
              {
                valuesMatch = false;
              }
            }
            
            if (!valuesMatch)
            {
              double x = inputPoints(pointOrdinal,0);
              double y = (spaceDim > 1) ? inputPoints(pointOrdinal,1) : -2.0;
              double z = (spaceDim > 2) ? inputPoints(pointOrdinal,2) : -2.0;
              
              if (spaceDim == 1)
                out << "values for "  << x  << " differ for field ordinal " << fieldOrdinal;
              else if (spaceDim == 2)
                out << "values for ("  << x << "," << y << ") differ for field ordinal " << fieldOrdinal;
              else
                out << "values for ("  << x << "," << y << "," << z << ") differ for field ordinal " << fieldOrdinal;
              out << ": expected value (lower-order basis fieldOrdinal " << subBasisFieldOrdinal << "): (";
              for (int d=0; d<spaceDim; d++)
              {
                out << subBasisOutputValues(subBasisFieldOrdinal,pointOrdinal,d);
                if (d<spaceDim-1) out << ",";
              }
              out << "); actual (larger basis fieldOrdinal " << fieldOrdinal << ") was (";
              for (int d=0; d<spaceDim; d++)
              {
                out << outputValues(fieldOrdinal,pointOrdinal,d);
                if (d<spaceDim-1) out << ",";
              }
              out << ")" << std::endl;
              success = false;
            }
          }
          subBasisFieldOrdinal++;
        }
      }
    }
    
    return success;
  }
  
  void runSubBasisTests(shards::CellTopology &cellTopo, Teuchos::FancyOStream &out, bool &success)
  {
    const double tol = TEST_TOLERANCE_TIGHT;
    const int maxDegree = 4;
    
    std::vector<Intrepid2::EFunctionSpace> functionSpaces_1D = {FUNCTION_SPACE_HGRAD,FUNCTION_SPACE_HVOL};
    std::vector<Intrepid2::EFunctionSpace> functionSpaces_2D = {FUNCTION_SPACE_HGRAD,FUNCTION_SPACE_HCURL,FUNCTION_SPACE_HDIV,FUNCTION_SPACE_HVOL};
    std::vector<Intrepid2::EFunctionSpace> functionSpaces_3D = {FUNCTION_SPACE_HGRAD,FUNCTION_SPACE_HCURL,FUNCTION_SPACE_HDIV,FUNCTION_SPACE_HVOL};
    
    std::vector<std::vector<Intrepid2::EFunctionSpace> > functionSpacesForDimension = {functionSpaces_1D,functionSpaces_2D,functionSpaces_3D};
    
    const int spaceDim = cellTopo.getDimension();
    
    auto functionSpaces = functionSpacesForDimension[spaceDim-1];
    
    for (auto fs : functionSpaces)
    {
      std::vector<bool> continuousBasisValues;
      if (fs != FUNCTION_SPACE_HVOL)
      {
        continuousBasisValues = {true,false};
      }
      else
      {
        continuousBasisValues = {true}; // false case not supported by the dof tag stuff that testSubBasis() does
      }
      for (auto continuousBasis : continuousBasisValues) // corresponds to "defineVertexFunctions" in line basis definitions
      {
        for (int degree=1; degree<=maxDegree; degree++)
        {
          testSubBasis(cellTopo, fs, tol, out, success, degree,degree,degree,continuousBasis);
        }
      }
    }
  }
  
  TEUCHOS_UNIT_TEST( SubBasisInclusion, Line )
  {
    shards::CellTopology lineTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >() );
    runSubBasisTests(lineTopo, out, success);
  }
  
  TEUCHOS_UNIT_TEST( SubBasisInclusion, Quadrilateral )
  {
    shards::CellTopology quadTopo = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<> >() );
    runSubBasisTests(quadTopo, out, success);
  }
  
  TEUCHOS_UNIT_TEST( SubBasisInclusion, Hexahedron )
  {
    shards::CellTopology hexTopo = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<> >() );
    runSubBasisTests(hexTopo, out, success);
  }
  
} // namespace
