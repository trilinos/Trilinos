// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    using DeviceType = DefaultTestDeviceType;
    using OutputScalar = double;
    using PointScalar  = double;
    using namespace Intrepid2;
    using Basis = Basis<DeviceType,OutputScalar,PointScalar>;
    
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
    
    bool isLine = cellTopo.getKey() == shards::Line<>::key;
    bool isQuad = cellTopo.getKey() == shards::Quadrilateral<>::key;
    bool isHex  = cellTopo.getKey() == shards::Hexahedron<>::key;
    bool isWedge = cellTopo.getKey() == shards::Wedge<>::key;
    int polyOrderDim = -1; // the number of dimensions of p-anisotropy allowed
    if (isLine || isQuad || isHex)
    {
      polyOrderDim = spaceDim;
    }
    else if (isWedge)
    {
      polyOrderDim = 2; // line x tri
    }
    else
    {
      polyOrderDim = 1;
    }
    auto subBasisDegreeTestCases = getBasisTestCasesUpToDegree(polyOrderDim, minDegree, polyOrder_x, polyOrder_y, polyOrder_z);
    
    std::vector<int> degrees(polyOrderDim);
    degrees[0] = polyOrder_x;
    if (polyOrderDim > 1) degrees[1] = polyOrder_y;
    if (polyOrderDim > 2) degrees[2] = polyOrder_z;
    
    int numPoints_1D = 5;
    auto inputPoints = getInputPointsView<PointScalar,DeviceType>(cellTopo, numPoints_1D);
    int numPoints = inputPoints.extent_int(0);
    
    auto op = Intrepid2::OPERATOR_VALUE;
    auto outputValues = getOutputView<OutputScalar,DeviceType>(fs, op, basis->getCardinality(), numPoints, spaceDim);
    
    basis->getValues(outputValues, inputPoints, op);
    out << "Testing sub-basis inclusion in degree ";
    for (int i=0; i<polyOrderDim; i++)
    {
      out << degrees[i];
      if (i<polyOrderDim-1) out << " x ";
    }
    out << " basis " << basis->getName() << std::endl;
    
    for (auto testCase : subBasisDegreeTestCases)
    {
      out << "testing sub-basis of degree ";
      for (int i=0; i<polyOrderDim; i++)
      {
        out << testCase[i];
        if (i<polyOrderDim-1) out << " x ";
      }
      out << std::endl;
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
      
      auto subBasisOutputValues = getOutputView<OutputScalar,DeviceType>(fs, op, subBasis->getCardinality(), numPoints, spaceDim);
      subBasis->getValues(subBasisOutputValues, inputPoints, op);
      Kokkos::fence();
      bool vectorValued = (outputValues.rank() == 3); // F,P,D -- if scalar-valued, F,P
      
      auto inputPointsHost          = getHostCopy(inputPoints);
      auto outputValuesHost         = getHostCopy(outputValues);
      auto subBasisOutputValuesHost = getHostCopy(subBasisOutputValues);
      
      for (int pointOrdinal=0; pointOrdinal<numPoints; pointOrdinal++)
      {
        // by construction, the sub-basis should have fields in the same order as the original basis
        // (the original basis just has extra fields interspersed)
        int subBasisFieldOrdinal = 0;
        for (int fieldOrdinal : lowerOrderFieldOrdinals)
        {
          if (!vectorValued)
          {
            double originalValue = outputValuesHost(fieldOrdinal,pointOrdinal);
            double subBasisValue = subBasisOutputValuesHost(subBasisFieldOrdinal,pointOrdinal);
            
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
              double x = inputPointsHost(pointOrdinal,0);
              double y = (spaceDim > 1) ? inputPointsHost(pointOrdinal,1) : -2.0;
              double z = (spaceDim > 2) ? inputPointsHost(pointOrdinal,2) : -2.0;
              
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
              double originalValue = outputValuesHost(fieldOrdinal,pointOrdinal,d);
              double subBasisValue = subBasisOutputValuesHost(subBasisFieldOrdinal,pointOrdinal,d);
              
              if (!essentiallyEqual(originalValue, subBasisValue, tol))
              {
                valuesMatch = false;
              }
            }
            
            if (!valuesMatch)
            {
              double x = inputPointsHost(pointOrdinal,0);
              double y = (spaceDim > 1) ? inputPointsHost(pointOrdinal,1) : -2.0;
              double z = (spaceDim > 2) ? inputPointsHost(pointOrdinal,2) : -2.0;
              
              if (spaceDim == 1)
                out << "values for "  << x  << " differ for field ordinal " << fieldOrdinal;
              else if (spaceDim == 2)
                out << "values for ("  << x << "," << y << ") differ for field ordinal " << fieldOrdinal;
              else
                out << "values for ("  << x << "," << y << "," << z << ") differ for field ordinal " << fieldOrdinal;
              out << ": expected value (lower-order basis fieldOrdinal " << subBasisFieldOrdinal << "): (";
              for (int d=0; d<spaceDim; d++)
              {
                out << subBasisOutputValuesHost(subBasisFieldOrdinal,pointOrdinal,d);
                if (d<spaceDim-1) out << ",";
              }
              out << "); actual (larger basis fieldOrdinal " << fieldOrdinal << ") was (";
              for (int d=0; d<spaceDim; d++)
              {
                out << outputValuesHost(fieldOrdinal,pointOrdinal,d);
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
        for (int degree=2; degree<=maxDegree; degree++)
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

  TEUCHOS_UNIT_TEST( SubBasisInclusion, Pyramid )
  {
    shards::CellTopology pyrTopo = shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<> >() );
//    runSubBasisTests(pyrTopo, out, success);
    
    // so far, only HGRAD, HVOL implemented for hierarchical pyramid.  Once we have full exact sequence, we can
    // switch to calling runSubBasisTests(tetTopo, out, success).
    // NOTE: due to the way that lower-degree polynomials get added in the higher-degree bases as p increases (combined with the way that the bases are ordered), the strategy within runSubBasisTests for mapping basis ordinals from the "sub-basis" to the full basis WILL NOT WORK for H(curl) and H(div) pyramids.  We could specify the mapping in some fairly manual way (with awareness of implementation details).  Better still would be to provide a mapping within Basis that would specify how one basis is included in the other.  For now, we simply omit these tests.
    std::vector<Intrepid2::EFunctionSpace> functionSpaces = {FUNCTION_SPACE_HGRAD, FUNCTION_SPACE_HVOL};
    auto cellTopo = pyrTopo;
    
    const int maxDegree = 5;
    const double tol = TEST_TOLERANCE_TIGHT;
    
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
      for (auto continuousBasis : continuousBasisValues) // corresponds to "defineVertexFunctions" in basis definitions
      {
        for (int degree=1; degree<=maxDegree; degree++)
        {
          testSubBasis(cellTopo, fs, tol, out, success, degree, -1, -1, continuousBasis);
        }
      }
    }
  }
  
  TEUCHOS_UNIT_TEST( SubBasisInclusion, Tetrahedron )
  {
    shards::CellTopology tetTopo = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<> >() );
    runSubBasisTests(tetTopo, out, success);
  }
  
  TEUCHOS_UNIT_TEST( SubBasisInclusion, Triangle )
  {
    shards::CellTopology triTopo = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<> >() );
    
    // so far, only HGRAD implemented for hierarchical triangle.  Once we have full exact sequence, we can
    // switch to calling runSubBasisTests(triTopo, out, success).
    std::vector<Intrepid2::EFunctionSpace> functionSpaces = {FUNCTION_SPACE_HGRAD};
    auto cellTopo = triTopo;
    
    const int maxDegree = 5;
    const double tol = TEST_TOLERANCE_TIGHT;
    
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
          testSubBasis(cellTopo, fs, tol, out, success, degree, -1, -1,continuousBasis);
        }
      }
    }
  }
  
} // namespace
