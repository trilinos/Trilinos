// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   AnalyticPolynomialsMatchTests.cpp
    \brief  Tests to verify that basis implementations match the analytic polynomials on which they are based.
    \author Created by N.V. Roberts.
 */
 
#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_Sacado.hpp" // important to include this prior to other Intrepid2 headers.

#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"
#include "Intrepid2_Polynomials.hpp"
#include "Intrepid2_Types.hpp"

#include "Intrepid2_TestUtils.hpp"

#include "Kokkos_Core.hpp"

namespace
{
  using namespace Intrepid2;
  
  template<typename Scalar>
  Scalar integratedJacobi(Scalar x, Scalar t, double alpha, const int n, const int derivativeOrder = 0)
  {
    // ideally, we would have closed-form expressions somewhere in here, as we do for integrated Legendre below
    // for now, though, this is just a thin wrapper around the call to Intrepid2::Polynomials::shiftedScaledIntegratedJacobiValues()
    auto values = getView<Scalar,Kokkos::HostSpace>("integrated Jacobi values", n+1);
    Polynomials::shiftedScaledIntegratedJacobiValues(values, alpha, n, x, t);
    return values(n);
  }
  
  template<typename Scalar>
  Scalar integratedLegendreAnalytic(Scalar x, const int n, bool useMinusOneToOne, const int derivativeOrder = 0)
  {
    // formulas below are for x in [-1,1]; if we are using [0,1], need to remap appropriately:
    const double derivativeScaling = useMinusOneToOne ? 1.0 : pow(2.0, derivativeOrder);
    if (!useMinusOneToOne)
    {
      x = 2.0 * x - 1.0;
    }
    Scalar value;
    switch (derivativeOrder)
    {
      case 0:
        switch (n)
      {
        case 0:
          value = (1.0-x)/2.0;                   // left vertex function (node at -1)
          break;
        case 1:
          value = (1.0+x)/2.0;                   // right vertex function (node at 1)
          break;
        case 2:
          value = (x*x-1.0)/4.0;                 // L_2 : (x^2 - 1) / 4
          break;
        case 3:
          value = (x*x*x-x)/4.0;                 // L_3 : (x^3 - x) / 4
          break;
        case 4:
          value = (5.*x*x*x*x - 6.*x*x+1.)/16.0; // L_4 : (5x^4-6x^2+1) / 16
          break;
        default:
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "unsupported n");
      }
        break;
      case 1:
        switch (n)
      {
        case 0:
          value = -1.0/2.0;              // left vertex function (node at -1)
          break;
        case 1:
          value = 1.0/2.0;               // right vertex function (node at 1)
          break;
        case 2:
          value = x/2.0;                 // L_2 : (x^2 - 1) / 4
          break;
        case 3:
          value = (3.0*x*x-1.0)/4.0;     // L_3 : (x^3 - x) / 4
          break;
        case 4:
          value = (5.*x*x*x - 3.*x)/4.0; // L_4 : (5x^4-6x^2+1) / 16
          break;
        default:
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "unsupported n");
      }
        break;
      case 2:
        switch (n)
      {
        case 0:
          value = 0.0;               // left vertex function (node at -1)
          break;
        case 1:
          value = 0.0;               // right vertex function (node at 1)
          break;
        case 2:
          value = 1.0/2.0;           // L_2 : (x^2 - 1) / 4
          break;
        case 3:
          value = 3.0*x/2.0;         // L_3 : (x^3 - x) / 4
          break;
        case 4:
          value = (15.*x*x - 3.)/4.; // L_4 : (5x^4-6x^2+1) / 16
          break;
        default:
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "unsupported n");
      }
        break;
      case 3:
        switch (n)
      {
        case 0:
          value = 0.0;               // left vertex function (node at -1)
          break;
        case 1:
          value = 0.0;               // right vertex function (node at 1)
          break;
        case 2:
          value = 0.0;               // L_2 : (x^2 - 1) / 4
          break;
        case 3:
          value = 3.0/2.0;           // L_3 : (x^3 - x) / 4
          break;
        case 4:
          value = (15.*x)/2.;        // L_4 : (5x^4-6x^2+1) / 16
          break;
        default:
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "unsupported n");
      }
        break;
      case 4:
        switch (n)
      {
        case 0:
          value = 0.0;           // left vertex function (node at -1)
          break;
        case 1:
          value = 0.0;           // right vertex function (node at 1)
          break;
        case 2:
          value = 0.0;           // L_2 : (x^2 - 1) / 4
          break;
        case 3:
          value = 0.0;           // L_3 : (x^3 - x) / 4
          break;
        case 4:
          value = 15./2.;        // L_4 : (5x^4-6x^2+1) / 16
          break;
        default:
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "unsupported n");
      }
        break;
      default:
        value = 0.0;
    }
    return value * derivativeScaling;
  }
  
  template<typename Scalar>
  Scalar shiftedScaledIntegratedLegendreAnalytic(Scalar x, Scalar t, const int n)
  {
    bool useMinusOneToOne = false; // our domain is [0,1]
    const int derivativeOrder = 0;
    double tol = 1e-15;
    using std::abs;
    if ((n == 0) || (n == 1))
    {
      // linear in x: scaling by t cancels
      return integratedLegendreAnalytic(x, n, useMinusOneToOne, derivativeOrder);
    }
    else if (abs(t) < tol)
    {
      // define a scaling by 0 to be 0 (does this match what Fuentes et al do??)
      // even if this isn't perfectly general, it works for our tests (at least for OPERATOR_VALUE):
      // we only ever will end up with t=0 for edge functions on the triangle, evaluated at a vertex.  These vanish.
      return 0.0;
    }
    else
    {
      Scalar tPower = 1.0;
      for (int i=0; i<n; i++)
      {
        tPower *= t;
      }
      return tPower * integratedLegendreAnalytic(Scalar(x/t), n, useMinusOneToOne);
    }
  }
  
  template<typename DeviceType,
           typename OutputScalar = double,
           typename PointScalar  = double>
  void testHierarchicalHGRAD_LINE_MatchesAnalyticValues(Intrepid2::EOperator op, const double tol, Teuchos::FancyOStream &out, bool &success)
  {
    using namespace Intrepid2;
    using BasisFamily = HierarchicalBasisFamily<DeviceType,OutputScalar,PointScalar>;
    
    const int polyOrder = 4;
    auto hgradBasis = getLineBasis<BasisFamily>(FUNCTION_SPACE_HGRAD, polyOrder);
    
    int numPoints_1D = 5;
    shards::CellTopology lineTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >() );
    auto inputPointsView = getInputPointsView<PointScalar,DeviceType>(lineTopo, numPoints_1D);
    
    auto hgradOutputView = getOutputView<OutputScalar,DeviceType>(FUNCTION_SPACE_HGRAD, op, hgradBasis->getCardinality(), numPoints_1D, 1);
    
    hgradBasis->getValues(hgradOutputView, inputPointsView, op);
    
    auto hgradOutputViewHost = getHostCopy(hgradOutputView);
    auto inputPointsViewHost = getHostCopy(inputPointsView);
    
    auto expectedValuesView     = getView<OutputScalar,DeviceType>("expected values", hgradBasis->getCardinality());
    auto expectedValuesViewHost = getHostCopy(expectedValuesView);
    
    for (int pointOrdinal=0; pointOrdinal<numPoints_1D; pointOrdinal++)
    {
      int pointPassed = true;
      PointScalar x = inputPointsViewHost(pointOrdinal,0);
      
      const bool useMinusOneToOne = true; // Intrepid2's reference element
      int derivativeOrder;
      
      switch (op)
      {
        case Intrepid2::OPERATOR_VALUE:
          derivativeOrder = 0;
          break;
        case Intrepid2::OPERATOR_GRAD:
          derivativeOrder = 1;
          break;
        case Intrepid2::OPERATOR_D1:
        case Intrepid2::OPERATOR_D2:
        case Intrepid2::OPERATOR_D3:
        case Intrepid2::OPERATOR_D4:
        case Intrepid2::OPERATOR_D5:
        case Intrepid2::OPERATOR_D6:
        case Intrepid2::OPERATOR_D7:
        case Intrepid2::OPERATOR_D8:
        case Intrepid2::OPERATOR_D9:
        case Intrepid2::OPERATOR_D10:
          derivativeOrder = op - OPERATOR_D1 + 1;
          break;
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported operator");
      }
      for (int n=0; n<polyOrder+1; n++)
      {
        expectedValuesViewHost(n) = integratedLegendreAnalytic(x, n, useMinusOneToOne, derivativeOrder);
      }
      
      for (int fieldOrdinal=0; fieldOrdinal<hgradBasis->getCardinality(); fieldOrdinal++)
      {
        OutputScalar actual    = hgradOutputViewHost.access(fieldOrdinal,pointOrdinal,0);
        OutputScalar expected  = expectedValuesViewHost(fieldOrdinal);
        
        bool valuesMatch = true;
        bool valuesAreBothSmall = valuesAreSmall(actual, expected, tol);
        if (!valuesAreBothSmall)
        {
          TEUCHOS_TEST_FLOATING_EQUALITY(actual, expected, tol, out, valuesMatch);
        }
        
        if (!valuesMatch)
        {
          pointPassed = false;
          PointScalar x = inputPointsViewHost(pointOrdinal,0);
          if (op == OPERATOR_VALUE) out << "values";
          else
          {
            int derivativeOrder = getOperatorOrder(op);
            if (derivativeOrder == 1)
            {
              out << "first ";
            }
            else if (derivativeOrder == 2)
            {
              out << "second ";
            }
            else if (derivativeOrder == 3)
            {
              out << "third ";
            }
            else
            {
              out << derivativeOrder << "th ";
            }
            out << "derivatives";
          }
          out << " for "  << x  << " differ for field ordinal " << fieldOrdinal;
          out << ": expected " << expected << "; actual " << actual;
          out << " (diff: " << expected-actual << ")" << std::endl;
          success = false;
        }
      }
      if (!pointPassed)
      {
        out << "point " << pointOrdinal << " failed.\n";
      }
    }
  }
  
  void testHierarchicalHVOL_LINE_MatchesAnalyticValues(Intrepid2::EOperator op, const double tol, Teuchos::FancyOStream &out, bool &success)
  {
    using namespace Intrepid2;
    using DeviceType   = DefaultTestDeviceType;
    using BasisFamily  = HierarchicalBasisFamily<DeviceType>;
    using PointScalar  = BasisFamily::PointValueType;
    using OutputScalar = BasisFamily::OutputValueType;
    
    const int polyOrder = 4;
    auto hvolBasis = getLineBasis<BasisFamily>(FUNCTION_SPACE_HVOL, polyOrder);
    
    int numPoints_1D = 5;
    shards::CellTopology lineTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >() );
    auto inputPointsView = getInputPointsView<PointScalar,DeviceType>(lineTopo, numPoints_1D);
    
    auto hvolOutputView = getOutputView<OutputScalar,DeviceType>(FUNCTION_SPACE_HVOL, op, hvolBasis->getCardinality(), numPoints_1D, 1);
    
    hvolBasis->getValues(hvolOutputView, inputPointsView, op);
    
    auto hvolOutputViewHost  = getHostCopy(hvolOutputView);
    auto inputPointsViewHost = getHostCopy(inputPointsView);
    
    for (int pointOrdinal=0; pointOrdinal<numPoints_1D; pointOrdinal++)
    {
      int pointPassed = true;
      double x = inputPointsViewHost(pointOrdinal,0);
      
      std::vector<double> expectedValues(hvolBasis->getCardinality());
      switch (op)
      {
        case Intrepid2::OPERATOR_VALUE:
          expectedValues[0] = 1;                              // P_0 : 1
          expectedValues[1] = x;                              // P_1 : x
          expectedValues[2] = (3.*x*x-1.0)/2.0;               // P_2 : (3x^2 - 1) / 2
          expectedValues[3] = (5.*x*x*x-3.*x)/2.0;            // P_3 : (5x^3 - 3x) / 2
          expectedValues[4] = (35.*x*x*x*x - 30.*x*x+3.)/8.0; // P_4 : (35x^4-30x^2+3) / 8
          break;
        case Intrepid2::OPERATOR_GRAD:
        case Intrepid2::OPERATOR_D1:
          // first derivatives of the above:
          expectedValues[0] = 0.0;                      // P_0 : 1
          expectedValues[1] = 1.0;                      // P_1 : x
          expectedValues[2] = 3.*x;                     // P_2 : (3x^2 - 1) / 2
          expectedValues[3] = (15.0*x*x-3.0)/2.0;       // P_3 : (5x^3 - 3x) / 2
          expectedValues[4] = (35.*x*x*x - 15.*x)/2.0;  // P_4 : (35x^4-30x^2+3) / 8
          break;
        case Intrepid2::OPERATOR_D2:
          // second derivatives:
          expectedValues[0] = 0.0;                 // P_0 : 1
          expectedValues[1] = 0.0;                 // P_1 : x
          expectedValues[2] = 3.0;                 // P_2 : (3x^2 - 1) / 2
          expectedValues[3] = 15.0*x;              // P_3 : (5x^3 - 3x) / 2
          expectedValues[4] = (105.*x*x - 15.)/2.; // P_4 : (35x^4-30x^2+3) / 8
          break;
        case Intrepid2::OPERATOR_D3:
          // third derivatives:
          expectedValues[0] = 0.0;           // P_0 : 1
          expectedValues[1] = 0.0;           // P_1 : x
          expectedValues[2] = 0.0;           // P_2 : (3x^2 - 1) / 2
          expectedValues[3] = 15.0;          // P_3 : (5x^3 - 3x) / 2
          expectedValues[4] = 105.*x;        // P_4 : (35x^4-30x^2+3) / 8
          break;
        case Intrepid2::OPERATOR_D4:
          // fourth derivatives:
          expectedValues[0] = 0.0;           // P_0 : 1
          expectedValues[1] = 0.0;           // P_1 : x
          expectedValues[2] = 0.0;           // P_2 : (3x^2 - 1) / 2
          expectedValues[3] = 0.0;           // P_3 : (5x^3 - 3x) / 2
          expectedValues[4] = 105.;          // P_4 : (35x^4-30x^2+3) / 8
          break;
        case Intrepid2::OPERATOR_D5:
        case Intrepid2::OPERATOR_D6:
        case Intrepid2::OPERATOR_D7:
        case Intrepid2::OPERATOR_D8:
        case Intrepid2::OPERATOR_D9:
        case Intrepid2::OPERATOR_D10:
          // nth (n≥5) derivatives are all 0:
          expectedValues[0] = 0.0;           // P_0 : 1
          expectedValues[1] = 0.0;           // P_1 : x
          expectedValues[2] = 0.0;           // P_2 : (3x^2 - 1) / 2
          expectedValues[3] = 0.0;           // P_3 : (5x^3 - 3x) / 2
          expectedValues[4] = 0.0;           // P_4 : (35x^4-30x^2+3) / 8
          break;
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported operator");
      }
      
      for (int fieldOrdinal=0; fieldOrdinal<hvolBasis->getCardinality(); fieldOrdinal++)
      {
        double actual    = hvolOutputViewHost.access(fieldOrdinal,pointOrdinal,0);
        double expected  = expectedValues[fieldOrdinal];
        
        bool valuesMatch = true;
        bool valuesAreBothSmall = valuesAreSmall(actual, expected, tol);
        if (!valuesAreBothSmall)
        {
          TEUCHOS_TEST_FLOATING_EQUALITY(actual, expected, tol, out, valuesMatch);
        }
        
        if (!valuesMatch)
        {
          pointPassed = false;
          double x = inputPointsViewHost(pointOrdinal,0);
          if (op == OPERATOR_VALUE) out << "values";
          else
          {
            int derivativeOrder = getOperatorOrder(op);
            if (derivativeOrder == 1)
            {
              out << "first ";
            }
            else if (derivativeOrder == 2)
            {
              out << "second ";
            }
            else if (derivativeOrder == 3)
            {
              out << "third ";
            }
            else
            {
              out << derivativeOrder << "th ";
            }
            out << "derivatives";
          }
          out << " for "  << x  << " differ for field ordinal " << fieldOrdinal;
          out << ": expected " << expected << "; actual " << actual;
          out << " (diff: " << expected-actual << ")" << std::endl;
          success = false;
        }
      }
      if (!pointPassed)
      {
        out << "point " << pointOrdinal << " failed.\n";
      }
    }
  }
  
  template<typename DeviceType,
           typename OutputScalar = double,
           typename PointScalar  = double>
  void testHierarchicalHGRAD_TRIANGLE_MatchesAnalyticValues(Intrepid2::EOperator op, const double tol, Teuchos::FancyOStream &out, bool &success)
  {
    using namespace Intrepid2;
    using BasisFamily = HierarchicalBasisFamily<DeviceType,OutputScalar,PointScalar>;
    
    const int spaceDim  = 2;
    const int polyOrder = 4;
    auto hgradBasis = getTriangleBasis<BasisFamily>(FUNCTION_SPACE_HGRAD, polyOrder);
    
    const int numVertices                  = 3;
    const int numFunctionsPerVertex        = 1;
    const int numVertexFunctionsExpected   = numVertices * numFunctionsPerVertex;
    const int numEdges                     = 3;
    const int num1DEdgeFunctions           = (polyOrder + 1) - 2; // line basis cardinality, minus two vertex functions
    const int numEdgeFunctionsExpected     = num1DEdgeFunctions * numEdges;
    const int numInteriorFunctionsExpected = (num1DEdgeFunctions-1)*num1DEdgeFunctions/2; // triangular sum
    const int expectedCardinality = numVertexFunctionsExpected + numEdgeFunctionsExpected + numInteriorFunctionsExpected;
    
    TEST_EQUALITY(expectedCardinality, hgradBasis->getCardinality());
    // could also test vertex/edge/face function count individually (worth doing, I think)
    
    int numPoints_1D = 5;
    shards::CellTopology triangleTopo = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<> >() );
    auto inputPointsView = getInputPointsView<PointScalar,DeviceType>(triangleTopo, numPoints_1D);
    
    const int numPoints = inputPointsView.extent_int(0);
        
    auto hgradOutputView = getOutputView<OutputScalar,DeviceType>(FUNCTION_SPACE_HGRAD, op, hgradBasis->getCardinality(), numPoints, spaceDim);
    hgradBasis->getValues(hgradOutputView, inputPointsView, op);
    
    auto hgradOutputViewHost  = getHostCopy(hgradOutputView);
    auto inputPointsViewHost = getHostCopy(inputPointsView);
    
    auto expectedValuesView = getOutputView<OutputScalar,DeviceType>(FUNCTION_SPACE_HGRAD, op, hgradBasis->getCardinality(), numPoints, spaceDim);
    auto expectedValuesViewHost = getHostCopy(expectedValuesView);
    
    auto legendreValuesAtPoint = getView<OutputScalar,DeviceType>("Legendre values temporary storage", polyOrder+1);
    auto legendreValuesAtPointHost = getHostCopy(legendreValuesAtPoint);
    
    for (int pointOrdinal=0; pointOrdinal<numPoints; pointOrdinal++)
    {
      int pointPassed = true;
      PointScalar x = inputPointsViewHost(pointOrdinal,0);
      PointScalar y = inputPointsViewHost(pointOrdinal,1);
      
      out << "Checking point (" << x << "," << y << ")\n";
      
      // write as barycentric coordinates:
      const PointScalar lambda[3] = {1. - x - y, x, y};
      const int edge_start[3] = {0,1,0};
      const int edge_end[3]   = {1,2,2};
      
      switch (op)
      {
        case Intrepid2::OPERATOR_VALUE:
        {
          // vertex polynomials come first, according to vertex ordering: (0,0), (1,0), (0,1)
          for (int vertexOrdinal=0; vertexOrdinal<numVertexFunctionsExpected; vertexOrdinal++)
          {
            expectedValuesViewHost(vertexOrdinal,pointOrdinal) = lambda[vertexOrdinal];
          }
          
          int fieldOrdinalOffset = 3;
          for (int edgeOrdinal=0; edgeOrdinal<numEdges; edgeOrdinal++)
          {
            const auto & s0 = lambda[edge_start[edgeOrdinal]];
            const auto & s1 = lambda[edge_end[edgeOrdinal]];
            const PointScalar t = s0 + s1;
            
            for (int edgeFunctionOrdinal=0; edgeFunctionOrdinal<num1DEdgeFunctions; edgeFunctionOrdinal++)
            {
              expectedValuesViewHost(edgeFunctionOrdinal+fieldOrdinalOffset,pointOrdinal) = shiftedScaledIntegratedLegendreAnalytic(s1, t, edgeFunctionOrdinal+2);
            }
            fieldOrdinalOffset += num1DEdgeFunctions;
          }
          // face functions
          for (int i=2; i<polyOrder; i++)
          {
            // we use edge function values from the 01 edge and blend with integrated Jacobi
            // the 01 edge function values start just after the vertex functions; since i starts at 2, our offset is therefore 1:
            const int offset = numVertexFunctionsExpected - 2;
            const OutputScalar & edgeFunctionValue = expectedValuesViewHost(i+offset,pointOrdinal);
            double alpha = i * 2.0;
            for (int j=1; i+j<=polyOrder; j++)
            {
              const PointScalar  &x = lambda[2];
              const PointScalar   t = 1.0;
              const OutputScalar jacobiValue = integratedJacobi(x, t, alpha, j);
              expectedValuesViewHost(fieldOrdinalOffset,pointOrdinal) = edgeFunctionValue * jacobiValue;
              fieldOrdinalOffset++;
            }
          }
        }
          break;
//        case Intrepid2::OPERATOR_GRAD:
//        case Intrepid2::OPERATOR_D1:
//          // first derivatives of the above:
//          expectedValues[0] = 0.0;                      // P_0 : 1
//          expectedValues[1] = 1.0;                      // P_1 : x
//          expectedValues[2] = 3.*x;                     // P_2 : (3x^2 - 1) / 2
//          expectedValues[3] = (15.0*x*x-3.0)/2.0;       // P_3 : (5x^3 - 3x) / 2
//          expectedValues[4] = (35.*x*x*x - 15.*x)/2.0;  // P_4 : (35x^4-30x^2+3) / 8
//          break;
//        case Intrepid2::OPERATOR_D2:
//          // second derivatives:
//          expectedValues[0] = 0.0;                 // P_0 : 1
//          expectedValues[1] = 0.0;                 // P_1 : x
//          expectedValues[2] = 3.0;                 // P_2 : (3x^2 - 1) / 2
//          expectedValues[3] = 15.0*x;              // P_3 : (5x^3 - 3x) / 2
//          expectedValues[4] = (105.*x*x - 15.)/2.; // P_4 : (35x^4-30x^2+3) / 8
//          break;
//        case Intrepid2::OPERATOR_D3:
//          // third derivatives:
//          expectedValues[0] = 0.0;           // P_0 : 1
//          expectedValues[1] = 0.0;           // P_1 : x
//          expectedValues[2] = 0.0;           // P_2 : (3x^2 - 1) / 2
//          expectedValues[3] = 15.0;          // P_3 : (5x^3 - 3x) / 2
//          expectedValues[4] = 105.*x;        // P_4 : (35x^4-30x^2+3) / 8
//          break;
//        case Intrepid2::OPERATOR_D4:
//          // fourth derivatives:
//          expectedValues[0] = 0.0;           // P_0 : 1
//          expectedValues[1] = 0.0;           // P_1 : x
//          expectedValues[2] = 0.0;           // P_2 : (3x^2 - 1) / 2
//          expectedValues[3] = 0.0;           // P_3 : (5x^3 - 3x) / 2
//          expectedValues[4] = 105.;          // P_4 : (35x^4-30x^2+3) / 8
//          break;
//        case Intrepid2::OPERATOR_D5:
//        case Intrepid2::OPERATOR_D6:
//        case Intrepid2::OPERATOR_D7:
//        case Intrepid2::OPERATOR_D8:
//        case Intrepid2::OPERATOR_D9:
//        case Intrepid2::OPERATOR_D10:
//          // nth (n≥5) derivatives are all 0:
//          expectedValues[0] = 0.0;           // P_0 : 1
//          expectedValues[1] = 0.0;           // P_1 : x
//          expectedValues[2] = 0.0;           // P_2 : (3x^2 - 1) / 2
//          expectedValues[3] = 0.0;           // P_3 : (5x^3 - 3x) / 2
//          expectedValues[4] = 0.0;           // P_4 : (35x^4-30x^2+3) / 8
//          break;
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported operator");
      }
      
      for (int fieldOrdinal=0; fieldOrdinal<hgradBasis->getCardinality(); fieldOrdinal++)
      {
        for (int d=0; d<hgradOutputViewHost.extent_int(2); d++)
        {
          OutputScalar actual    = hgradOutputViewHost.access(fieldOrdinal,pointOrdinal,d);
          OutputScalar expected  = expectedValuesViewHost.access(fieldOrdinal,pointOrdinal,d);
        
          bool valuesMatch = true;
          bool valuesAreBothSmall = valuesAreSmall(actual, expected, tol);
          if (!valuesAreBothSmall)
          {
            TEUCHOS_TEST_FLOATING_EQUALITY(actual, expected, tol, out, valuesMatch);
          }
        
          if (!valuesMatch)
          {
            pointPassed = false;
            PointScalar x = inputPointsViewHost(pointOrdinal,0);
            PointScalar y = inputPointsViewHost(pointOrdinal,1);
            if (op == OPERATOR_VALUE) out << "values";
            else
            {
              int derivativeOrder = getOperatorOrder(op);
              if (derivativeOrder == 1)
              {
                out << "first ";
              }
              else if (derivativeOrder == 2)
              {
                out << "second ";
              }
              else if (derivativeOrder == 3)
              {
                out << "third ";
              }
              else
              {
                out << derivativeOrder << "th ";
              }
              out << "derivatives";
            }
            out << " for ("  << x << "," << y << ") differ for field ordinal " << fieldOrdinal;
            out << ": expected " << expected << "; actual " << actual;
            out << " (diff: " << expected-actual << ")" << std::endl;
            success = false;
          }
        }
      }
      if (!pointPassed)
      {
        out << "point " << pointOrdinal << " failed.\n";
      }
    }
  }
  
  // compare derivatives of order derivativeOrder in H(grad) with derivatives of order (derivativeOrder-1) in H(vol)
  template<class LineBasisFamily>
  void testDerivativesMatch(int polyOrder, int derivativeOrder, const double tol, Teuchos::FancyOStream &out, bool &success)
  {
    using namespace Intrepid2;
    using Teuchos::rcp;
    
    using DeviceType   = typename LineBasisFamily::DeviceType;
    using PointScalar  = typename LineBasisFamily::PointValueType;
    using OutputScalar = typename LineBasisFamily::OutputValueType;
    
    auto hgradBasis = getLineBasis<LineBasisFamily>(FUNCTION_SPACE_HGRAD, polyOrder);
    auto hvolBasis  = getLineBasis<LineBasisFamily>(FUNCTION_SPACE_HVOL, polyOrder);
    
    auto opGrad = Intrepid2::EOperator(OPERATOR_D1 - 1 + derivativeOrder);
    auto opVol  = (derivativeOrder > 1) ? Intrepid2::EOperator(OPERATOR_D1 - 2 + derivativeOrder) : OPERATOR_VALUE;
    
    int numPoints_1D = 5;
    shards::CellTopology lineTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >() );
    auto inputPointsView = getInputPointsView<PointScalar,DeviceType>(lineTopo, numPoints_1D);
    const int spaceDim = 1;
    
    auto hgradOutputView = getOutputView<OutputScalar,DeviceType>(FUNCTION_SPACE_HGRAD, opGrad, hgradBasis->getCardinality(), numPoints_1D, spaceDim);
    auto hvolOutputView  = getOutputView<OutputScalar,DeviceType>(FUNCTION_SPACE_HVOL, opVol, hvolBasis->getCardinality(), numPoints_1D, spaceDim);
    
    hgradBasis->getValues(hgradOutputView, inputPointsView, opGrad);
    hvolBasis->getValues(  hvolOutputView, inputPointsView, opVol);
    
    auto hgradOutputViewHost = getHostCopy(hgradOutputView);
    auto hvolOutputViewHost  = getHostCopy(hvolOutputView);
    auto inputPointsViewHost = getHostCopy(inputPointsView);
    
    for (int pointOrdinal=0; pointOrdinal<numPoints_1D; pointOrdinal++)
    {
      int pointPassed = true;
      for (int fieldOrdinal=2; fieldOrdinal<hgradBasis->getCardinality(); fieldOrdinal++)
      {
        // relationship between our H^1 basis and our L^2 (see the "Analytic" tests above for explicit polynomial expressions)
        // - nth derivative of H^1 fieldOrdinal i is one half of the (n-1)th derivative of L^2 fieldOrdinal i-1 (for i>=2)
        OutputScalar hgradValue = hgradOutputViewHost(fieldOrdinal,  pointOrdinal,0);
        OutputScalar hvolValue  =  hvolOutputViewHost(fieldOrdinal-1,pointOrdinal,0);
        
        OutputScalar actual = hgradValue;
        OutputScalar expected = OutputScalar(0.5 * hvolValue);
        
        bool valuesMatch = true;
        
        bool valuesAreBothSmall = valuesAreSmall(actual, expected, tol);
        if (!valuesAreBothSmall)
        {
          TEUCHOS_TEST_FLOATING_EQUALITY(actual, expected, tol, out, valuesMatch);
        }
        
        if (!valuesMatch)
        {
          pointPassed = false;
          PointScalar x = inputPointsViewHost(pointOrdinal,0);
          
          std::string gradDerivative, volDerivative;
          if (derivativeOrder == 1)
          {
            gradDerivative = "st";
            volDerivative  = "th";
          }
          else if (derivativeOrder == 2)
          {
            gradDerivative = "nd";
            volDerivative  = "st";
          }
          else if (derivativeOrder == 3)
          {
            gradDerivative = "rd";
            volDerivative  = "nd";
          }
          else if (derivativeOrder == 4)
          {
            gradDerivative = "th";
            volDerivative  = "rd";
          }
          else
          {
            gradDerivative = "th";
            volDerivative  = "th";
          }
          
          out << "values for "  << x  << " differ for field ordinal " << fieldOrdinal;
          out << ": h(grad) " << derivativeOrder << gradDerivative << " derivative is " << hgradValue << "; h(vol) " << derivativeOrder-1 << volDerivative << " derivative is " << hvolValue;
          out << " (diff: " << hgradValue-hvolValue << ")" << std::endl;
          success = false;
        }
        
      }
      if (!pointPassed)
      {
        out << "point " << pointOrdinal << " failed.\n";
      }
    }
  }

  // compare basis values between standard Intrepid2 nodal basis and our derived nodal basis family members
  template<class DerivedNodalBasisFamily, class StandardNodalBasisFamily>
  void testNodalBasisMatches(shards::CellTopology &cellTopo, Intrepid2::EFunctionSpace fs, Intrepid2::EOperator op, int polyOrder, const double tol,
                             Teuchos::FancyOStream &out, bool &success)
  {
    using namespace Intrepid2;
    using BasisPtr       = typename DerivedNodalBasisFamily::BasisPtr;
    using ExecutionSpace = typename DerivedNodalBasisFamily::ExecutionSpace;
    using DeviceType     = typename DerivedNodalBasisFamily::DeviceType;
    using Teuchos::rcp;
    
    BasisPtr derivedBasis, standardBasis;
    
    if (cellTopo.getKey() == shards::Line<>::key)
    {
      derivedBasis  = getLineBasis<DerivedNodalBasisFamily> (fs, polyOrder);
      standardBasis = getLineBasis<StandardNodalBasisFamily>(fs, polyOrder);
    }
    else if (cellTopo.getKey() == shards::Quadrilateral<>::key)
    {
      derivedBasis  = getQuadrilateralBasis<DerivedNodalBasisFamily> (fs, polyOrder); // derived basis supports both isotropic and anisotropic polyOrder
      standardBasis = getQuadrilateralBasis<StandardNodalBasisFamily>(fs, polyOrder); // standard basis is isotropic
    }
    else if (cellTopo.getKey() == shards::Hexahedron<>::key)
    {
      derivedBasis  = getHexahedronBasis<DerivedNodalBasisFamily> (fs, polyOrder); // derived basis supports both isotropic and anisotropic polyOrder
      standardBasis = getHexahedronBasis<StandardNodalBasisFamily>(fs, polyOrder); // isotropic
    }
    
    int standardCardinality = standardBasis->getCardinality();
    int derivedCardinality = derivedBasis->getCardinality();
    
    if (standardCardinality != derivedCardinality)
    {
      success = false;
      out << "FAILURE: standard basis cardinality (" << standardCardinality;
      out << ") does not match derived basis cardinality (" << derivedCardinality << ")\n";
    }
    
    int spaceDim = cellTopo.getDimension();
    // we do allow ordering to be different.  dofMapToDerived maps from standard field ordinal to the derived.
    // we use getDofCoords() to perform the mapping.
    using ScalarViewType         = typename DerivedNodalBasisFamily::Basis::ScalarViewType;
    using OrdinalTypeArray1D     = typename DerivedNodalBasisFamily::Basis::OrdinalTypeArray1D;
    OrdinalTypeArray1D     dofMapToStandard     = OrdinalTypeArray1D("dofMapToStandard",standardCardinality);
//    OrdinalTypeArray1DHost dofMapToStandardHost = Kokkos::create_mirror_view(dofMapToStandard);
    
    using ValueType    = typename ScalarViewType::value_type;
    using ResultLayout = typename DeduceLayout< ScalarViewType >::result_layout;
    using AllocatableScalarViewType = Kokkos::DynRankView<ValueType, ResultLayout, DeviceType >;
    
    AllocatableScalarViewType dofCoordsStandard ("dofCoordsStandard", standardCardinality, spaceDim);
    AllocatableScalarViewType dofCoordsDerived  ("dofCoordsDerived",  standardCardinality, spaceDim);
    standardBasis->getDofCoords(dofCoordsStandard);
    derivedBasis-> getDofCoords(dofCoordsDerived );
    
    Kokkos::deep_copy(dofMapToStandard, -1);
    
    Kokkos::parallel_for(standardCardinality, KOKKOS_LAMBDA (const int fieldOrdinalStandard)
    {
      // search for a fieldOrdinalDerived that matches
      bool matchFound = false;
      for (int fieldOrdinalDerived=0; fieldOrdinalDerived<standardCardinality; fieldOrdinalDerived++)
      {
        bool matches = true;
        for (int d=0; d<spaceDim; d++)
        {
          if (dofCoordsStandard(fieldOrdinalStandard,d) != dofCoordsDerived(fieldOrdinalDerived,d))
          {
            matches = false;
            break;
          }
        }
        
        if (matches)
        {
          dofMapToStandard(fieldOrdinalDerived) = fieldOrdinalStandard;
          matchFound = true;
          break;
        }
      }
      if (!matchFound)
      {
        // failure; abort
        device_assert(matchFound);
      }
    });
    // copy dofMapToDerived to host view
//    Kokkos::deep_copy(dofMapToStandardHost, dofMapToStandard);
    
    using PointScalar  = typename DerivedNodalBasisFamily::PointValueType;
    using OutputScalar = typename DerivedNodalBasisFamily::OutputValueType;
    
    int numPoints_1D = 5;
    auto inputPointsView = getInputPointsView<PointScalar,DeviceType>(cellTopo, numPoints_1D);
    int numPoints = inputPointsView.extent_int(0);
    auto standardOutputView = getOutputView<OutputScalar,DeviceType>(fs, op, standardCardinality, numPoints, spaceDim);
    auto derivedOutputView  = getOutputView<OutputScalar,DeviceType>(fs, op, standardCardinality, numPoints, spaceDim);
    
    standardBasis->getValues(standardOutputView, inputPointsView, op);
    derivedBasis->getValues(derivedOutputView, inputPointsView, op);
    
    // derived may be in a different order.  Remap so it's in the same order as standard basis
    auto derivedOutputViewRemapped = getOutputView<OutputScalar,DeviceType>(fs, op, standardCardinality, numPoints, spaceDim);
    
    using ViewIteratorScalar = Intrepid2::ViewIterator<decltype(derivedOutputView), OutputScalar>;
    const int entryCount = standardOutputView.size();

    Kokkos::RangePolicy < ExecutionSpace > policy(0,entryCount);
    Kokkos::parallel_for( policy,
    KOKKOS_LAMBDA (const int &enumerationIndex )
    {
      ViewIteratorScalar vi1(derivedOutputView);
      vi1.setEnumerationIndex(enumerationIndex);
      
      auto location = vi1.getLocation();
      
      const int derivedFieldOrdinal = location[0];
      const int standardFieldOrdinal = dofMapToStandard(derivedFieldOrdinal);
      
      location[0] = standardFieldOrdinal;
      
      ViewIteratorScalar vi2(derivedOutputViewRemapped);
      vi2.setLocation(location);
      vi2.set(vi1.get());
    }
    );
    
    testViewFloatingEquality(standardOutputView, derivedOutputViewRemapped, tol, tol, out, success);
  }
  
  template<class DerivedNodalBasisFamily, class StandardNodalBasisFamily>
  void runNodalBasisComparisonTests(const int polyOrderToTest, shards::CellTopology &cellTopo, std::vector<Intrepid2::EFunctionSpace> functionSpaces, std::vector<Intrepid2::EOperator> ops,
                                    const double tol, Teuchos::FancyOStream &out, bool &success)
  {
    using namespace Intrepid2;
    
    for (auto fs : functionSpaces)
    {
      for (auto op : ops)
      {
        testNodalBasisMatches<DerivedNodalBasisFamily,StandardNodalBasisFamily>(cellTopo, fs, op, polyOrderToTest, tol, out, success);
      }
    }
  }
  
//  TEUCHOS_UNIT_TEST( AnalyticPolynomialsMatch, Hierarchical_HGRAD_LINE )
//  {
//    const double tol = TEST_TOLERANCE_TIGHT;
//
//    std::vector<Intrepid2::EOperator> operators = {{OPERATOR_VALUE, OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5, OPERATOR_D6, OPERATOR_D7, OPERATOR_D8, OPERATOR_D9, OPERATOR_D10}};
//    for (auto op : operators)
//    {
//      testHierarchicalHGRAD_LINE_MatchesAnalyticValues(op, tol, out, success);
//    }
//  }
//
  TEUCHOS_UNIT_TEST( AnalyticPolynomialsMatch, Hierarchical_HVOL_LINE )
  {
    const double tol = TEST_TOLERANCE_TIGHT;

    std::vector<Intrepid2::EOperator> operators = {{OPERATOR_VALUE, OPERATOR_D1, OPERATOR_D2, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5, OPERATOR_D6, OPERATOR_D7, OPERATOR_D8, OPERATOR_D9, OPERATOR_D10}};
    for (auto op : operators)
    {
      testHierarchicalHVOL_LINE_MatchesAnalyticValues(op, tol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( AnalyticPolynomialsMatch, Hierarchical_HGRAD_LINE, OutputScalar, PointScalar )
  {
    const double tol = TEST_TOLERANCE_TIGHT;
    using DeviceType = DefaultTestDeviceType;
    
    std::vector<Intrepid2::EOperator> operators = {{OPERATOR_VALUE, OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5, OPERATOR_D6, OPERATOR_D7, OPERATOR_D8, OPERATOR_D9, OPERATOR_D10}};
    for (auto op : operators)
    {
      testHierarchicalHGRAD_LINE_MatchesAnalyticValues<DeviceType,OutputScalar,PointScalar>(op, tol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( AnalyticPolynomialsMatch, Hierarchical_HGRAD_TRI, OutputScalar, PointScalar )
  {
    const double tol = TEST_TOLERANCE_TIGHT;
    using DeviceType = DefaultTestDeviceType;
    
//    std::vector<Intrepid2::EOperator> operators = {{OPERATOR_VALUE, OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5, OPERATOR_D6, OPERATOR_D7, OPERATOR_D8, OPERATOR_D9, OPERATOR_D10}};
    std::vector<Intrepid2::EOperator> operators = {OPERATOR_VALUE};
    for (auto op : operators)
    {
      testHierarchicalHGRAD_TRIANGLE_MatchesAnalyticValues<DeviceType,OutputScalar,PointScalar>(op, tol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( AnalyticPolynomialsMatch, Hierarchical_LineBasisDerivativesAgree, OutputScalar, PointScalar )
  {
    const int maxPolyOrder = 10;
    const double tol = TEST_TOLERANCE_TIGHT;
    
    using DeviceType = DefaultTestDeviceType;
    using HierarchicalBasisFamily = HierarchicalBasisFamily<DeviceType,OutputScalar,PointScalar>;
    
    for (int derivativeOrder=1; derivativeOrder<=5; derivativeOrder++)
    {
      // compare derivatives of order derivativeOrder in H(grad) with derivatives of order (derivativeOrder-1) in H(vol)
      testDerivativesMatch<HierarchicalBasisFamily>(maxPolyOrder, derivativeOrder, tol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( AnalyticPolynomialsMatch, HierarchicalNodalComparisons, OutputScalar, PointScalar )
  {
    // the following should match in H(grad) and H(vol), but are not expected to match in H(curl) and H(div)
    // (the latter differ in that Standard uses another H(grad) instance to represent the L^2 component of H(curl) and H(div), while Derived uses H(vol))
    using DeviceType = DefaultTestDeviceType;
    using DerivedNodalBasisFamily  = Intrepid2::DerivedNodalBasisFamily<DeviceType,OutputScalar,PointScalar>;
    using StandardNodalBasisFamily = Intrepid2::NodalBasisFamily       <DeviceType,OutputScalar,PointScalar>;
    
    const double tol = TEST_TOLERANCE_TIGHT;
    
    shards::CellTopology lineTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >() );
    shards::CellTopology quadTopo = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<> >() );
    shards::CellTopology hexTopo = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<> >() );
    
    std::vector<EOperator> operators_dk = {OPERATOR_D1,OPERATOR_D2,OPERATOR_D3,OPERATOR_D4,OPERATOR_D5,OPERATOR_D6,OPERATOR_D7,OPERATOR_D8,OPERATOR_D9,OPERATOR_D10};
    // for Fad types under CUDA, we sometimes run into allocation errors with the highest-k dk operators in 3D.  Not sure why (are we actually running out of memory? It seems possible, but still surprises me a little), but for now we don't test beyond D8.
    std::vector<EOperator> operators_dk_3D = {OPERATOR_D1,OPERATOR_D2,OPERATOR_D3,OPERATOR_D4,OPERATOR_D5,OPERATOR_D6,OPERATOR_D7,OPERATOR_D8};
    
    // "harder" test is probably higher-order, but it's more expensive in higher dimensions.
    // the below are compromises according to spatial dimension.
    const int polyOrder_1D = 4;
    const int polyOrder_2D = 3;
    const int polyOrder_3D = 2;
    
    out << "Running 1D nodal basis comparison tests…\n";
    runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(polyOrder_1D, lineTopo, {FUNCTION_SPACE_HVOL}, {OPERATOR_VALUE}, tol, out, success);
    runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(polyOrder_1D, lineTopo, {FUNCTION_SPACE_HGRAD}, {OPERATOR_VALUE,OPERATOR_GRAD}, tol, out, success);
    runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(polyOrder_1D, lineTopo, {FUNCTION_SPACE_HGRAD}, operators_dk, tol, out, success);

    out << "Running 2D nodal quadrilateral basis comparison tests…\n";
    runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(polyOrder_2D, quadTopo, {FUNCTION_SPACE_HVOL}, {OPERATOR_VALUE}, tol, out, success);
    runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(polyOrder_2D, quadTopo, {FUNCTION_SPACE_HGRAD}, {OPERATOR_VALUE,OPERATOR_GRAD}, tol, out, success);
    runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(polyOrder_2D, quadTopo, {FUNCTION_SPACE_HGRAD}, operators_dk, tol, out, success);

    out << "Running 3D nodal hexahedron basis comparison tests…\n";
    runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(polyOrder_3D, hexTopo, {FUNCTION_SPACE_HVOL}, {OPERATOR_VALUE}, tol, out, success);
    runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(polyOrder_3D, hexTopo, {FUNCTION_SPACE_HGRAD}, {OPERATOR_VALUE,OPERATOR_GRAD}, tol, out, success);
    runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(polyOrder_3D, hexTopo, {FUNCTION_SPACE_HGRAD}, operators_dk_3D, tol, out, success);
  }
                                                        
  INTREPID2_OUTPUTSCALAR_POINTSCALAR_TEST_INSTANT( AnalyticPolynomialsMatch, Hierarchical_HGRAD_LINE )
  INTREPID2_OUTPUTSCALAR_POINTSCALAR_TEST_INSTANT( AnalyticPolynomialsMatch, Hierarchical_LineBasisDerivativesAgree )
  INTREPID2_OUTPUTSCALAR_POINTSCALAR_TEST_INSTANT( AnalyticPolynomialsMatch, Hierarchical_HGRAD_TRI )
  INTREPID2_OUTPUTSCALAR_POINTSCALAR_TEST_INSTANT( AnalyticPolynomialsMatch, HierarchicalNodalComparisons )
} // namespace
