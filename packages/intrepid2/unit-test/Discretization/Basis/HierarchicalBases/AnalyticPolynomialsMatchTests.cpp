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

/** \file   AnalyticPolynomialsMatchTests.cpp
    \brief  Tests to verify that basis implementations match the analytic polynomials on which they are based.
    \author Created by N.V. Roberts.
 */
 
#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_Sacado.hpp" // important to include this prior to other Intrepid2 headers.

#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"
#include "Intrepid2_Types.hpp"

#include "Intrepid2_TestUtils.hpp"

#include "Kokkos_Core.hpp"

namespace
{
  using namespace Intrepid2;
  
  template<typename ExecutionSpace=Kokkos::DefaultExecutionSpace,
           typename OutputScalar = double,
           typename PointScalar  = double>
  void testHierarchicalHGRAD_LINE_MatchesAnalyticValues(Intrepid2::EOperator op, const double tol, Teuchos::FancyOStream &out, bool &success)
  {
    using namespace Intrepid2;
    using BasisFamily = HierarchicalBasisFamily<ExecutionSpace,OutputScalar,PointScalar>;
    
    const int polyOrder = 4;
    auto hgradBasis = getLineBasis<BasisFamily>(FUNCTION_SPACE_HGRAD, polyOrder);
    
    int numPoints_1D = 5;
    shards::CellTopology lineTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >() );
    auto inputPointsView = getInputPointsView<PointScalar>(lineTopo, numPoints_1D);
    
    auto hgradOutputView = getOutputView<OutputScalar>(FUNCTION_SPACE_HGRAD, op, hgradBasis->getCardinality(), numPoints_1D, 1);
    
    hgradBasis->getValues(hgradOutputView, inputPointsView, op);
    
    auto hgradOutputViewHost = getHostCopy(hgradOutputView);
    auto inputPointsViewHost = getHostCopy(inputPointsView);
    
    auto expectedValuesView     = getView<OutputScalar>("expected values", hgradBasis->getCardinality());
    auto expectedValuesViewHost = getHostCopy(expectedValuesView);
    
    for (int pointOrdinal=0; pointOrdinal<numPoints_1D; pointOrdinal++)
    {
      int pointPassed = true;
      PointScalar x = inputPointsViewHost(pointOrdinal,0);
      
      switch (op)
      {
        case Intrepid2::OPERATOR_VALUE:
          expectedValuesViewHost(0) = (1.0-x)/2.0;                   // left vertex function (node at -1)
          expectedValuesViewHost(1) = (1.0+x)/2.0;                   // right vertex function (node at 1)
          expectedValuesViewHost(2) = (x*x-1.0)/4.0;                 // L_2 : (x^2 - 1) / 4
          expectedValuesViewHost(3) = (x*x*x-x)/4.0;                 // L_3 : (x^3 - x) / 4
          expectedValuesViewHost(4) = (5.*x*x*x*x - 6.*x*x+1.)/16.0; // L_4 : (5x^4-6x^2+1) / 16
          break;
        case Intrepid2::OPERATOR_GRAD:
        case Intrepid2::OPERATOR_D1:
          // first derivatives of the above:
          expectedValuesViewHost(0) = -1.0/2.0;              // left vertex function (node at -1)
          expectedValuesViewHost(1) = 1.0/2.0;               // right vertex function (node at 1)
          expectedValuesViewHost(2) = x/2.0;                 // L_2 : (x^2 - 1) / 4
          expectedValuesViewHost(3) = (3.0*x*x-1.0)/4.0;     // L_3 : (x^3 - x) / 4
          expectedValuesViewHost(4) = (5.*x*x*x - 3.*x)/4.0; // L_4 : (5x^4-6x^2+1) / 16
          break;
        case Intrepid2::OPERATOR_D2:
          // second derivatives:
          expectedValuesViewHost(0) = 0.0;               // left vertex function (node at -1)
          expectedValuesViewHost(1) = 0.0;               // right vertex function (node at 1)
          expectedValuesViewHost(2) = 1.0/2.0;           // L_2 : (x^2 - 1) / 4
          expectedValuesViewHost(3) = 3.0*x/2.0;         // L_3 : (x^3 - x) / 4
          expectedValuesViewHost(4) = (15.*x*x - 3.)/4.; // L_4 : (5x^4-6x^2+1) / 16
          break;
        case Intrepid2::OPERATOR_D3:
          // third derivatives:
          expectedValuesViewHost(0) = 0.0;               // left vertex function (node at -1)
          expectedValuesViewHost(1) = 0.0;               // right vertex function (node at 1)
          expectedValuesViewHost(2) = 0.0;               // L_2 : (x^2 - 1) / 4
          expectedValuesViewHost(3) = 3.0/2.0;           // L_3 : (x^3 - x) / 4
          expectedValuesViewHost(4) = (15.*x)/2.;        // L_4 : (5x^4-6x^2+1) / 16
          break;
        case Intrepid2::OPERATOR_D4:
          // fourth derivatives:
          expectedValuesViewHost(0) = 0.0;           // left vertex function (node at -1)
          expectedValuesViewHost(1) = 0.0;           // right vertex function (node at 1)
          expectedValuesViewHost(2) = 0.0;           // L_2 : (x^2 - 1) / 4
          expectedValuesViewHost(3) = 0.0;           // L_3 : (x^3 - x) / 4
          expectedValuesViewHost(4) = 15./2.;        // L_4 : (5x^4-6x^2+1) / 16
          break;
        case Intrepid2::OPERATOR_D5:
        case Intrepid2::OPERATOR_D6:
        case Intrepid2::OPERATOR_D7:
        case Intrepid2::OPERATOR_D8:
        case Intrepid2::OPERATOR_D9:
        case Intrepid2::OPERATOR_D10:
          // fourth derivatives:
          expectedValuesViewHost(0) = 0.0;           // left vertex function (node at -1)
          expectedValuesViewHost(1) = 0.0;           // right vertex function (node at 1)
          expectedValuesViewHost(2) = 0.0;           // L_2 : (x^2 - 1) / 4
          expectedValuesViewHost(3) = 0.0;           // L_3 : (x^3 - x) / 4
          expectedValuesViewHost(4) = 0.0;           // L_4 : (5x^4-6x^2+1) / 16
          break;
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported operator");
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
    using BasisFamily  = HierarchicalBasisFamily<>;
    using PointScalar  = BasisFamily::PointValueType;
    using OutputScalar = BasisFamily::OutputValueType;
    
    const int polyOrder = 4;
    auto hvolBasis = getLineBasis<BasisFamily>(FUNCTION_SPACE_HVOL, polyOrder);
    
    int numPoints_1D = 5;
    shards::CellTopology lineTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >() );
    auto inputPointsView = getInputPointsView<PointScalar>(lineTopo, numPoints_1D);
    
    auto hvolOutputView = getOutputView<OutputScalar>(FUNCTION_SPACE_HVOL, op, hvolBasis->getCardinality(), numPoints_1D, 1);
    
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
  
  // compare derivatives of order derivativeOrder in H(grad) with derivatives of order (derivativeOrder-1) in H(vol)
  template<class LineBasisFamily>
  void testDerivativesMatch(int polyOrder, int derivativeOrder, const double tol, Teuchos::FancyOStream &out, bool &success)
  {
    using namespace Intrepid2;
    using Teuchos::rcp;
    
    using PointScalar  = typename LineBasisFamily::PointValueType;
    using OutputScalar = typename LineBasisFamily::OutputValueType;
    
    auto hgradBasis = getLineBasis<LineBasisFamily>(FUNCTION_SPACE_HGRAD, polyOrder);
    auto hvolBasis  = getLineBasis<LineBasisFamily>(FUNCTION_SPACE_HVOL, polyOrder);
    
    auto opGrad = Intrepid2::EOperator(OPERATOR_D1 - 1 + derivativeOrder);
    auto opVol  = (derivativeOrder > 1) ? Intrepid2::EOperator(OPERATOR_D1 - 2 + derivativeOrder) : OPERATOR_VALUE;
    
    int numPoints_1D = 5;
    shards::CellTopology lineTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >() );
    auto inputPointsView = getInputPointsView<PointScalar>(lineTopo, numPoints_1D);
    const int spaceDim = 1;
    
    auto hgradOutputView = getOutputView<OutputScalar>(FUNCTION_SPACE_HGRAD, opGrad, hgradBasis->getCardinality(), numPoints_1D, spaceDim);
    auto hvolOutputView  = getOutputView<OutputScalar>(FUNCTION_SPACE_HVOL, opVol, hvolBasis->getCardinality(), numPoints_1D, spaceDim);
    
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
    using BasisPtr = typename DerivedNodalBasisFamily::BasisPtr;
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
      derivedBasis  = getHexahedralBasis<DerivedNodalBasisFamily> (fs, polyOrder); // derived basis supports both isotropic and anisotropic polyOrder
      standardBasis = getHexahedralBasis<StandardNodalBasisFamily>(fs, polyOrder); // isotropic
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
    using ExecutionSpace         = typename DerivedNodalBasisFamily::Basis::ExecutionSpace;
    using ScalarViewType         = typename DerivedNodalBasisFamily::Basis::ScalarViewType;
    using OrdinalTypeArray1D     = typename DerivedNodalBasisFamily::Basis::OrdinalTypeArray1D;
    using OrdinalTypeArray1DHost = typename OrdinalTypeArray1D::HostMirror;
    OrdinalTypeArray1D     dofMapToDerived     = OrdinalTypeArray1D("dofMapToDerived",standardCardinality);
    OrdinalTypeArray1DHost dofMapToDerivedHost = Kokkos::create_mirror_view(dofMapToDerived);
    
    using ValueType    = typename ScalarViewType::value_type;
    using ResultLayout = typename DeduceLayout< ScalarViewType >::result_layout;
    using DeviceType = typename ScalarViewType::device_type;
    using ViewType = Kokkos::DynRankView<ValueType, ResultLayout, DeviceType >;
    
    ViewType dofCoordsStandard("dofCoordsStandard", standardCardinality, spaceDim);
    ViewType dofCoordsDerived ("dofCoordsDerived",  standardCardinality, spaceDim);
    standardBasis->getDofCoords(dofCoordsStandard);
    derivedBasis-> getDofCoords(dofCoordsDerived );
    
    Kokkos::deep_copy(dofMapToDerived, -1);
    
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
          dofMapToDerived(fieldOrdinalStandard) = fieldOrdinalDerived;
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
    Kokkos::deep_copy(dofMapToDerivedHost, dofMapToDerived);
    
    using PointScalar  = typename DerivedNodalBasisFamily::PointValueType;
    using OutputScalar = typename DerivedNodalBasisFamily::OutputValueType;
    
    int numPoints_1D = 5;
    auto inputPointsView = getInputPointsView<PointScalar>(cellTopo, numPoints_1D);
    int numPoints = inputPointsView.extent_int(0);
    auto standardOutputView = getOutputView<OutputScalar>(fs, op, standardCardinality, numPoints, spaceDim);
    auto derivedOutputView  = getOutputView<OutputScalar>(fs, op, standardCardinality, numPoints, spaceDim);
    
    standardBasis->getValues(standardOutputView, inputPointsView, op);
    derivedBasis->getValues(derivedOutputView, inputPointsView, op);
    
    auto standardOutputViewHost = getHostCopy(standardOutputView);
    auto derivedOutputViewHost  = getHostCopy(derivedOutputView);
    auto inputPointsViewHost    = getHostCopy(inputPointsView);
    
    bool scalarValued = (standardOutputView.rank() == 2); // F,P -- if vector-valued, F,P,D or F,P,Dk
    
    for (int pointOrdinal=0; pointOrdinal<numPoints; pointOrdinal++)
    {
      int pointPassed = true;
      for (int fieldOrdinalStandard=0; fieldOrdinalStandard<derivedCardinality; fieldOrdinalStandard++)
      {
        int fieldOrdinalDerived = dofMapToDerivedHost(fieldOrdinalStandard);
        if (scalarValued)
        {
          OutputScalar standardValue = standardOutputViewHost(fieldOrdinalStandard,pointOrdinal);
          OutputScalar derivedValue  =  derivedOutputViewHost(fieldOrdinalDerived, pointOrdinal);
          
          bool valuesMatch = true;
          TEUCHOS_TEST_FLOATING_EQUALITY(standardValue, derivedValue, tol, out, valuesMatch);
          
          if (!valuesMatch)
          {
            pointPassed = false;
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
            PointScalar x = inputPointsViewHost(pointOrdinal,0);
            PointScalar y = (spaceDim > 1) ? PointScalar(inputPointsViewHost(pointOrdinal,1)) : PointScalar(-2.0);
            PointScalar z = (spaceDim > 2) ? PointScalar(inputPointsViewHost(pointOrdinal,2)) : PointScalar(-2.0);
            
            if (spaceDim == 1)
              out << "values for "  << x  << " differ for field ordinal " << fieldOrdinalStandard;
            else if (spaceDim == 2)
              out << "values for ("  << x << "," << y << ") differ for field ordinal " << fieldOrdinalStandard;
            else
              out << "values for ("  << x << "," << y << "," << z << ") differ for field ordinal " << fieldOrdinalStandard;
            out << ": expected " << standardValue << "; actual " << derivedValue;
            out << " (diff: " << standardValue-derivedValue << ")" << std::endl;
            success = false;
          }
        }
        else // vector-valued
        {
          bool valuesMatch = true;
          int dkcard = standardOutputView.extent_int(2);
          for (int d=0; d<dkcard; d++)
          {
            OutputScalar standardValue = standardOutputViewHost(fieldOrdinalStandard,pointOrdinal,d);
            OutputScalar derivedValue  =  derivedOutputViewHost(fieldOrdinalDerived, pointOrdinal,d);
            
            TEUCHOS_TEST_FLOATING_EQUALITY(standardValue, derivedValue, tol, out, valuesMatch);
          }
          
          if (!valuesMatch)
          {
            pointPassed = false;
            PointScalar x = inputPointsViewHost(pointOrdinal,0);
            PointScalar y = (spaceDim > 1) ? PointScalar(inputPointsViewHost(pointOrdinal,1)) : PointScalar(-2.0);
            PointScalar z = (spaceDim > 2) ? PointScalar(inputPointsViewHost(pointOrdinal,2)) : PointScalar(-2.0);
            
            if (spaceDim == 1)
              out << "values for "  << x  << " differ for field ordinal " << fieldOrdinalStandard;
            else if (spaceDim == 2)
              out << "values for ("  << x << "," << y << ") differ for field ordinal " << fieldOrdinalStandard;
            else
              out << "values for ("  << x << "," << y << "," << z << ") differ for field ordinal " << fieldOrdinalStandard;
            out << ": expected (";
            for (int d=0; d<dkcard; d++)
            {
              out << standardOutputViewHost(fieldOrdinalStandard,pointOrdinal,d);
              if (d<dkcard-1) out << ",";
            }
            out << "); actual was (";
            for (int d=0; d<dkcard; d++)
            {
              out << derivedOutputViewHost(fieldOrdinalStandard,pointOrdinal,d);
              if (d<dkcard-1) out << ",";
            }
            out << ")" << std::endl;
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
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    
    std::vector<Intrepid2::EOperator> operators = {{OPERATOR_VALUE, OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5, OPERATOR_D6, OPERATOR_D7, OPERATOR_D8, OPERATOR_D9, OPERATOR_D10}};
    for (auto op : operators)
    {
      testHierarchicalHGRAD_LINE_MatchesAnalyticValues<ExecSpace,OutputScalar,PointScalar>(op, tol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( AnalyticPolynomialsMatch, Hierarchical_LineBasisDerivativesAgree, OutputScalar, PointScalar )
  {
    const int maxPolyOrder = 10;
    const double tol = TEST_TOLERANCE_TIGHT;
    
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using HierarchicalBasisFamily = HierarchicalBasisFamily<ExecSpace,OutputScalar,PointScalar>;
    
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
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using DerivedNodalBasisFamily  = Intrepid2::DerivedNodalBasisFamily<ExecSpace,OutputScalar,PointScalar>;
    using StandardNodalBasisFamily = Intrepid2::NodalBasisFamily       <ExecSpace,OutputScalar,PointScalar>;
    
    const double tol = TEST_TOLERANCE_TIGHT;
    
    shards::CellTopology lineTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >() );
    shards::CellTopology quadTopo = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<> >() );
    shards::CellTopology hexTopo = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<> >() );
    
    std::vector<EOperator> operators_dk = {OPERATOR_D1,OPERATOR_D2,OPERATOR_D3,OPERATOR_D4,OPERATOR_D5,OPERATOR_D6,OPERATOR_D7,OPERATOR_D8,OPERATOR_D9,OPERATOR_D10};
    
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
    
    out << "Running 3D nodal hexahedral basis comparison tests…\n";
    runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(polyOrder_3D, hexTopo, {FUNCTION_SPACE_HVOL}, {OPERATOR_VALUE}, tol, out, success);
    runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(polyOrder_3D, hexTopo, {FUNCTION_SPACE_HGRAD}, {OPERATOR_VALUE,OPERATOR_GRAD}, tol, out, success);
    runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(polyOrder_3D, hexTopo, {FUNCTION_SPACE_HGRAD}, operators_dk, tol, out, success);
  }
                                                        
  INTREPID2_OUTPUTSCALAR_POINTSCALAR_TEST_INSTANT( AnalyticPolynomialsMatch, Hierarchical_HGRAD_LINE )
  INTREPID2_OUTPUTSCALAR_POINTSCALAR_TEST_INSTANT( AnalyticPolynomialsMatch, Hierarchical_LineBasisDerivativesAgree )
  INTREPID2_OUTPUTSCALAR_POINTSCALAR_TEST_INSTANT( AnalyticPolynomialsMatch, HierarchicalNodalComparisons )
} // namespace
