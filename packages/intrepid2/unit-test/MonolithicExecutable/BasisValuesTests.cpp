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

/** \file   BasisValuesTests.cpp
    \brief  Tests to verify that the version of Basis::getValues() that takes the new BasisValues and TensorPoints arguments produces the same results as the one that takes raw Kokkos DynRankViews as arguments.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_CellGeometry.hpp"
#include "Intrepid2_CellGeometryTestUtils.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"
#include "Intrepid2_ProjectedGeometry.hpp"
#include "Intrepid2_ProjectedGeometryExamples.hpp"
#include "Intrepid2_ScalarView.hpp"
#include "Intrepid2_TransformedVectorData.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_TestUtils.hpp"

#include <Kokkos_Core.hpp>

using namespace Intrepid2;

namespace
{
  using namespace Intrepid2;
  
  //! tests that basis returns the same values for the getValues() that takes raw Kokkos containers as arguments and for the getValues() that takes
  //! BasisValues and TensorPoints containers as arguments.
  template<class Basis>
  void testGetValuesEquality(Basis &basis, const std::vector<EOperator> &opsToTest,
                             const double relTol, const double absTol, Teuchos::FancyOStream &out, bool &success)
  {
    // get quadrature points for integrating up to polyOrder
    const int quadratureDegree = 2*basis.getDegree();
    using ExecutionSpace = typename Basis::ExecutionSpace;
    using PointScalar = typename Basis::PointValueType;
    using WeightScalar = typename Basis::OutputValueType;
    DefaultCubatureFactory cub_factory;
    auto cellTopoKey = basis.getBaseCellTopology().getKey();
    auto quadrature = cub_factory.create<ExecutionSpace, PointScalar, WeightScalar>(cellTopoKey, quadratureDegree);
    ordinal_type numRefPoints = quadrature->getNumPoints();
    const int spaceDim = basis.getBaseCellTopology().getDimension();
    ViewType<PointScalar> points = ViewType<PointScalar>("quadrature points ref cell", numRefPoints, spaceDim);
    ViewType<WeightScalar> weights = ViewType<WeightScalar>("quadrature weights ref cell", numRefPoints);
    quadrature->getCubature(points, weights);
    
    TensorPoints<PointScalar> tensorPoints;
    TensorData<WeightScalar>  tensorWeights;
    
    using HostPointViewType = Kokkos::DynRankView<PointScalar,Kokkos::HostSpace>;
    using HostWeightViewType = Kokkos::DynRankView<WeightScalar,Kokkos::HostSpace>;
    auto hostQuadrature = cub_factory.create<Kokkos::HostSpace, PointScalar, WeightScalar>(cellTopoKey, quadratureDegree);
    HostPointViewType hostPoints("quadrature points ref cell - host", numRefPoints, spaceDim);
    HostWeightViewType hostWeights("quadrature weights ref cell - host", numRefPoints);
    hostQuadrature->getCubature(hostPoints, hostWeights);
    
    using CubatureTensorType = CubatureTensor<ExecutionSpace,PointScalar,WeightScalar>;
    CubatureTensorType* tensorQuadrature = dynamic_cast<CubatureTensorType*>(quadrature.get());

    if (tensorQuadrature)
    {
      tensorPoints  = tensorQuadrature->allocateCubaturePoints();
      tensorWeights = tensorQuadrature->allocateCubatureWeights();
      tensorQuadrature->getCubature(tensorPoints, tensorWeights);
    }
    else
    {
      std::vector<ViewType<PointScalar>> pointComponents {points};
      tensorPoints = TensorPoints<PointScalar>(pointComponents);
      Data<WeightScalar> weightData(weights);
      std::vector<Data<WeightScalar>> weightComponents {weightData};
      tensorWeights = TensorData<WeightScalar>(weightComponents);
    }
    
    out << "Points being tested:\n";
    for (int pointOrdinal=0; pointOrdinal<numRefPoints; pointOrdinal++)
    {
      out << pointOrdinal << ": " << "(";
      for (int d=0; d<spaceDim; d++)
      {
        out << points(pointOrdinal,d);
        if (d < spaceDim-1) out << ",";
      }
      out << ")\n";
    }
    
    printFunctor2(tensorPoints, out, "tensorPoints");
    
    testFloatingEquality2(points,tensorPoints,  relTol, absTol, out, success, "points", "tensorPoints");
        
    auto hostBasisPtr = basis.getHostBasis();
        
    for (const auto &op : opsToTest)
    {
      auto basisValuesView = basis.allocateOutputView(numRefPoints, op);
      auto basisValues     = basis.allocateBasisValues(tensorPoints, op);
      
      auto hostBasisView   = hostBasisPtr->allocateOutputView(numRefPoints, op);
      
      basis.getValues(basisValuesView, points, op);
      basis.getValues(basisValues, tensorPoints, op);
      
      // copy basisValuesView to host for hostBasis comparison
      auto basisValuesViewHostMirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), basisValuesView);
      hostBasisPtr->getValues(hostBasisView, hostPoints, op);
      
      out << "Comparing getValues() results for " << EOperatorToString(op) << std::endl;
      
      TEST_EQUALITY(basisValues.rank(), basisValuesView.rank());

      if (basisValues.rank() == basisValuesView.rank())
      {
        if (basisValuesView.rank() == 2)
        {
          bool localSuccess = true;
          testFloatingEquality2(basisValuesView, basisValues, relTol, absTol, out, localSuccess, "DynRankView path", "BasisValues path");
          success = success && localSuccess;
          
          if (!localSuccess)
          {
            printFunctor2(basisValues,     out, "basisValues");
            printFunctor2(basisValuesView, out, "basisValuesView");
          }
          
          localSuccess = true;
          testFloatingEquality2(basisValuesView, basisValues, relTol, absTol, out, localSuccess, "DynRankView path - device", "DynRankView path - host");
          success = success && localSuccess;
          
          if (!localSuccess)
          {
            printFunctor2(hostBasisView,             out, "hostBasisView");
            printFunctor2(basisValuesViewHostMirror, out, "basisValuesViewHostMirror");
          }
        }
        else if (basisValuesView.rank() == 3)
        {
          bool localSuccess = true;
          testFloatingEquality3(basisValuesView, basisValues, relTol, absTol, out, localSuccess, "DynRankView path", "BasisValues path");
          success = success && localSuccess;
          
          if (!localSuccess)
          {
            printFunctor3(basisValues,     out, "basisValues");
            printFunctor3(basisValuesView, out, "basisValuesView");
          }
          
          localSuccess = true;
          testFloatingEquality3(basisValuesView, basisValues, relTol, absTol, out, localSuccess, "DynRankView path - device", "DynRankView path - host");
          success = success && localSuccess;
          
          if (!localSuccess)
          {
            printFunctor3(hostBasisView,             out, "hostBasisView");
            printFunctor3(basisValuesViewHostMirror, out, "basisValuesViewHostMirror");
          }
        }
        else
        {
          INTREPID2_TEST_FOR_EXCEPTION(true, std::logic_error, "Unexpected rank encountered for basisValuesView");
        }
      }
      else
      {
        out << "FAILURE: basisValues.rank() does not match basisValuesView.rank().\n";
        if (basisValues.rank() == 2)
        {
          printFunctor2(basisValues, out, "basisValues");
        }
        else if (basisValues.rank() == 3)
        {
          printFunctor3(basisValues, out, "basisValues");
        }
      }
    }
  }

  TEUCHOS_UNIT_TEST( BasisValues, DefaultConstructor )
  {
    // test of default-constructed basis values object.
    BasisValues<double> emptyBasisValues;
    TEST_EQUALITY(0, emptyBasisValues.rank());
    for (int d=0; d<8; d++)
    {
      TEST_EQUALITY(0, emptyBasisValues.extent(d));
      TEST_EQUALITY(0, emptyBasisValues.extent_int(d));
    }
    TEST_EQUALITY(0, emptyBasisValues.numFamilies());
  }


  TEUCHOS_UNIT_TEST( BasisValues, HierarchicalHGRAD_LINE )
  {
    using Basis = HierarchicalBasisFamily<>::HGRAD_LINE;
    
    // for now, the BasisValues path only supports the standard exact-sequence operators
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_GRAD};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      Basis basis(polyOrder);
      testGetValuesEquality(basis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisValues, HierarchicalHGRAD_LINE_DeviceType )
  {
    using DeviceType = Intrepid2::DefaultTestDeviceType;
    using Basis = HierarchicalBasisFamily<DeviceType>::HGRAD_LINE;
    
    // for now, the BasisValues path only supports the standard exact-sequence operators
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_GRAD};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      Basis basis(polyOrder);
      testGetValuesEquality(basis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisValues, HierarchicalHGRAD_QUAD )
  {
    using Basis = HierarchicalBasisFamily<>::HGRAD_QUAD;
    
    // for now, the BasisValues path only supports the standard exact-sequence operators
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_GRAD};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      Basis basis(polyOrder);
      testGetValuesEquality(basis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisValues, HierarchicalHGRAD_TRI )
  {
    using Basis = HierarchicalBasisFamily<>::HGRAD_TRI;
    
    // for now, the BasisValues path only supports the standard exact-sequence operators
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_GRAD};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      Basis basis(polyOrder);
      testGetValuesEquality(basis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisValues, HierarchicalHGRAD_HEX )
  {
    using Basis = HierarchicalBasisFamily<>::HGRAD_HEX;
    
    // for now, the BasisValues path only supports the standard exact-sequence operators
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_GRAD};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      Basis basis(polyOrder);
      testGetValuesEquality(basis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisValues, HierarchicalHGRAD_TET )
  {
    using Basis = HierarchicalBasisFamily<>::HGRAD_TET;
    
    // for now, the BasisValues path only supports the standard exact-sequence operators
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_GRAD};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      Basis basis(polyOrder);
      testGetValuesEquality(basis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisValues, HierarchicalHDIV_QUAD )
  {
    using Basis = HierarchicalBasisFamily<>::HDIV_QUAD;
    
    // for now, the BasisValues path only supports the standard exact-sequence operators
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_DIV};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      Basis basis(polyOrder);
      testGetValuesEquality(basis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisValues, NodalHDIV_TRI )
  {
    using Basis = NodalBasisFamily<>::HDIV_TRI; // Hierarchical basis family does not yet support HDIV_TRI
    
    // for now, the BasisValues path only supports the standard exact-sequence operators
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_DIV};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      Basis basis(polyOrder);
      testGetValuesEquality(basis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisValues, HierarchicalHDIV_HEX )
  {
    using Basis = HierarchicalBasisFamily<>::HDIV_HEX;
    
    // for now, the BasisValues path only supports the standard exact-sequence operators
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_DIV};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      Basis basis(polyOrder);
      testGetValuesEquality(basis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisValues, NodalHDIV_TET )
  {
    using Basis = NodalBasisFamily<>::HDIV_TET;  // Hierarchical basis family does not yet support HDIV_TET
    
    // for now, the BasisValues path only supports the standard exact-sequence operators
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_DIV};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      Basis basis(polyOrder);
      testGetValuesEquality(basis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisValues, HierarchicalHCURL_QUAD )
  {
    using Basis = HierarchicalBasisFamily<>::HCURL_QUAD;
    
    // for now, the BasisValues path only supports the standard exact-sequence operators
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_CURL};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      Basis basis(polyOrder);
      testGetValuesEquality(basis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisValues, NodalHCURL_TRI )
  {
    using Basis = NodalBasisFamily<>::HCURL_TRI;  // Hierarchical basis family does not yet support HCURL_TRI
    
    // for now, the BasisValues path only supports the standard exact-sequence operators
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_CURL};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      Basis basis(polyOrder);
      testGetValuesEquality(basis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisValues, HierarchicalHCURL_HEX )
  {
    using Basis = HierarchicalBasisFamily<>::HCURL_HEX;
    
    // for now, the BasisValues path only supports the standard exact-sequence operators
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_CURL};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      Basis basis(polyOrder);
      testGetValuesEquality(basis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisValues, NodalHCURL_TET )
  {
    using Basis = NodalBasisFamily<>::HCURL_TET;  // Hierarchical basis family does not yet support HCURL_TET
    
    // for now, the BasisValues path only supports the standard exact-sequence operators
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_CURL};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      Basis basis(polyOrder);
      testGetValuesEquality(basis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisValues, NodalHGRAD_LINE )
  {
    using Basis = NodalBasisFamily<>::HGRAD_LINE;
    
    // for now, the BasisValues path only supports the standard exact-sequence operators
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_GRAD};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      Basis basis(polyOrder);
      testGetValuesEquality(basis, opsToTest, relTol, absTol, out, success);
    }
  }
} // namespace
