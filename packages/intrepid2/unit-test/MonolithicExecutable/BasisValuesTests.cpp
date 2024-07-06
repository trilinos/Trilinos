// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "Intrepid2_TransformedBasisValues.hpp"
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
    using PointScalar  = typename Basis::PointValueType;
    using WeightScalar = typename Basis::OutputValueType;
    using DeviceType   = typename Basis::DeviceType;
    DefaultCubatureFactory cub_factory;
    auto cellTopoKey = basis.getBaseCellTopology().getKey();
    
    using Cubature       = Intrepid2::Cubature<DeviceType, PointScalar, WeightScalar>;
    using CubatureTensor = Intrepid2::CubatureTensor<DeviceType, PointScalar, WeightScalar>;
    using CubatureDirect = Intrepid2::CubatureDirect<DeviceType, PointScalar, WeightScalar>;
    
    const int quadratureDegree = 2*basis.getDegree();
    Teuchos::RCP<Cubature> lineTopoQuadrature = cub_factory.create<DeviceType, PointScalar, WeightScalar>(shards::Line<>::key, quadratureDegree);
    Teuchos::RCP<Cubature> baseTopoQuadrature = cub_factory.create<DeviceType, PointScalar, WeightScalar>(cellTopoKey, quadratureDegree);
    
    Teuchos::RCP<Intrepid2::Cubature<DeviceType, PointScalar, WeightScalar>> quadrature;
    
    const int numTensorialExtrusions = basis.getNumTensorialExtrusions();
    if (numTensorialExtrusions == 0)
    {
      quadrature = baseTopoQuadrature;
    }
    else
    {
      const CubatureDirect* baseTopoQuadratureDirect = dynamic_cast<CubatureDirect*>(baseTopoQuadrature.get());
      const CubatureDirect* lineTopoQuadratureDirect = dynamic_cast<CubatureDirect*>(lineTopoQuadrature.get());
      
      Teuchos::RCP<CubatureTensor> tensorCubature = Teuchos::rcp( new CubatureTensor(*baseTopoQuadratureDirect, *lineTopoQuadratureDirect));
      
      for (int d=1; d<numTensorialExtrusions; d++)
      {
        tensorCubature = Teuchos::rcp( new CubatureTensor(*tensorCubature, *lineTopoQuadratureDirect) );
      }
      
      quadrature = tensorCubature;
    }
    
    ordinal_type numRefPoints = quadrature->getNumPoints();
    const int spaceDim = basis.getBaseCellTopology().getDimension() + basis.getNumTensorialExtrusions();
    
    auto tensorPoints  = quadrature->allocateCubaturePoints();
    auto tensorWeights = quadrature->allocateCubatureWeights();
    
    quadrature->getCubature(tensorPoints, tensorWeights);
    
    using HostExecSpace = Kokkos::HostSpace::execution_space;
    TensorPoints<PointScalar,HostExecSpace> pointsHost(tensorPoints);
  
    using PointViewType = Kokkos::DynRankView<PointScalar>;
    using WeightViewType = Kokkos::DynRankView<WeightScalar>;
    PointViewType points("quadrature points ref cell - view", numRefPoints, spaceDim);
    WeightViewType weights("quadrature weights ref cell - view", numRefPoints);
    
    // copy from tensorPoints/Weights to points/weights
    Kokkos::parallel_for(numRefPoints, KOKKOS_LAMBDA(const int pointOrdinal)
    {
      weights(pointOrdinal) = tensorWeights(pointOrdinal);
      for (ordinal_type d=0; d<spaceDim; d++)
      {
        points(pointOrdinal,d) = tensorPoints(pointOrdinal,d);
      }
    });
    Kokkos::fence();
    
    auto hostPoints  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), points);
    auto hostWeights = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), weights);

    printFunctor2(points, out, "points being tested");
    printFunctor2(tensorPoints, out, "tensorPoints");
        
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
    BasisValues<double,DefaultTestDeviceType> emptyBasisValues;
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
    using DeviceType = Intrepid2::DefaultTestDeviceType;
    using Basis = HierarchicalBasisFamily<DeviceType>::HGRAD_QUAD;
    
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
    using DeviceType = Intrepid2::DefaultTestDeviceType;
    using Basis = HierarchicalBasisFamily<DeviceType>::HGRAD_TRI;
    
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
    using DeviceType = Intrepid2::DefaultTestDeviceType;
    using Basis = HierarchicalBasisFamily<DeviceType>::HGRAD_HEX;
    
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
    using DeviceType = Intrepid2::DefaultTestDeviceType;
    using Basis = HierarchicalBasisFamily<DeviceType>::HGRAD_TET;
    
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
    using DeviceType = Intrepid2::DefaultTestDeviceType;
    using Basis = HierarchicalBasisFamily<DeviceType>::HDIV_QUAD;
    
    // for now, the BasisValues path only supports the standard exact-sequence operators
//    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_DIV};
    std::vector<EOperator> opsToTest {OPERATOR_VALUE};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<2; polyOrder++)
    {
      Basis basis(polyOrder);
      testGetValuesEquality(basis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisValues, NodalHDIV_TRI )
  {
    using DeviceType = Intrepid2::DefaultTestDeviceType;
    using Basis = NodalBasisFamily<DeviceType>::HDIV_TRI; // Hierarchical basis family does not yet support HDIV_TRI
    
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
    using DeviceType = Intrepid2::DefaultTestDeviceType;
    using Basis = HierarchicalBasisFamily<DeviceType>::HDIV_HEX;
    
    // for now, the BasisValues path only supports the standard exact-sequence operators
//    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_DIV};
    std::vector<EOperator> opsToTest {OPERATOR_VALUE};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<2; polyOrder++)
    {
      Basis basis(polyOrder);
      testGetValuesEquality(basis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisValues, NodalHDIV_TET )
  {
    using DeviceType = Intrepid2::DefaultTestDeviceType;
    using Basis = NodalBasisFamily<DeviceType>::HDIV_TET;  // Hierarchical basis family does not yet support HDIV_TET
    
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
    using DeviceType = Intrepid2::DefaultTestDeviceType;
    using Basis = HierarchicalBasisFamily<DeviceType>::HCURL_QUAD;
    
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
    using DeviceType = Intrepid2::DefaultTestDeviceType;
    using Basis = NodalBasisFamily<DeviceType>::HCURL_TRI;  // Hierarchical basis family does not yet support HCURL_TRI
    
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
    using DeviceType = Intrepid2::DefaultTestDeviceType;
    using Basis = HierarchicalBasisFamily<DeviceType>::HCURL_HEX;
    
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

  TEUCHOS_UNIT_TEST( BasisValues, HierarchicalHVOL_QUAD )
  {
    using DeviceType = Intrepid2::DefaultTestDeviceType;
    using Basis = HierarchicalBasisFamily<DeviceType>::HVOL_QUAD;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      Basis basis(polyOrder);
      testGetValuesEquality(basis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisValues, HierarchicalHVOL_HEX )
  {
    using DeviceType = Intrepid2::DefaultTestDeviceType;
    using Basis = HierarchicalBasisFamily<DeviceType>::HVOL_HEX;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE};
    
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
    using DeviceType = Intrepid2::DefaultTestDeviceType;
    using Basis = NodalBasisFamily<DeviceType>::HCURL_TET;  // Hierarchical basis family does not yet support HCURL_TET
    
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
    using DeviceType = Intrepid2::DefaultTestDeviceType;
    using Basis = NodalBasisFamily<DeviceType>::HGRAD_LINE;
    
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
