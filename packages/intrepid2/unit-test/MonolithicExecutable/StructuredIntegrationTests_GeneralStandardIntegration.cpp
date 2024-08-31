// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Tests against structured integration facilities.
    \author Nathan V. Roberts
*/

#include "Kokkos_Core.hpp"

#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_UnitTestRepository.hpp>

#include <Intrepid2_CellGeometryTestUtils.hpp>
#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_DefaultCubatureFactory.hpp>
#include <Intrepid2_FunctionSpaceTools.hpp>
#include <Intrepid2_IntegrationTools.hpp>
#include <Intrepid2_Kernels.hpp>
#include <Intrepid2_NodalBasisFamily.hpp>
#include <Intrepid2_TensorArgumentIterator.hpp>
#include <Intrepid2_TestUtils.hpp>

#include <Intrepid2_CellGeometry.hpp>
#include "Intrepid2_Data.hpp"
#include "Intrepid2_TensorData.hpp"
#include "Intrepid2_TensorPoints.hpp"
#include "Intrepid2_TransformedBasisValues.hpp"
#include "Intrepid2_VectorData.hpp"

#include "Intrepid2_ScalarView.hpp"

#include "GRADGRADStandardAssembly.hpp"
#include "GRADGRADStructuredAssembly.hpp"
#include "H1StandardAssembly.hpp"
#include "H1StructuredAssembly.hpp"
#include "HDIVStandardAssembly.hpp"
#include "HDIVStructuredAssembly.hpp"
#include "HCURLStandardAssembly.hpp"
#include "HCURLStructuredAssembly.hpp"
#include "HVOLStandardAssembly.hpp"
#include "HVOLStructuredAssembly.hpp"

#include "StandardAssembly.hpp"
#include "StructuredAssembly.hpp"

#include "StructuredIntegrationTests_TagDefs.hpp"
#include "StructuredIntegrationTests_Utils.hpp"

namespace
{
  using namespace Intrepid2;

template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType>
void testStandardIntegration(int meshWidth, int polyOrder, int worksetSize,
                             const FormulationChoice &formulation,
                             const double &relTol, const double &absTol,
                             Teuchos::FancyOStream &out, bool &success)
{
  // compare the general integration in StandardAssembly.hpp, which takes two bases, two function spaces, and two ops,
  // with the specific implementations for (grad, grad), etc. formulations.
  
  using namespace std;
  
  Kokkos::Array<int,spaceDim> gridCellCounts;
  for (int d=0; d<spaceDim; d++)
  {
    gridCellCounts[d] = meshWidth;
  }
  
  auto geometry = getMesh<PointScalar, spaceDim, DeviceType>(Standard, gridCellCounts);
  shards::CellTopology cellTopo = geometry.cellTopology();
  
  EFunctionSpace fs;
  EOperator op1, op2;
  int numOps = 0; // can be 1 or 2
  Teuchos::RCP<Kokkos::Array<Scalar,spaceDim>> vectorWeight1, vectorWeight2;
  switch (formulation)
  {
    case Poisson:
      numOps = 1;
      op1 = EOperator::OPERATOR_GRAD;
      fs = EFunctionSpace::FUNCTION_SPACE_HGRAD;
      break;
    case Hgrad:
      numOps = 2;
      op1 = EOperator::OPERATOR_GRAD;
      op2 = EOperator::OPERATOR_VALUE;
      fs = EFunctionSpace::FUNCTION_SPACE_HGRAD;
      break;
    case Hdiv:
      numOps = 2;
      op1 = EOperator::OPERATOR_DIV;
      op2 = EOperator::OPERATOR_VALUE;
      fs = EFunctionSpace::FUNCTION_SPACE_HDIV;
      break;
    case Hcurl:
      numOps = 2;
      op1 = EOperator::OPERATOR_CURL;
      op2 = EOperator::OPERATOR_VALUE;
      fs = EFunctionSpace::FUNCTION_SPACE_HCURL;
      break;
    case L2:
      numOps = 1;
      op1 = EOperator::OPERATOR_VALUE;
      fs = EFunctionSpace::FUNCTION_SPACE_HDIV;
      break;
    case VectorWeightedPoisson:
      numOps = 1;
      op1 = EOperator::OPERATOR_GRAD;
      fs = EFunctionSpace::FUNCTION_SPACE_HGRAD;
      vectorWeight1 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
      vectorWeight2 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
      double weight = 1.0;
      for (int d=0; d<spaceDim; d++)
      {
        (*vectorWeight1)[d] = weight;
        weight /= 2.0;
      }
      
      weight = 0.5;
      for (int d=0; d<spaceDim; d++)
      {
        (*vectorWeight2)[d] = weight;
        weight *= 2.0;
      }
      break;
  }
    
  double flopCountIntegration = 0, flopCountJacobian = 0;
  auto generalIntegrals = performStandardAssembly<Scalar,BasisFamily>(geometry, worksetSize,
                                                                      polyOrder, fs, op1, vectorWeight1,
                                                                      polyOrder, fs, op1, vectorWeight2,
                                                                      flopCountIntegration, flopCountJacobian);
  if (numOps == 2)
  {
    auto generalIntegrals2 = performStandardAssembly<Scalar,BasisFamily>(geometry, worksetSize,
                                                                         polyOrder, fs, op2,
                                                                         polyOrder, fs, op2,
                                                                         flopCountIntegration, flopCountJacobian);
    
    using ExecutionSpace = typename DeviceType::execution_space;
    auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{generalIntegrals.extent_int(0),generalIntegrals.extent_int(1),generalIntegrals.extent_int(2)});
    Kokkos::parallel_for("sum integrals", policy,
    KOKKOS_LAMBDA (const int &C, const int &F1, const int &F2)
    {
      generalIntegrals(C,F1,F2) += generalIntegrals2(C,F1,F2);
    });
  }
  
  auto specificIntegrals = performStandardQuadrature<Scalar,BasisFamily>(formulation, geometry, polyOrder, worksetSize, flopCountIntegration, flopCountJacobian, vectorWeight1, vectorWeight2);
    
  out << "Comparing new general standard assembly implementation to previous formulation-specific integration path…\n";
  testFloatingEquality3(generalIntegrals, specificIntegrals, relTol, absTol, out, success, "general integral", "specific formulation integral");
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(StructuredIntegration, GeneralStandardIntegration, FormulationTag, DimTag, PolyOrderTag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 2;
  const int spaceDim = DimTag::spaceDim;
  const int polyOrder = PolyOrderTag::polyOrder;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const FormulationChoice formulation = FormulationTag::formulation;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardIntegration<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, polyOrder, worksetSize, formulation,
                                                                                      relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(StructuredIntegration, GeneralStandardIntegration, PoissonFormulation, D1, P1)
TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(StructuredIntegration, GeneralStandardIntegration, PoissonFormulation, D2, P3)
TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(StructuredIntegration, GeneralStandardIntegration, PoissonFormulation, D3, P3)

TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(StructuredIntegration, GeneralStandardIntegration, VectorWeightedPoissonFormulation, D1, P1)
TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(StructuredIntegration, GeneralStandardIntegration, VectorWeightedPoissonFormulation, D2, P3)
TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(StructuredIntegration, GeneralStandardIntegration, VectorWeightedPoissonFormulation, D3, P3)

} // anonymous namespace
