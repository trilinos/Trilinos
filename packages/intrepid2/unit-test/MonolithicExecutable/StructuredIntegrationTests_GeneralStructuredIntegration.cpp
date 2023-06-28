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
void testStructuredIntegration(int meshWidth, int polyOrder, int worksetSize,
                               const FormulationChoice &formulation,
                               const double &relTol, const double &absTol,
                               Teuchos::FancyOStream &out, bool &success)
{
  // compare the general integration in StructuredAssembly.hpp, which takes two bases, two function spaces, and two ops,
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
  int numOps; // can be 1 or 2
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
  }
    
  double flopCountIntegration = 0, flopCountJacobian = 0;
  auto generalIntegrals = performStructuredAssembly<Scalar,BasisFamily>(geometry, worksetSize,
                                                                        polyOrder, fs, op1,
                                                                        polyOrder, fs, op1,
                                                                        flopCountIntegration, flopCountJacobian);
  if (numOps == 2)
  {
    auto generalIntegrals2 = performStructuredAssembly<Scalar,BasisFamily>(geometry, worksetSize,
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
  
  auto specificIntegrals = performStructuredQuadrature<Scalar, BasisFamily>(formulation, geometry, polyOrder, worksetSize, flopCountIntegration, flopCountJacobian);
    
  out << "Comparing new general standard assembly implementation to previous formulation-specific integration path…\n";
  testFloatingEquality3(generalIntegrals, specificIntegrals, relTol, absTol, out, success, "general integral", "specific formulation integral");
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(StructuredIntegration, GeneralStructuredIntegration, FormulationTag, DimTag, PolyOrderTag)
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
  
  testStructuredIntegration<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, polyOrder, worksetSize, formulation,
                                                                                        relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(StructuredIntegration, GeneralStructuredIntegration, PoissonFormulation, D1, P1)
TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(StructuredIntegration, GeneralStructuredIntegration, PoissonFormulation, D2, P3)
TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(StructuredIntegration, GeneralStructuredIntegration, PoissonFormulation, D3, P3)

} // anonymous namespace
