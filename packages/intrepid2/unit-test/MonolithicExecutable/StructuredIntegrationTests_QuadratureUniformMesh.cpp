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
    \brief  Tests against structured integration facilities - tests against uniform quadrature for symmetric formulations, comparing standard assembly to structured assembly.
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
  void testQuadratureHypercube(int meshWidth, int polyOrder, int worksetSize,
                               const FormulationChoice &formulation, const AlgorithmChoice &algorithm,
                               const double &relTol, const double &absTol,
                               Teuchos::FancyOStream &out, bool &success)
  {
    using namespace std;
    
    Kokkos::Array<int,spaceDim> gridCellCounts;
    for (int d=0; d<spaceDim; d++)
    {
      gridCellCounts[d] = meshWidth;
    }
    
    auto geometry = getMesh<PointScalar, spaceDim, DeviceType>(algorithm, gridCellCounts);
    double flopCountIntegration = 0, flopCountJacobian = 0;
    auto standardIntegrals = performStandardQuadrature<Scalar, BasisFamily>(formulation, geometry, polyOrder, worksetSize, flopCountIntegration, flopCountJacobian);
    
    auto structuredIntegrals = performStructuredQuadrature<Scalar, BasisFamily>(formulation, geometry, polyOrder, worksetSize, flopCountIntegration, flopCountJacobian);
    
    out << "Comparing standard Intrepid2 integration to new integration path…\n";
    testFloatingEquality3(standardIntegrals, structuredIntegrals, relTol, absTol, out, success, "standard Intrepid2 integral", "structured integral");
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredIntegration, QuadratureUniformMesh, FormulationTag, AlgorithmTag, DimTag, PolyOrderTag)
  {
    using DataScalar  = double;
    using PointScalar = double;
    const int meshWidth = 2;
    const int spaceDim = DimTag::spaceDim;
    const int polyOrder = PolyOrderTag::polyOrder;
    const int worksetSize = meshWidth;

    using DeviceType = DefaultTestDeviceType;
    using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
    
    const AlgorithmChoice algorithm = AlgorithmTag::algorithm;
    const FormulationChoice formulation = FormulationTag::formulation;
    
    double relTol = 1e-12;
    double absTol = 1e-12;
    
    testQuadratureHypercube<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, polyOrder, worksetSize, formulation, algorithm,
                                                                                        relTol, absTol, out, success);
  }

  // comparisons are to Standard algorithm, so we don't instantiate with Standard:
  // 1D, p=1 tests:
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, AffineTensorAlgorithm,    D1, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, NonAffineTensorAlgorithm, D1, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, AffineNonTensorAlgorithm, D1, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, UniformAlgorithm,         D1, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   AffineTensorAlgorithm,    D1, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   NonAffineTensorAlgorithm, D1, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   AffineNonTensorAlgorithm, D1, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   UniformAlgorithm,         D1, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      AffineTensorAlgorithm,    D1, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      NonAffineTensorAlgorithm, D1, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      AffineNonTensorAlgorithm, D1, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      UniformAlgorithm,         D1, P1)
  // 1D, p=2 tests:
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, AffineTensorAlgorithm,    D1, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, NonAffineTensorAlgorithm, D1, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, AffineNonTensorAlgorithm, D1, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, UniformAlgorithm,         D1, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   AffineTensorAlgorithm,    D1, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   NonAffineTensorAlgorithm, D1, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   AffineNonTensorAlgorithm, D1, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   UniformAlgorithm,         D1, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      AffineTensorAlgorithm,    D1, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      NonAffineTensorAlgorithm, D1, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      AffineNonTensorAlgorithm, D1, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      UniformAlgorithm,         D1, P2)
  // 1D, p=4 tests:
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, AffineTensorAlgorithm,    D1, P4)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, NonAffineTensorAlgorithm, D1, P4)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, AffineNonTensorAlgorithm, D1, P4)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, UniformAlgorithm,         D1, P4)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   AffineTensorAlgorithm,    D1, P4)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   NonAffineTensorAlgorithm, D1, P4)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   AffineNonTensorAlgorithm, D1, P4)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   UniformAlgorithm,         D1, P4)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      AffineTensorAlgorithm,    D1, P4)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      NonAffineTensorAlgorithm, D1, P4)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      AffineNonTensorAlgorithm, D1, P4)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      UniformAlgorithm,         D1, P4)

  // 2D, p=1 tests:
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, AffineTensorAlgorithm,    D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, NonAffineTensorAlgorithm, D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, AffineNonTensorAlgorithm, D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, UniformAlgorithm,         D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   AffineTensorAlgorithm,    D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   NonAffineTensorAlgorithm, D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   AffineNonTensorAlgorithm, D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   UniformAlgorithm,         D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    AffineTensorAlgorithm,    D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    NonAffineTensorAlgorithm, D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    AffineNonTensorAlgorithm, D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    UniformAlgorithm,         D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   AffineTensorAlgorithm,    D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   NonAffineTensorAlgorithm, D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   AffineNonTensorAlgorithm, D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   UniformAlgorithm,         D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      AffineTensorAlgorithm,    D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      NonAffineTensorAlgorithm, D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      AffineNonTensorAlgorithm, D2, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      UniformAlgorithm,         D2, P1)
  // 2D, p=2 tests:
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, AffineTensorAlgorithm,    D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, NonAffineTensorAlgorithm, D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, AffineNonTensorAlgorithm, D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, UniformAlgorithm,         D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   AffineTensorAlgorithm,    D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   NonAffineTensorAlgorithm, D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   AffineNonTensorAlgorithm, D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   UniformAlgorithm,         D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    AffineTensorAlgorithm,    D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    NonAffineTensorAlgorithm, D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    AffineNonTensorAlgorithm, D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    UniformAlgorithm,         D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   AffineTensorAlgorithm,    D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   NonAffineTensorAlgorithm, D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   AffineNonTensorAlgorithm, D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   UniformAlgorithm,         D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      AffineTensorAlgorithm,    D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      NonAffineTensorAlgorithm, D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      AffineNonTensorAlgorithm, D2, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      UniformAlgorithm,         D2, P2)
  // 2D, p=3 tests:
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, AffineTensorAlgorithm,    D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, NonAffineTensorAlgorithm, D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, AffineNonTensorAlgorithm, D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, UniformAlgorithm,         D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   AffineTensorAlgorithm,    D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   NonAffineTensorAlgorithm, D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   AffineNonTensorAlgorithm, D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   UniformAlgorithm,         D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    AffineTensorAlgorithm,    D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    NonAffineTensorAlgorithm, D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    AffineNonTensorAlgorithm, D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    UniformAlgorithm,         D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   AffineTensorAlgorithm,    D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   NonAffineTensorAlgorithm, D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   AffineNonTensorAlgorithm, D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   UniformAlgorithm,         D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      AffineTensorAlgorithm,    D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      NonAffineTensorAlgorithm, D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      AffineNonTensorAlgorithm, D2, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      UniformAlgorithm,         D2, P3)

  // 3D, p=1 tests:
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, AffineTensorAlgorithm,    D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, NonAffineTensorAlgorithm, D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, AffineNonTensorAlgorithm, D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, UniformAlgorithm,         D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   AffineTensorAlgorithm,    D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   NonAffineTensorAlgorithm, D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   AffineNonTensorAlgorithm, D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   UniformAlgorithm,         D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    AffineTensorAlgorithm,    D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    NonAffineTensorAlgorithm, D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    AffineNonTensorAlgorithm, D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    UniformAlgorithm,         D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   AffineTensorAlgorithm,    D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   NonAffineTensorAlgorithm, D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   AffineNonTensorAlgorithm, D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   UniformAlgorithm,         D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      AffineTensorAlgorithm,    D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      NonAffineTensorAlgorithm, D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      AffineNonTensorAlgorithm, D3, P1)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      UniformAlgorithm,         D3, P1)
  // 3D, p=2 tests:
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, AffineTensorAlgorithm,    D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, NonAffineTensorAlgorithm, D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, AffineNonTensorAlgorithm, D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, UniformAlgorithm,         D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   AffineTensorAlgorithm,    D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   NonAffineTensorAlgorithm, D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   AffineNonTensorAlgorithm, D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   UniformAlgorithm,         D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    AffineTensorAlgorithm,    D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    NonAffineTensorAlgorithm, D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    AffineNonTensorAlgorithm, D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    UniformAlgorithm,         D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   AffineTensorAlgorithm,    D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   NonAffineTensorAlgorithm, D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   AffineNonTensorAlgorithm, D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   UniformAlgorithm,         D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      AffineTensorAlgorithm,    D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      NonAffineTensorAlgorithm, D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      AffineNonTensorAlgorithm, D3, P2)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      UniformAlgorithm,         D3, P2)
  // 3D, p=3 tests:
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, AffineTensorAlgorithm,    D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, NonAffineTensorAlgorithm, D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, AffineNonTensorAlgorithm, D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, PoissonFormulation, UniformAlgorithm,         D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   AffineTensorAlgorithm,    D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   NonAffineTensorAlgorithm, D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   AffineNonTensorAlgorithm, D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HgradFormulation,   UniformAlgorithm,         D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    AffineTensorAlgorithm,    D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    NonAffineTensorAlgorithm, D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    AffineNonTensorAlgorithm, D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HdivFormulation,    UniformAlgorithm,         D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   AffineTensorAlgorithm,    D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   NonAffineTensorAlgorithm, D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   AffineNonTensorAlgorithm, D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, HcurlFormulation,   UniformAlgorithm,         D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      AffineTensorAlgorithm,    D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      NonAffineTensorAlgorithm, D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      AffineNonTensorAlgorithm, D3, P3)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, QuadratureUniformMesh, L2Formulation,      UniformAlgorithm,         D3, P3)
} // anonymous namespace
