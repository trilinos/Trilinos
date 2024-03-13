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
// Questions? Contact Mauro Perego  (mperego@sandia.gov) or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file
    \brief  Tests against structured integration facilities - nonsymmetric tests (potentially different function spaces and ops).
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
#include <Intrepid2_HGRAD_TET_Cn_FEM.hpp>
#include <Intrepid2_HGRAD_HEX_Cn_FEM.hpp>
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
void testStandardVersusStructuredIntegration(const int &meshWidth, const int &worksetSize,
                                             const EFunctionSpace &fs1, const EOperator &op1, const int &p1,
                                             const EFunctionSpace &fs2, const EOperator &op2, const int &p2,
                                             const double &relTol, const double &absTol,
                                             Teuchos::FancyOStream &out, bool &success)
{
  // compare the integration in StructuredAssembly.hpp with that in StandardAssembly.hpp
  
  using namespace std;
  
  Kokkos::Array<int,spaceDim> gridCellCounts;
  for (int d=0; d<spaceDim; d++)
  {
    gridCellCounts[d] = meshWidth;
  }
  
  auto geometry = getMesh<PointScalar, spaceDim, DeviceType>(Standard, gridCellCounts);
  shards::CellTopology cellTopo = geometry.cellTopology();
  
  double flopCountIntegration = 0, flopCountJacobian = 0;
  auto structuredIntegrals = performStructuredAssembly<Scalar,BasisFamily>(geometry, worksetSize,
                                                                           p1, fs1, op1,
                                                                           p2, fs2, op2,
                                                                           flopCountIntegration, flopCountJacobian);
  
  auto standardIntegrals = performStandardAssembly<Scalar,BasisFamily>(geometry, worksetSize,
                                                                       p1, fs1, op1,
                                                                       p2, fs2, op2,
                                                                       flopCountIntegration, flopCountJacobian);
    
  out << "Comparing general standard assembly to structured integration path…\n";
  testFloatingEquality3(standardIntegrals, structuredIntegrals, relTol, absTol, out, success, "standard integral", "structured formulation integral");
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredIntegration, StructuredVersusStandard_D1_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 2;
  const int spaceDim = 1;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusStructuredIntegration<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredIntegration, StructuredVersusStandard_D1_P1_P2, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 2;
  const int spaceDim = 1;
  const int p1 = 1;
  const int p2 = 2;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusStructuredIntegration<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredIntegration, StructuredVersusStandard_D1_P2_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 2;
  const int spaceDim = 1;
  const int p1 = 2;
  const int p2 = 1;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusStructuredIntegration<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredIntegration, StructuredVersusStandard_D2_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 2;
  const int spaceDim = 2;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusStructuredIntegration<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredIntegration, StructuredVersusStandard_D2_P1_P2, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 2;
  const int spaceDim = 2;
  const int p1 = 1;
  const int p2 = 2;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusStructuredIntegration<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredIntegration, StructuredVersusStandard_D2_P2_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 2;
  const int spaceDim = 2;
  const int p1 = 2;
  const int p2 = 1;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusStructuredIntegration<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredIntegration, StructuredVersusStandard_D3_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 2;
  const int spaceDim = 3;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusStructuredIntegration<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredIntegration, StructuredVersusStandard_D3_P1_P2, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 3;
  const int p1 = 1;
  const int p2 = 2;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusStructuredIntegration<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredIntegration, StructuredVersusStandard_D3_P2_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 3;
  const int p1 = 2;
  const int p2 = 1;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusStructuredIntegration<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

// asymmetric tests (mostly -- a couple symmetric ones tossed in as sanity checks on the test itself)

// 1D tests: H(grad) and H(vol) bases defined
// p1, p1:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D1_P1_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D1_P1_P1, HGRAD, VALUE, HGRAD, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D1_P1_P1, HVOL,  VALUE, HGRAD, VALUE)
// p1, p2:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D1_P1_P2, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D1_P1_P2, HGRAD, VALUE, HGRAD, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D1_P1_P2, HVOL,  VALUE, HGRAD, VALUE)
// p2, p1:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D1_P2_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D1_P2_P1, HGRAD, VALUE, HGRAD, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D1_P2_P1, HVOL,  VALUE, HGRAD, VALUE)

// 2D tests: curls of H(curl) are scalars.
// p1, p1:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P1_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P1_P1, HGRAD, GRAD,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P1_P1, HGRAD, GRAD,  HCURL, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P1_P1, HGRAD, VALUE, HDIV,  DIV)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P1_P1, HGRAD, VALUE, HCURL, CURL)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P1_P1, HDIV,  DIV,   HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P1_P1, HCURL, CURL,  HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P1_P1, HVOL,  VALUE, HGRAD, VALUE)
// p2, p1:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P2_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P2_P1, HGRAD, GRAD,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P2_P1, HGRAD, GRAD,  HCURL, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P2_P1, HGRAD, VALUE, HDIV,  DIV)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P2_P1, HGRAD, VALUE, HCURL, CURL)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P2_P1, HDIV,  DIV,   HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P2_P1, HCURL, CURL,  HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P2_P1, HVOL,  VALUE, HGRAD, VALUE)
// p1, p2:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P1_P2, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P1_P2, HGRAD, GRAD,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P1_P2, HGRAD, GRAD,  HCURL, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P1_P2, HGRAD, VALUE, HDIV,  DIV)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P1_P2, HGRAD, VALUE, HCURL, CURL)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P1_P2, HDIV,  DIV,   HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P1_P2, HCURL, CURL,  HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D2_P1_P2, HVOL,  VALUE, HGRAD, VALUE)

// 3D tests: curls of H(curl) are vectors
// p1, p1:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P1_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P1_P1, HGRAD, GRAD,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P1_P1, HGRAD, GRAD,  HCURL, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P1_P1, HGRAD, GRAD,  HCURL, CURL)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P1_P1, HGRAD, VALUE, HDIV,  DIV)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P1_P1, HDIV,  DIV,   HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P1_P1, HCURL, CURL,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P1_P1, HVOL,  VALUE, HGRAD, VALUE)
// p2, p1:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P2_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P2_P1, HGRAD, GRAD,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P2_P1, HGRAD, GRAD,  HCURL, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P2_P1, HGRAD, GRAD,  HCURL, CURL)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P2_P1, HGRAD, VALUE, HDIV,  DIV)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P2_P1, HDIV,  DIV,   HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P2_P1, HCURL, CURL,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P2_P1, HVOL,  VALUE, HGRAD, VALUE)
// p1, p2:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P1_P2, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P1_P2, HGRAD, GRAD,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P1_P2, HGRAD, GRAD,  HCURL, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P1_P2, HGRAD, GRAD,  HCURL, CURL)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P1_P2, HGRAD, VALUE, HDIV,  DIV)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P1_P2, HDIV,  DIV,   HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P1_P2, HCURL, CURL,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredIntegration, StructuredVersusStandard_D3_P1_P2, HVOL,  VALUE, HGRAD, VALUE)


} // anonymous namespace
