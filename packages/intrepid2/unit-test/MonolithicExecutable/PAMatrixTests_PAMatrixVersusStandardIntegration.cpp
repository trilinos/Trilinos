// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Comparing PAMatrix versus standard matrix assembly - nonsymmetric tests (potentially different function spaces and ops).
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

#include "StandardAssembly.hpp"
#include "PAMatrixAssembly.hpp"

#include "StructuredIntegrationTests_TagDefs.hpp"
#include "StructuredIntegrationTests_Utils.hpp"

namespace
{
  using namespace Intrepid2;

template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType>
void testStandardVersusPAMatrix(const int &meshWidth, const int &worksetSize,
                                const EFunctionSpace &fs1, const EOperator &op1, const int &p1, Teuchos::RCP< Kokkos::Array<PointScalar,spaceDim> > vectorWeight1,
                                const EFunctionSpace &fs2, const EOperator &op2, const int &p2, Teuchos::RCP< Kokkos::Array<PointScalar,spaceDim> > vectorWeight2,
                                const double &relTol, const double &absTol, Teuchos::FancyOStream &out, bool &success)
{
  // compare the matrix-free action (etc.) specified by PAMatrixAssembly.hpp with the explicit integration in StandardAssembly.hpp
  
  using namespace std;
  
  Kokkos::Array<int,spaceDim> gridCellCounts;
  for (int d=0; d<spaceDim; d++)
  {
    gridCellCounts[d] = meshWidth;
  }
  
  auto geometry = getMesh<PointScalar, spaceDim, DeviceType>(Standard, gridCellCounts);
  shards::CellTopology cellTopo = geometry.cellTopology();
  
  double flopCountIntegration = 0, flopCountJacobian = 0;
  auto paMatrix = constructPAMatrix<Scalar,BasisFamily>(geometry,
                                                        p1, fs1, op1, vectorWeight1,
                                                        p2, fs2, op2, vectorWeight2);
  
  auto standardIntegrals = performStandardAssembly<Scalar,BasisFamily>(geometry, worksetSize,
                                                                       p1, fs1, op1, vectorWeight1,
                                                                       p2, fs2, op2, vectorWeight2,
                                                                       flopCountIntegration, flopCountJacobian);
  
  auto fullMatrix = paMatrix.allocateMatrixStorage();
  paMatrix.assemble(fullMatrix);
  
  out << "Comparing general standard assembly to PAMatrix full-assembly path…\n";
  testFloatingEquality3(standardIntegrals, fullMatrix, relTol, absTol, out, success, "standard integral", "PAMatrix full assembly");
}

template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType>
void testStandardVersusPAMatrix(const int &meshWidth, const int &worksetSize,
                                             const EFunctionSpace &fs1, const EOperator &op1, const int &p1,
                                             const EFunctionSpace &fs2, const EOperator &op2, const int &p2,
                                             const double &relTol, const double &absTol,
                                             Teuchos::FancyOStream &out, bool &success)
{
  testStandardVersusPAMatrix<Scalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, worksetSize,
                                                                                                  fs1, op1, p1, Teuchos::null,
                                                                                                  fs2, op2, p2, Teuchos::null,
                                                                                                  relTol, absTol, out, success);
}

// MARK: PAMatrixVersusStandard_D1_P1_P1
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandard_D1_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
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
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

// MARK: PAMatrixVersusStandard_D1_P1_P2
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandard_D1_P1_P2, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
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
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandard_D1_P2_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
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
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandard_D2_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
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
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandard_D2_P1_P2, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
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
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandard_D2_P2_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
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
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandard_D3_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
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
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandard_D3_P1_P2, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
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
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandard_D3_P2_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
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
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandardVectorWeighted_D1_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 1;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  auto vectorWeight1 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  auto vectorWeight2 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  
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

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandardVectorWeighted_D2_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 2;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  auto vectorWeight1 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  auto vectorWeight2 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  
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

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandardVectorWeighted_D2_P2_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 2;
  const int p1 = 2;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  auto vectorWeight1 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  auto vectorWeight2 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  
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

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandardScalarAgainstVectorDotVector_D1_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 1;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  Teuchos::RCP<Kokkos::Array<double,spaceDim> > vectorWeight1; // no vector weight on scalar term
  auto vectorWeight2 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  
  double weight = 1.0;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight2)[d] = weight;
    weight /= 2.0;
  }

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandardScalarAgainstVectorDotVector_D2_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 2;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  Teuchos::RCP<Kokkos::Array<double,spaceDim> > vectorWeight1; // no vector weight on scalar term
  auto vectorWeight2 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  
  double weight = 1.0;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight2)[d] = weight;
    weight /= 2.0;
  }

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandardScalarAgainstVectorDotVector_D3_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 3;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  Teuchos::RCP<Kokkos::Array<double,spaceDim> > vectorWeight1; // no vector weight on scalar term
  auto vectorWeight2 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  
  double weight = 1.0;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight2)[d] = weight;
    weight /= 2.0;
  }

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandardVectorDotVectorAgainstScalar_D1_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 1;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  auto vectorWeight1 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  double weight = 1.0;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight1)[d] = weight;
    weight /= 2.0;
  }
  Teuchos::RCP<Kokkos::Array<double,spaceDim> > vectorWeight2; // no vector weight on scalar term
  
  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandardVectorDotVectorAgainstScalar_D2_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 2;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  auto vectorWeight1 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  double weight = 1.0;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight1)[d] = weight;
    weight /= 2.0;
  }
  Teuchos::RCP<Kokkos::Array<double,spaceDim> > vectorWeight2; // no vector weight on scalar term
  
  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandardVectorDotVectorAgainstScalar_D3_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 3;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  auto vectorWeight1 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  double weight = 1.0;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight1)[d] = weight;
    weight /= 2.0;
  }
  Teuchos::RCP<Kokkos::Array<double,spaceDim> > vectorWeight2; // no vector weight on scalar term
  
  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, PAMatrixVersusStandardVectorWeighted_D3_P2_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 3;
  const int p1 = 2;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  auto vectorWeight1 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  auto vectorWeight2 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  
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

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testStandardVersusPAMatrix<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

// asymmetric tests (mostly -- a couple symmetric ones tossed in as sanity checks on the test itself)

// 1D tests: H(grad) and H(vol) bases defined
// p1, p1:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D1_P1_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D1_P1_P1, HGRAD, VALUE, HGRAD, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D1_P1_P1, HVOL,  VALUE, HGRAD, VALUE)
// p1, p2:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D1_P1_P2, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D1_P1_P2, HGRAD, VALUE, HGRAD, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D1_P1_P2, HVOL,  VALUE, HGRAD, VALUE)
// p2, p1:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D1_P2_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D1_P2_P1, HGRAD, VALUE, HGRAD, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D1_P2_P1, HVOL,  VALUE, HGRAD, VALUE)

// 1D vector-weighted test
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardVectorWeighted_D1_P1_P1, HGRAD, GRAD, HGRAD, GRAD)

// 1D scalar against vector-weighted tests
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardScalarAgainstVectorDotVector_D1_P1_P1, HVOL, VALUE, HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardScalarAgainstVectorDotVector_D1_P1_P1, HGRAD, VALUE, HGRAD, GRAD)

// 1D vector-weighted against scalar tests
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardVectorDotVectorAgainstScalar_D1_P1_P1, HGRAD, GRAD, HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardVectorDotVectorAgainstScalar_D1_P1_P1, HGRAD, GRAD, HGRAD, VALUE)

// 2D tests: curls of H(curl) are scalars.
// p1, p1:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P1_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P1_P1, HGRAD, GRAD,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P1_P1, HGRAD, GRAD,  HCURL, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P1_P1, HGRAD, VALUE, HDIV,  DIV)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P1_P1, HGRAD, VALUE, HCURL, CURL)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P1_P1, HDIV,  DIV,   HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P1_P1, HCURL, CURL,  HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P1_P1, HVOL,  VALUE, HGRAD, VALUE)
// p2, p1:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P2_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P2_P1, HGRAD, GRAD,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P2_P1, HGRAD, GRAD,  HCURL, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P2_P1, HGRAD, VALUE, HDIV,  DIV)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P2_P1, HGRAD, VALUE, HCURL, CURL)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P2_P1, HDIV,  DIV,   HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P2_P1, HCURL, CURL,  HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P2_P1, HVOL,  VALUE, HGRAD, VALUE)
// p1, p2:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P1_P2, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P1_P2, HGRAD, GRAD,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P1_P2, HGRAD, GRAD,  HCURL, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P1_P2, HGRAD, VALUE, HDIV,  DIV)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P1_P2, HGRAD, VALUE, HCURL, CURL)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P1_P2, HDIV,  DIV,   HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P1_P2, HCURL, CURL,  HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D2_P1_P2, HVOL,  VALUE, HGRAD, VALUE)

// 2D vector-weighted tests
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardVectorWeighted_D2_P1_P1, HGRAD, GRAD, HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardVectorWeighted_D2_P2_P1, HGRAD, GRAD, HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardVectorWeighted_D2_P1_P1, HCURL, VALUE, HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardVectorWeighted_D2_P2_P1, HCURL, VALUE, HDIV,  VALUE)

// 2D scalar against vector-weighted tests
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardScalarAgainstVectorDotVector_D2_P1_P1, HVOL, VALUE, HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardScalarAgainstVectorDotVector_D2_P1_P1, HGRAD, VALUE, HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardScalarAgainstVectorDotVector_D2_P1_P1, HGRAD, VALUE, HDIV, VALUE)

// 2D vector-weighted against scalar tests
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardVectorDotVectorAgainstScalar_D2_P1_P1, HGRAD, GRAD, HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardVectorDotVectorAgainstScalar_D2_P1_P1, HGRAD, GRAD, HGRAD, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardVectorDotVectorAgainstScalar_D2_P1_P1, HDIV, VALUE, HGRAD, VALUE)

// 3D tests: curls of H(curl) are vectors
// p1, p1:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P1_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P1_P1, HGRAD, GRAD,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P1_P1, HGRAD, GRAD,  HCURL, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P1_P1, HGRAD, GRAD,  HCURL, CURL)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P1_P1, HGRAD, VALUE, HDIV,  DIV)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P1_P1, HDIV,  DIV,   HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P1_P1, HCURL, CURL,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P1_P1, HVOL,  VALUE, HGRAD, VALUE)
// p2, p1:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P2_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P2_P1, HGRAD, GRAD,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P2_P1, HGRAD, GRAD,  HCURL, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P2_P1, HGRAD, GRAD,  HCURL, CURL)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P2_P1, HGRAD, VALUE, HDIV,  DIV)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P2_P1, HDIV,  DIV,   HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P2_P1, HCURL, CURL,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P2_P1, HVOL,  VALUE, HGRAD, VALUE)
// p1, p2:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P1_P2, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P1_P2, HGRAD, GRAD,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P1_P2, HGRAD, GRAD,  HCURL, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P1_P2, HGRAD, GRAD,  HCURL, CURL)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P1_P2, HGRAD, VALUE, HDIV,  DIV)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P1_P2, HDIV,  DIV,   HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P1_P2, HCURL, CURL,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandard_D3_P1_P2, HVOL,  VALUE, HGRAD, VALUE)

// 3D vector-weighted tests
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardVectorWeighted_D3_P2_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardVectorWeighted_D3_P2_P1, HCURL, VALUE, HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardVectorWeighted_D3_P2_P1, HCURL, CURL,  HGRAD, GRAD)

// 3D scalar against vector-weighted tests
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardScalarAgainstVectorDotVector_D3_P1_P1, HVOL, VALUE, HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardScalarAgainstVectorDotVector_D3_P1_P1, HGRAD, VALUE, HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardScalarAgainstVectorDotVector_D3_P1_P1, HGRAD, VALUE, HDIV, VALUE)

// 3D vector-weighted against scalar tests
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardVectorDotVectorAgainstScalar_D3_P1_P1, HGRAD, GRAD, HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardVectorDotVectorAgainstScalar_D3_P1_P1, HGRAD, GRAD, HGRAD, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, PAMatrixVersusStandardVectorDotVectorAgainstScalar_D3_P1_P1, HDIV, VALUE, HGRAD, VALUE)

} // anonymous namespace
