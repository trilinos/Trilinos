// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   PAMatrixTests_Serendipity.cpp
    \brief  Tests against PAMatrix::assemble() for Serendipity bases (which PAMatrix treats as ordinal-filtered tensor bases),
            comparing with standard assembly.
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
#include <Intrepid2_HierarchicalBasisFamily.hpp>
#include <Intrepid2_NodalBasisFamily.hpp>
#include <Intrepid2_SerendipityBasisFamily.hpp>
#include <Intrepid2_TestUtils.hpp>

#include <Intrepid2_CellGeometry.hpp>
#include "Intrepid2_Data.hpp"
#include "Intrepid2_TensorData.hpp"
#include "Intrepid2_TensorPoints.hpp"
#include "Intrepid2_TransformedBasisValues.hpp"
#include "Intrepid2_VectorData.hpp"
#include "Intrepid2_ScalarView.hpp"

#include "PAMatrixAssembly.hpp"

#include "StructuredIntegrationTests_TagDefs.hpp"
#include "StructuredIntegrationTests_Utils.hpp"

namespace
{
  using namespace Intrepid2;

  //! Serendipity field ordinal --> full tensor basis field ordinal, on host.
  template<class SerendipityFamily>
  std::vector<int> serendipityOrdinals(const shards::CellTopology &cellTopo, const EFunctionSpace &fs, const int &p)
  {
    using DeviceType = typename SerendipityFamily::DeviceType;
    using BasisBase  = typename SerendipityFamily::Basis;
    auto basis = getBasis< SerendipityFamily >(cellTopo, fs, p);
    auto serendipityBasis = Teuchos::rcp_dynamic_cast< SerendipityBasis<BasisBase> >(basis, true);
    auto ordinalMapHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), serendipityBasis->ordinalMap());
    std::vector<int> ordinals(ordinalMapHost.extent(0));
    for (unsigned i=0; i<ordinals.size(); i++) ordinals[i] = ordinalMapHost(i);
    return ordinals;
  }

  //! Compares PAMatrix::assemble() for Serendipity bases with the rows and columns of the full-tensor matrix that belong to
  //! the Serendipity subspaces.  The full-tensor reference matrix is computed by referenceAssembly(), which uses the generic
  //! IntegrationTools::integrate() path.  (SerendipityBasis does not support the getValues() variant used by standard assembly.)
  //! PAMatrix does not support ordinal filters together with orientations, so neither matrix applies orientations.
  template<class Scalar, class PointScalar, int spaceDim, typename DeviceType>
  void testFilteredFullTensorVersusPAMatrixAssemble2(const int &meshWidth,
                                                     const EFunctionSpace &fs1, const EOperator &op1, const int &p1,
                                                     Teuchos::RCP< Kokkos::Array<PointScalar,spaceDim> > vectorWeight1,
                                                     const EFunctionSpace &fs2, const EOperator &op2, const int &p2,
                                                     Teuchos::RCP< Kokkos::Array<PointScalar,spaceDim> > vectorWeight2,
                                                     const double &relTol, const double &absTol, Teuchos::FancyOStream &out, bool &success)
  {
    using FullFamily        = HierarchicalBasisFamily<DeviceType>;
    using SerendipityFamily = SerendipityBasisFamily<DeviceType>;

    Kokkos::Array<int,spaceDim> gridCellCounts;
    for (int d=0; d<spaceDim; d++) gridCellCounts[d] = meshWidth;
    auto geometry = getMesh<PointScalar, spaceDim, DeviceType>(Standard, gridCellCounts);

    const bool applyOrientations = false;
    auto fullPAMatrix = constructPAMatrix<Scalar,FullFamily>(geometry, p1, fs1, op1, vectorWeight1, p2, fs2, op2, vectorWeight2, applyOrientations);
    auto fullIntegrals = fullPAMatrix.allocateMatrixStorage();
    referenceAssembly(fullPAMatrix, fullIntegrals);

    const auto cellTopo = geometry.cellTopology();
    const auto leftOrdinals  = serendipityOrdinals<SerendipityFamily>(cellTopo, fs1, p1);
    const auto rightOrdinals = serendipityOrdinals<SerendipityFamily>(cellTopo, fs2, p2);

    const int numCells = geometry.numCells();
    using HostView = Kokkos::View<Scalar***,Kokkos::HostSpace>;
    HostView standardIntegrals("filtered full-tensor integrals", numCells, leftOrdinals.size(), rightOrdinals.size());
    for (int cell=0; cell<numCells; cell++)
      for (unsigned i=0; i<leftOrdinals.size(); i++)
        for (unsigned j=0; j<rightOrdinals.size(); j++)
          standardIntegrals(cell,i,j) = fullIntegrals(cell,leftOrdinals[i],rightOrdinals[j]); // Data access is host-accessible on the (serial) test device

    auto paMatrix = constructPAMatrix<Scalar,SerendipityFamily>(geometry, p1, fs1, op1, vectorWeight1, p2, fs2, op2, vectorWeight2, applyOrientations);
    auto paIntegrals = paMatrix.allocateMatrixStorage();
    paMatrix.assemble(paIntegrals);
    auto paIntegralsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), paIntegrals.getUnderlyingView());

    out << "Comparing filtered full-tensor assembly to PAMatrix::assemble() for Serendipity bases…\n";
    testFloatingEquality3(standardIntegrals, paIntegralsHost, relTol, absTol, out, success, "filtered full-tensor integral", "PAMatrix assemble");
  }

  template<int spaceDim>
  Teuchos::RCP< Kokkos::Array<double,spaceDim> > testVectorWeight(const double firstWeight, const double factor)
  {
    auto weight = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
    double w = firstWeight;
    for (int d=0; d<spaceDim; d++)
    {
      (*weight)[d] = w;
      w *= factor;
    }
    return weight;
  }

  template<int spaceDim, class FS1Tag, class Op1Tag, class FS2Tag, class Op2Tag>
  void runSerendipityCase(const int p1, const int p2, const bool vectorWeighted, Teuchos::FancyOStream &out, bool &success)
  {
    using DeviceType  = DefaultTestDeviceType;
    const int meshWidth = 2;
    // copy the tag values (the tags' static members have in-class initializers but no out-of-line definitions)
    const EFunctionSpace fs1 = FS1Tag::functionSpace;
    const EFunctionSpace fs2 = FS2Tag::functionSpace;
    const EOperator op1 = Op1Tag::op;
    const EOperator op2 = Op2Tag::op;
    Teuchos::RCP< Kokkos::Array<double,spaceDim> > vectorWeight1, vectorWeight2;
    // a (value, gradient) pairing requires a vector weight on the gradient; when vectorWeighted is set, weight all gradients
    const bool mixedRanks = (op1 == OPERATOR_GRAD) != (op2 == OPERATOR_GRAD);
    if (vectorWeighted || mixedRanks)
    {
      if (op1 == OPERATOR_GRAD) vectorWeight1 = testVectorWeight<spaceDim>(1.0, 0.5);
      if (op2 == OPERATOR_GRAD) vectorWeight2 = testVectorWeight<spaceDim>(0.5, 2.0);
    }
    const double relTol = 1e-11;
    const double absTol = 1e-11;
    testFilteredFullTensorVersusPAMatrixAssemble2<double, double, spaceDim, DeviceType>
      (meshWidth, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2,
       relTol, absTol, out, success);
  }

  // MARK: Serendipity, 2D
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, Serendipity_D2_P2_P2, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
  {
    runSerendipityCase<2,FS1Tag,Op1Tag,FS2Tag,Op2Tag>(2, 2, false, out, success);
  }
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, Serendipity_D2_P3_P3, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
  {
    runSerendipityCase<2,FS1Tag,Op1Tag,FS2Tag,Op2Tag>(3, 3, false, out, success);
  }
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, Serendipity_D2_P2_P3, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
  {
    runSerendipityCase<2,FS1Tag,Op1Tag,FS2Tag,Op2Tag>(2, 3, false, out, success);
  }
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, SerendipityVectorWeighted_D2_P3_P3, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
  {
    runSerendipityCase<2,FS1Tag,Op1Tag,FS2Tag,Op2Tag>(3, 3, true, out, success);
  }

  // MARK: Serendipity, 3D
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, Serendipity_D3_P2_P2, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
  {
    runSerendipityCase<3,FS1Tag,Op1Tag,FS2Tag,Op2Tag>(2, 2, false, out, success);
  }
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, Serendipity_D3_P3_P3, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
  {
    runSerendipityCase<3,FS1Tag,Op1Tag,FS2Tag,Op2Tag>(3, 3, false, out, success);
  }
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, SerendipityVectorWeighted_D3_P3_P3, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
  {
    runSerendipityCase<3,FS1Tag,Op1Tag,FS2Tag,Op2Tag>(3, 3, true, out, success);
  }

#define SERENDIPITY_OPERATOR_PAIRS(TESTNAME) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, TESTNAME, HGRAD, VALUE, HGRAD, VALUE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, TESTNAME, HGRAD, GRAD,  HGRAD, GRAD)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, TESTNAME, HGRAD, VALUE, HGRAD, GRAD)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, TESTNAME, HGRAD, GRAD,  HGRAD, VALUE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, TESTNAME, HVOL,  VALUE, HVOL,  VALUE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, TESTNAME, HVOL,  VALUE, HGRAD, GRAD)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, TESTNAME, HGRAD, GRAD,  HVOL,  VALUE)

  SERENDIPITY_OPERATOR_PAIRS(Serendipity_D2_P2_P2)
  SERENDIPITY_OPERATOR_PAIRS(Serendipity_D2_P3_P3)
  SERENDIPITY_OPERATOR_PAIRS(Serendipity_D2_P2_P3)
  SERENDIPITY_OPERATOR_PAIRS(SerendipityVectorWeighted_D2_P3_P3)
  SERENDIPITY_OPERATOR_PAIRS(Serendipity_D3_P2_P2)
  SERENDIPITY_OPERATOR_PAIRS(Serendipity_D3_P3_P3)
  SERENDIPITY_OPERATOR_PAIRS(SerendipityVectorWeighted_D3_P3_P3)

} // anonymous namespace
