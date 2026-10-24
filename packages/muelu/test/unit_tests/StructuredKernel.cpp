// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include <Teuchos_ScalarTraits.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "MueLu_config.hpp"

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Galeri_XpetraMaps.hpp>

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_TripleMatrixMultiply.hpp>

#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_GeometricInterpolationPFactory.hpp"
#include "MueLu_NoFactory.hpp"
#include "MueLu_StructuredAggregationFactory.hpp"
#include "MueLu_StructuredRAPFactory.hpp"
#include "MueLu_StructuredRAPKernel.hpp"

namespace MueLuTests {

namespace {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
struct StructuredKernelData {
  using Matrix = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  Teuchos::RCP<Matrix> A;
  Teuchos::RCP<Matrix> P;
  Teuchos::Array<LocalOrdinal> fineNodesPerDim;
  Teuchos::Array<LocalOrdinal> coarseNodesPerDim;
  int numDimensions;
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
StructuredKernelData<Scalar, LocalOrdinal, GlobalOrdinal, Node>
buildStructuredKernelData(const std::string& matrixType,
                          const GlobalOrdinal nx,
                          const GlobalOrdinal ny,
                          const GlobalOrdinal nz,
                          const int coarseningRate) {
  using SC                    = Scalar;
  using LO                    = LocalOrdinal;
  using GO                    = GlobalOrdinal;
  using NO                    = Node;
  using Map                   = Xpetra::Map<LO, GO, NO>;
  using Matrix                = Xpetra::Matrix<SC, LO, GO, NO>;
  using CrsMatrixWrap         = Xpetra::CrsMatrixWrap<SC, LO, GO, NO>;
  using MultiVector           = Xpetra::MultiVector<SC, LO, GO, NO>;
  using coordinate_type       = typename Teuchos::ScalarTraits<SC>::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<coordinate_type, LO, GO, NO>;
  using AmalgamationFactory   = MueLu::AmalgamationFactory<SC, LO, GO, NO>;
  using CoalesceDropFactory   = MueLu::CoalesceDropFactory<SC, LO, GO, NO>;
  using AggregationFactory    = MueLu::StructuredAggregationFactory<SC, LO, GO, NO>;
  using ProlongatorFactory    = MueLu::GeometricInterpolationPFactory<SC, LO, GO, NO>;

  Teuchos::RCP<const Teuchos::Comm<int> > comm =
      TestHelpers::Parameters::getDefaultComm();

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  if (ny > 0)
    galeriList.set("ny", ny);
  if (nz > 0)
    galeriList.set("nz", nz);

  StructuredKernelData<SC, LO, GO, NO> data;
  data.fineNodesPerDim = Teuchos::Array<LO>(3, Teuchos::as<LO>(1));

  std::string mapType;
  std::string coordinateType;
  int dofsPerNode = 1;
  if (matrixType == "Laplace1D") {
    data.numDimensions = 1;
    mapType            = "Cartesian1D";
    coordinateType     = "1D";
  } else if (matrixType == "Laplace2D" || matrixType == "Elasticity2D") {
    data.numDimensions = 2;
    mapType            = "Cartesian2D";
    coordinateType     = "2D";
    dofsPerNode        = matrixType == "Elasticity2D" ? 2 : 1;
  } else if (matrixType == "Laplace3D" || matrixType == "Elasticity3D") {
    data.numDimensions = 3;
    mapType            = "Cartesian3D";
    coordinateType     = "3D";
    dofsPerNode        = matrixType == "Elasticity3D" ? 3 : 1;
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::runtime_error,
        "Unsupported StructuredRAPKernel test matrix type: " << matrixType);
  }

  Teuchos::RCP<const Map> nodalMap =
      Galeri::Xpetra::CreateMap<LO, GO, NO>(
          TestHelpers::Parameters::getLib(), mapType, comm, galeriList);
  Teuchos::RCP<RealValuedMultiVector> coordinates =
      Galeri::Xpetra::Utils::CreateCartesianCoordinates<
          coordinate_type, LO, GO, Map, RealValuedMultiVector>(
          coordinateType, nodalMap, galeriList);
  Teuchos::RCP<const Map> dofMap = nodalMap;
  if (dofsPerNode > 1)
    dofMap = Xpetra::MapFactory<LO, GO, NO>::Build(nodalMap, dofsPerNode);

  Teuchos::RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> >
      galeriProblem =
          Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap,
                                       MultiVector>(
              matrixType, dofMap, galeriList);

  data.A = galeriProblem->BuildMatrix();
  data.A->SetFixedBlockSize(dofsPerNode);
  data.fineNodesPerDim[0] = galeriList.get<LO>("lnx");
  if (data.numDimensions > 1)
    data.fineNodesPerDim[1] = galeriList.get<LO>("lny");
  if (data.numDimensions > 2)
    data.fineNodesPerDim[2] = galeriList.get<LO>("lnz");

  MueLu::Level fineLevel, coarseLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(
      fineLevel, coarseLevel);
  fineLevel.Set("A", data.A);
  fineLevel.Set("Coordinates", coordinates);
  fineLevel.Set("numDimensions", data.numDimensions);
  fineLevel.Set("lNodesPerDim", data.fineNodesPerDim);

  Teuchos::RCP<MultiVector> nullspace =
      Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(
          data.A->getRowMap(), 1);
  nullspace->putScalar(Teuchos::ScalarTraits<SC>::one());
  fineLevel.Set("Nullspace", nullspace);

  Teuchos::RCP<AmalgamationFactory> amalgamation =
      Teuchos::rcp(new AmalgamationFactory());
  Teuchos::RCP<CoalesceDropFactory> coalesceDrop =
      Teuchos::rcp(new CoalesceDropFactory());
  coalesceDrop->SetFactory("UnAmalgamationInfo", amalgamation);

  Teuchos::RCP<AggregationFactory> aggregation =
      Teuchos::rcp(new AggregationFactory());
  aggregation->SetParameter(
      "aggregation: mode", Teuchos::ParameterEntry(std::string("uncoupled")));
  aggregation->SetParameter(
      "aggregation: output type",
      Teuchos::ParameterEntry(std::string("CrsGraph")));
  aggregation->SetParameter(
      "aggregation: coarsening order", Teuchos::ParameterEntry(0));
  aggregation->SetParameter(
      "aggregation: coarsening rate",
      Teuchos::ParameterEntry(std::string("{") +
                              std::to_string(coarseningRate) + "}"));
  aggregation->SetFactory("Graph", coalesceDrop);
  aggregation->SetFactory("DofsPerNode", coalesceDrop);

  Teuchos::RCP<ProlongatorFactory> prolongator =
      Teuchos::rcp(new ProlongatorFactory());
  prolongator->SetFactory("A", MueLu::NoFactory::getRCP());
  prolongator->SetFactory("Coordinates", MueLu::NoFactory::getRCP());
  prolongator->SetFactory("Nullspace", MueLu::NoFactory::getRCP());
  prolongator->SetFactory("prolongatorGraph", aggregation);
  prolongator->SetFactory("coarseCoordinatesFineMap", aggregation);
  prolongator->SetFactory("coarseCoordinatesMap", aggregation);
  prolongator->SetFactory("numDimensions", aggregation);
  prolongator->SetFactory("lCoarseNodesPerDim", aggregation);
  prolongator->SetFactory("structuredInterpolationOrder", aggregation);

  coarseLevel.Request("P", prolongator.get());
  coarseLevel.Request(*prolongator);
  prolongator->Build(fineLevel, coarseLevel);

  data.P = coarseLevel.Get<Teuchos::RCP<Matrix> >(
      "P", prolongator.get());
  data.coarseNodesPerDim = fineLevel.Get<Teuchos::Array<LO> >(
      "lCoarseNodesPerDim", aggregation.get());
  return data;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
computeStructuredRAP(
    const StructuredKernelData<Scalar, LocalOrdinal, GlobalOrdinal, Node>& data,
    const int coarseningRate) {
  using SC         = Scalar;
  using LO         = LocalOrdinal;
  using GO         = GlobalOrdinal;
  using NO         = Node;
  using Matrix     = Xpetra::Matrix<SC, LO, GO, NO>;
  using RAPFactory = MueLu::StructuredRAPFactory<SC, LO, GO, NO>;

  RAPFactory rapFactory;
  const typename RAPFactory::FineStencilSpec fineStencil =
      rapFactory.DetectFineStencil(*data.A, data.numDimensions,
                                   data.fineNodesPerDim);
  const typename RAPFactory::StructuredGraphSpec graphSpec =
      rapFactory.DeriveCoarseRAPStencil(fineStencil, 0);

  Teuchos::RCP<Matrix> Ac;
  rapFactory.GetStructuredGraph(Ac, data.P, data.coarseNodesPerDim,
                                graphSpec);

  Teuchos::Array<int> rate(1, coarseningRate);
  Teuchos::RCP<Teuchos::ParameterList> kernelParams =
      Teuchos::rcp(new Teuchos::ParameterList());
  MueLu::Details::StructuredRAPKernel<SC, LO, GO, NO>::Compute(
      *data.A, *data.P, *Ac, fineStencil, 0, data.fineNodesPerDim,
      data.coarseNodesPerDim, rate, kernelParams);
  return Ac;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
computeXpetraRAP(
    const StructuredKernelData<Scalar, LocalOrdinal, GlobalOrdinal, Node>& data) {
  using SC     = Scalar;
  using LO     = LocalOrdinal;
  using GO     = GlobalOrdinal;
  using NO     = Node;
  using Matrix = Xpetra::Matrix<SC, LO, GO, NO>;

  Teuchos::RCP<Matrix> Ac = Xpetra::MatrixFactory<SC, LO, GO, NO>::Build(
      data.P->getDomainMap(), Teuchos::as<LO>(0));
  Teuchos::RCP<Teuchos::ParameterList> params =
      Teuchos::rcp(new Teuchos::ParameterList());
  Xpetra::TripleMatrixMultiply<SC, LO, GO, NO>::MultiplyRAP(
      *data.P, true, *data.A, false, *data.P, false, *Ac, true, true,
      "MueLu StructuredRAPKernel unit-test reference", params);
  return Ac;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
relativeFrobeniusError(
    const Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& actual,
    const Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& expected,
    Teuchos::FancyOStream& out) {
  using SC        = Scalar;
  using LO        = LocalOrdinal;
  using GO        = GlobalOrdinal;
  using NO        = Node;
  using Matrix    = Xpetra::Matrix<SC, LO, GO, NO>;
  using magnitude = typename Teuchos::ScalarTraits<SC>::magnitudeType;

  Teuchos::RCP<Matrix> difference;
  Xpetra::MatrixMatrix<SC, LO, GO, NO>::TwoMatrixAdd(
      *actual, false, Teuchos::ScalarTraits<SC>::one(),
      *expected, false, -Teuchos::ScalarTraits<SC>::one(),
      difference, out);
  if (!difference->isFillComplete())
    difference->fillComplete(expected->getDomainMap(),
                             expected->getRangeMap());

  const magnitude scale = std::max(
      expected->getFrobeniusNorm(), Teuchos::ScalarTraits<magnitude>::one());
  return difference->getFrobeniusNorm() / scale;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void compareConstantStructuredProblem(const std::string& matrixType,
                                      const GlobalOrdinal nx,
                                      const GlobalOrdinal ny,
                                      const GlobalOrdinal nz,
                                      const int coarseningRate,
                                      Teuchos::FancyOStream& out,
                                      bool& success) {
  using SC        = Scalar;
  using LO        = LocalOrdinal;
  using GO        = GlobalOrdinal;
  using NO        = Node;
  using magnitude = typename Teuchos::ScalarTraits<SC>::magnitudeType;

  const StructuredKernelData<SC, LO, GO, NO> data =
      buildStructuredKernelData<SC, LO, GO, NO>(
          matrixType, nx, ny, nz, coarseningRate);
  const Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO> > structuredAc =
      computeStructuredRAP<SC, LO, GO, NO>(data, coarseningRate);
  const Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO> > xpetraAc =
      computeXpetraRAP<SC, LO, GO, NO>(data);

  const magnitude error = relativeFrobeniusError<SC, LO, GO, NO>(
      structuredAc, xpetraAc, out);
  const magnitude tolerance =
      magnitude(1000) * Teuchos::ScalarTraits<magnitude>::eps();
  out << matrixType << " " << nx;
  if (data.numDimensions > 1)
    out << "x" << ny;
  if (data.numDimensions > 2)
    out << "x" << nz;
  out
      << ", coarsening rate " << coarseningRate
      << ": relative Frobenius error = " << error << std::endl;
  TEST_COMPARE(error, <=, tolerance);
}

}  // namespace

template <class Node>
bool canRunStructuredKernelTest(Teuchos::FancyOStream& out) {
  constexpr bool hostMemory =
      std::is_same<typename Node::device_type::memory_space,
                   Kokkos::HostSpace>::value;
  const bool canRun = hostMemory &&
                      TestHelpers::Parameters::getLib() == Xpetra::UseTpetra;
  if (!canRun)
    out << "The direct StructuredRAPKernel test requires host-memory Tpetra; "
           "skipping this test."
        << std::endl;
  return canRun;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
    StructuredKernel, ConstantLaplace1DMatchesXpetra,
    Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  if (!canRunStructuredKernelTest<NO>(out))
    return;

  compareConstantStructuredProblem<SC, LO, GO, NO>(
      "Laplace1D", 13, -1, -1, 3, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
    StructuredKernel, ConstantLaplace2DMatchesXpetra,
    Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  if (!canRunStructuredKernelTest<NO>(out))
    return;

  compareConstantStructuredProblem<SC, LO, GO, NO>(
      "Laplace2D", 13, 10, -1, 3, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
    StructuredKernel, ConstantLaplace3DMatchesXpetra,
    Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  if (!canRunStructuredKernelTest<NO>(out))
    return;

  compareConstantStructuredProblem<SC, LO, GO, NO>(
      "Laplace3D", 8, 7, 6, 2, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
    StructuredKernel, ConstantElasticity2DMatchesXpetra,
    Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  if (!canRunStructuredKernelTest<NO>(out))
    return;

  compareConstantStructuredProblem<SC, LO, GO, NO>(
      "Elasticity2D", 12, 12, -1, 3, out, success);
  compareConstantStructuredProblem<SC, LO, GO, NO>(
      "Elasticity2D", 13, 10, -1, 3, out, success);
  compareConstantStructuredProblem<SC, LO, GO, NO>(
      "Elasticity2D", 11, 12, -1, 2, out, success);
  compareConstantStructuredProblem<SC, LO, GO, NO>(
      "Elasticity2D", 14, 11, -1, 4, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
    StructuredKernel, ConstantElasticity3DMatchesXpetra,
    Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  if (!canRunStructuredKernelTest<NO>(out))
    return;

  compareConstantStructuredProblem<SC, LO, GO, NO>(
      "Elasticity3D", 7, 6, 5, 2, out, success);
}

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(                    \
      StructuredKernel, ConstantLaplace1DMatchesXpetra,    \
      Scalar, LO, GO, Node)                                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(                    \
      StructuredKernel, ConstantLaplace2DMatchesXpetra,    \
      Scalar, LO, GO, Node)                                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(                    \
      StructuredKernel, ConstantLaplace3DMatchesXpetra,    \
      Scalar, LO, GO, Node)                                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(                    \
      StructuredKernel, ConstantElasticity2DMatchesXpetra, \
      Scalar, LO, GO, Node)                                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(                    \
      StructuredKernel, ConstantElasticity3DMatchesXpetra, \
      Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
