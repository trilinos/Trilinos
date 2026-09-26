// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <algorithm>
#include <array>
#include <sstream>
#include <tuple>
#include <vector>

#include "MueLu_config.hpp"

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Galeri_XpetraMaps.hpp>

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_GeometricInterpolationPFactory.hpp"
#include "MueLu_NoFactory.hpp"
#include "MueLu_StructuredAggregationFactory.hpp"
#include "MueLu_StructuredRAPFactory.hpp"

namespace MueLuTests {

namespace {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
struct StructuredProblemData {
  using Matrix                = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using RealValuedMultiVector = Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node>;
  Teuchos::RCP<Matrix> A;
  Teuchos::RCP<RealValuedMultiVector> coordinates;
  Teuchos::Array<LocalOrdinal> lNodesPerDim;
  int numDimensions;
  int dofsPerNode;
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
StructuredProblemData<Scalar, LocalOrdinal, GlobalOrdinal, Node>
buildStructuredProblem(const std::string& matrixType,
                       const GlobalOrdinal nx,
                       const GlobalOrdinal ny,
                       const GlobalOrdinal nz,
                       const GlobalOrdinal mx,
                       const GlobalOrdinal my,
                       const GlobalOrdinal mz) {
  using SC                    = Scalar;
  using LO                    = LocalOrdinal;
  using GO                    = GlobalOrdinal;
  using NO                    = Node;
  using Map                   = Xpetra::Map<LO, GO, NO>;
  using Matrix                = Xpetra::Matrix<SC, LO, GO, NO>;
  using CrsMatrixWrap         = Xpetra::CrsMatrixWrap<SC, LO, GO, NO>;
  using MultiVector           = Xpetra::MultiVector<SC, LO, GO, NO>;
  using RealValuedMultiVector = typename StructuredProblemData<SC, LO, GO, NO>::RealValuedMultiVector;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  if (ny > 0)
    galeriList.set("ny", ny);
  if (nz > 0)
    galeriList.set("nz", nz);
  if (mx > 0)
    galeriList.set("mx", mx);
  if (my > 0)
    galeriList.set("my", my);
  if (mz > 0)
    galeriList.set("mz", mz);

  StructuredProblemData<SC, LO, GO, NO> problem;
  problem.lNodesPerDim = Teuchos::Array<LO>(3, Teuchos::as<LO>(1));
  problem.dofsPerNode  = 1;

  Teuchos::RCP<const Map> map;
  if (matrixType == "Laplace1D") {
    problem.numDimensions   = 1;
    map                     = Galeri::Xpetra::CreateMap<LO, GO, NO>(TestHelpers::Parameters::getLib(), "Cartesian1D", comm, galeriList);
    problem.coordinates     = Galeri::Xpetra::Utils::CreateCartesianCoordinates<typename RealValuedMultiVector::scalar_type, LO, GO, Map, RealValuedMultiVector>("1D", map, galeriList);
    problem.lNodesPerDim[0] = galeriList.get<LO>("lnx");
  } else if (matrixType == "Laplace2D" || matrixType == "Elasticity2D") {
    problem.numDimensions   = 2;
    map                     = Galeri::Xpetra::CreateMap<LO, GO, NO>(TestHelpers::Parameters::getLib(), "Cartesian2D", comm, galeriList);
    problem.coordinates     = Galeri::Xpetra::Utils::CreateCartesianCoordinates<typename RealValuedMultiVector::scalar_type, LO, GO, Map, RealValuedMultiVector>("2D", map, galeriList);
    problem.lNodesPerDim[0] = galeriList.get<LO>("lnx");
    problem.lNodesPerDim[1] = galeriList.get<LO>("lny");
  } else if (matrixType == "Laplace3D" || matrixType == "Elasticity3D") {
    problem.numDimensions   = 3;
    map                     = Galeri::Xpetra::CreateMap<LO, GO, NO>(TestHelpers::Parameters::getLib(), "Cartesian3D", comm, galeriList);
    problem.coordinates     = Galeri::Xpetra::Utils::CreateCartesianCoordinates<typename RealValuedMultiVector::scalar_type, LO, GO, Map, RealValuedMultiVector>("3D", map, galeriList);
    problem.lNodesPerDim[0] = galeriList.get<LO>("lnx");
    problem.lNodesPerDim[1] = galeriList.get<LO>("lny");
    problem.lNodesPerDim[2] = galeriList.get<LO>("lnz");
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Unsupported StructuredRAPFactory test matrix type: " << matrixType);
  }

  if (matrixType == "Elasticity2D") {
    problem.dofsPerNode = 2;
    map                 = Xpetra::MapFactory<LO, GO, NO>::Build(map, problem.dofsPerNode);
  } else if (matrixType == "Elasticity3D") {
    problem.dofsPerNode = 3;
    map                 = Xpetra::MapFactory<LO, GO, NO>::Build(map, problem.dofsPerNode);
  }

  Teuchos::RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > galeriProblem =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixType, map, galeriList);
  problem.A = galeriProblem->BuildMatrix();
  if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D")
    problem.A->SetFixedBlockSize(problem.dofsPerNode);

  return problem;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void checkDetectedFineStencil(
    const StructuredProblemData<Scalar, LocalOrdinal, GlobalOrdinal, Node>& problem,
    const std::string& matrixType) {
  using LO         = LocalOrdinal;
  using RAPFactory = MueLu::StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Offset     = std::array<int, 3>;
  using Coupling   = std::tuple<int, int, int, LO, LO>;

  RAPFactory rap;
  const typename RAPFactory::FineStencilSpec detected =
      rap.DetectFineStencil(*problem.A, problem.numDimensions, problem.lNodesPerDim);

  TEUCHOS_TEST_FOR_EXCEPTION(
      detected.numDimensions != problem.numDimensions, std::runtime_error,
      matrixType << ": detected " << detected.numDimensions
                 << " dimensions, expected " << problem.numDimensions << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(
      detected.dofsPerNode != problem.dofsPerNode, std::runtime_error,
      matrixType << ": detected " << detected.dofsPerNode
                 << " DOFs per node, expected " << problem.dofsPerNode << ".");

  std::vector<Offset> expectedOffsets;
  if (matrixType == "Laplace1D" || matrixType == "Laplace2D" ||
      matrixType == "Laplace3D") {
    expectedOffsets.push_back(Offset{{0, 0, 0}});
    for (int dim = 0; dim < problem.numDimensions; ++dim) {
      Offset lower{{0, 0, 0}};
      Offset upper{{0, 0, 0}};
      lower[dim] = -1;
      upper[dim] = 1;
      expectedOffsets.push_back(lower);
      expectedOffsets.push_back(upper);
    }
  } else if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D") {
    const int minZ = problem.numDimensions == 3 ? -1 : 0;
    const int maxZ = problem.numDimensions == 3 ? 1 : 0;
    for (int z = minZ; z <= maxZ; ++z)
      for (int y = -1; y <= 1; ++y)
        for (int x = -1; x <= 1; ++x)
          expectedOffsets.push_back(Offset{{x, y, z}});
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::runtime_error,
        "No expected fine stencil is defined for " << matrixType << ".");
  }

  std::vector<Offset> detectedOffsets;
  detectedOffsets.reserve(detected.stencilOffsets.size());
  for (const typename RAPFactory::StencilOffset& offset : detected.stencilOffsets)
    detectedOffsets.push_back(Offset{{offset.x, offset.y, offset.z}});
  std::sort(expectedOffsets.begin(), expectedOffsets.end());
  std::sort(detectedOffsets.begin(), detectedOffsets.end());

  TEUCHOS_TEST_FOR_EXCEPTION(
      detectedOffsets != expectedOffsets, std::runtime_error,
      matrixType << ": detected nodal stencil does not match the expected stencil."
                 << " Detected " << detectedOffsets.size() << " offsets, expected "
                 << expectedOffsets.size() << ".");

  std::vector<Coupling> expectedCouplings;
  expectedCouplings.reserve(expectedOffsets.size() * problem.dofsPerNode * problem.dofsPerNode);
  for (const Offset& offset : expectedOffsets)
    for (LO rowDof = 0; rowDof < Teuchos::as<LO>(problem.dofsPerNode); ++rowDof)
      for (LO columnDof = 0; columnDof < Teuchos::as<LO>(problem.dofsPerNode); ++columnDof)
        expectedCouplings.emplace_back(offset[0], offset[1], offset[2], rowDof, columnDof);

  std::vector<Coupling> detectedCouplings;
  detectedCouplings.reserve(detected.entries.size());
  for (const typename RAPFactory::FineStencilEntry& entry : detected.entries)
    detectedCouplings.emplace_back(entry.offset.x, entry.offset.y, entry.offset.z,
                                   entry.rowDof, entry.columnDof);
  std::sort(expectedCouplings.begin(), expectedCouplings.end());
  std::sort(detectedCouplings.begin(), detectedCouplings.end());

  TEUCHOS_TEST_FOR_EXCEPTION(
      detectedCouplings != expectedCouplings, std::runtime_error,
      matrixType << ": detected scalar stencil couplings do not match the expected stencil."
                 << " Detected " << detectedCouplings.size() << " couplings, expected "
                 << expectedCouplings.size() << ".");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
struct StructuredTransferData {
  using Matrix = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  Teuchos::RCP<Matrix> P;
  Teuchos::Array<LocalOrdinal> lFineNodesPerDim;
  Teuchos::Array<LocalOrdinal> lCoarseNodesPerDim;
  int numDimensions;
  int interpolationOrder;
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
StructuredTransferData<Scalar, LocalOrdinal, GlobalOrdinal, Node>
buildStructuredTransferData(const StructuredProblemData<Scalar, LocalOrdinal, GlobalOrdinal, Node>& problem,
                            const int interpolationOrder,
                            const std::string& coarseningRate) {
  using SC                  = Scalar;
  using LO                  = LocalOrdinal;
  using GO                  = GlobalOrdinal;
  using NO                  = Node;
  using Matrix              = Xpetra::Matrix<SC, LO, GO, NO>;
  using MultiVector         = Xpetra::MultiVector<SC, LO, GO, NO>;
  using AmalgamationFactory = MueLu::AmalgamationFactory<SC, LO, GO, NO>;
  using CoalesceDropFactory = MueLu::CoalesceDropFactory<SC, LO, GO, NO>;
  using AggregationFactory  = MueLu::StructuredAggregationFactory<SC, LO, GO, NO>;
  using ProlongatorFactory  = MueLu::GeometricInterpolationPFactory<SC, LO, GO, NO>;

  MueLu::Level fineLevel, coarseLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

  fineLevel.Set("A", problem.A);
  fineLevel.Set("Coordinates", problem.coordinates);
  fineLevel.Set("numDimensions", problem.numDimensions);
  fineLevel.Set("lNodesPerDim", problem.lNodesPerDim);

  Teuchos::RCP<MultiVector> nullspace =
      Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(problem.A->getRowMap(), 1);
  nullspace->putScalar(Teuchos::ScalarTraits<SC>::one());
  fineLevel.Set("Nullspace", nullspace);

  Teuchos::RCP<AmalgamationFactory> amalgamation = Teuchos::rcp(new AmalgamationFactory());
  Teuchos::RCP<CoalesceDropFactory> coalesceDrop = Teuchos::rcp(new CoalesceDropFactory());
  coalesceDrop->SetFactory("UnAmalgamationInfo", amalgamation);

  Teuchos::RCP<AggregationFactory> aggregation = Teuchos::rcp(new AggregationFactory());
  aggregation->SetParameter("aggregation: mode", Teuchos::ParameterEntry(std::string("uncoupled")));
  aggregation->SetParameter("aggregation: output type", Teuchos::ParameterEntry(std::string("CrsGraph")));
  aggregation->SetParameter("aggregation: coarsening order", Teuchos::ParameterEntry(interpolationOrder));
  aggregation->SetParameter("aggregation: coarsening rate", Teuchos::ParameterEntry(coarseningRate));
  aggregation->SetFactory("Graph", coalesceDrop);
  aggregation->SetFactory("DofsPerNode", coalesceDrop);

  Teuchos::RCP<ProlongatorFactory> prolongator = Teuchos::rcp(new ProlongatorFactory());
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

  StructuredTransferData<SC, LO, GO, NO> transferData;
  transferData.P                = coarseLevel.Get<Teuchos::RCP<Matrix> >("P", prolongator.get());
  transferData.lFineNodesPerDim = problem.lNodesPerDim;
  transferData.lCoarseNodesPerDim =
      fineLevel.Get<Teuchos::Array<LO> >("lCoarseNodesPerDim", aggregation.get());
  transferData.numDimensions      = problem.numDimensions;
  transferData.interpolationOrder = interpolationOrder;
  return transferData;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
buildCoarseMatrix(const StructuredProblemData<Scalar, LocalOrdinal, GlobalOrdinal, Node>& problem,
                  const StructuredTransferData<Scalar, LocalOrdinal, GlobalOrdinal, Node>& transferData,
                  const bool prebuildCoarseGraph) {
  using SC     = Scalar;
  using LO     = LocalOrdinal;
  using GO     = GlobalOrdinal;
  using NO     = Node;
  using Matrix = Xpetra::Matrix<SC, LO, GO, NO>;

  MueLu::Level fineLevel, coarseLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.Set("A", problem.A);
  fineLevel.Set("numDimensions", transferData.numDimensions);
  fineLevel.Set("lNodesPerDim", transferData.lFineNodesPerDim);
  fineLevel.Set("lCoarseNodesPerDim", transferData.lCoarseNodesPerDim);
  fineLevel.Set("structuredInterpolationOrder", transferData.interpolationOrder);
  coarseLevel.Set("P", transferData.P);

  MueLu::StructuredRAPFactory<SC, LO, GO, NO> rap;
  Teuchos::ParameterList rapParams = *rap.GetValidParameterList();
  rapParams.set("rap: triple product", true);
  rapParams.set("rap: prebuild coarse graph", prebuildCoarseGraph);
  rapParams.set("transpose: use implicit", true);
  rap.SetParameterList(rapParams);
  rap.SetFactory("A", MueLu::NoFactory::getRCP());
  rap.SetFactory("P", MueLu::NoFactory::getRCP());
  rap.SetFactory("numDimensions", MueLu::NoFactory::getRCP());
  rap.SetFactory("lNodesPerDim", MueLu::NoFactory::getRCP());
  rap.SetFactory("lCoarseNodesPerDim", MueLu::NoFactory::getRCP());
  rap.SetFactory("structuredInterpolationOrder", MueLu::NoFactory::getRCP());

  coarseLevel.Request("A", &rap);
  coarseLevel.Request(rap);
  rap.Build(fineLevel, coarseLevel);
  return coarseLevel.Get<Teuchos::RCP<Matrix> >("A", &rap);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void compareRAPMatrices(const Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& structuredAc,
                        const Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& referenceAc,
                        Teuchos::FancyOStream& out) {
  using SC        = Scalar;
  using LO        = LocalOrdinal;
  using GO        = GlobalOrdinal;
  using NO        = Node;
  using Map       = Xpetra::Map<LO, GO, NO>;
  using Matrix    = Xpetra::Matrix<SC, LO, GO, NO>;
  using TST       = Teuchos::ScalarTraits<SC>;
  using real_type = typename TST::magnitudeType;

  TEUCHOS_TEST_FOR_EXCEPTION(
      structuredAc.is_null() || referenceAc.is_null(),
      std::runtime_error,
      "Cannot compare RAP matrices because at least one matrix is null.");

  TEUCHOS_TEST_FOR_EXCEPTION(
      !structuredAc->isFillComplete() || !referenceAc->isFillComplete(),
      std::runtime_error,
      "Both RAP matrices must be fill complete before comparison.");

  TEUCHOS_TEST_FOR_EXCEPTION(
      !structuredAc->getRowMap()->isSameAs(*referenceAc->getRowMap()),
      std::runtime_error,
      "RAP matrices have different row maps.");

  TEUCHOS_TEST_FOR_EXCEPTION(
      !structuredAc->getDomainMap()->isSameAs(*referenceAc->getDomainMap()),
      std::runtime_error,
      "RAP matrices have different domain maps.");

  TEUCHOS_TEST_FOR_EXCEPTION(
      !structuredAc->getRangeMap()->isSameAs(*referenceAc->getRangeMap()),
      std::runtime_error,
      "RAP matrices have different range maps.");

  const size_t structuredNnz = structuredAc->getGlobalNumEntries();
  const size_t referenceNnz  = referenceAc->getGlobalNumEntries();
  TEUCHOS_TEST_FOR_EXCEPTION(
      structuredNnz < referenceNnz,
      std::runtime_error,
      "Prebuilt RAP graph has fewer global entries than the symbolic RAP graph: "
          << structuredNnz << " versus " << referenceNnz);

  const double maxGraphOverhead = 0.05;
  const double graphOverhead    = referenceNnz == 0
                                      ? 0.0
                                      : static_cast<double>(structuredNnz - referenceNnz) /
                                         static_cast<double>(referenceNnz);
  TEUCHOS_TEST_FOR_EXCEPTION(
      graphOverhead > maxGraphOverhead,
      std::runtime_error,
      "Prebuilt RAP graph has too many extra entries: "
          << structuredNnz << " versus " << referenceNnz
          << " (overhead = " << graphOverhead
          << ", limit = " << maxGraphOverhead << ")");

  const Teuchos::RCP<const Map> structuredRowMap     = structuredAc->getRowMap();
  const Teuchos::RCP<const Map> structuredColMap     = structuredAc->getColMap();
  const Teuchos::RCP<const Map> referenceColMap      = referenceAc->getColMap();
  const Teuchos::RCP<const Teuchos::Comm<int> > comm = structuredRowMap->getComm();

  int localGraphMismatch = 0;
  std::string localMismatchReason;

  const size_t localNumRows = structuredRowMap->getLocalNumElements();
  for (size_t row = 0; row < localNumRows; ++row) {
    const LO rowLid = Teuchos::as<LO>(row);
    const GO rowGid = structuredRowMap->getGlobalElement(rowLid);

    Teuchos::ArrayView<const LO> structuredIndices;
    Teuchos::ArrayView<const SC> structuredValues;
    Teuchos::ArrayView<const LO> referenceIndices;
    Teuchos::ArrayView<const SC> referenceValues;
    structuredAc->getLocalRowView(rowLid, structuredIndices, structuredValues);
    referenceAc->getLocalRowView(rowLid, referenceIndices, referenceValues);

    std::vector<GO> structuredGids(structuredIndices.size());
    std::vector<GO> referenceGids(referenceIndices.size());
    for (int entry = 0; entry < structuredIndices.size(); ++entry)
      structuredGids[entry] = structuredColMap->getGlobalElement(structuredIndices[entry]);
    for (int entry = 0; entry < referenceIndices.size(); ++entry)
      referenceGids[entry] = referenceColMap->getGlobalElement(referenceIndices[entry]);

    std::sort(structuredGids.begin(), structuredGids.end());
    std::sort(referenceGids.begin(), referenceGids.end());

    const bool structuredHasDuplicates =
        std::adjacent_find(structuredGids.begin(), structuredGids.end()) != structuredGids.end();
    const bool referenceHasDuplicates =
        std::adjacent_find(referenceGids.begin(), referenceGids.end()) != referenceGids.end();

    const bool graphMatches =
        std::includes(structuredGids.begin(), structuredGids.end(),
                      referenceGids.begin(), referenceGids.end());

    if (structuredHasDuplicates || referenceHasDuplicates || !graphMatches) {
      localGraphMismatch = 1;
      std::ostringstream reason;
      reason << "rank " << comm->getRank() << ", row GID " << rowGid
             << ": prebuilt columns = {";
      for (size_t entry = 0; entry < structuredGids.size(); ++entry)
        reason << (entry == 0 ? "" : ", ") << structuredGids[entry];
      reason << "}, symbolic columns = {";
      for (size_t entry = 0; entry < referenceGids.size(); ++entry)
        reason << (entry == 0 ? "" : ", ") << referenceGids[entry];
      reason << "}";
      if (structuredHasDuplicates)
        reason << "; prebuilt row contains duplicate columns";
      if (referenceHasDuplicates)
        reason << "; symbolic row contains duplicate columns";
      localMismatchReason = reason.str();
      break;
    }
  }

  int globalGraphMismatch = 0;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &localGraphMismatch, &globalGraphMismatch);
  const std::string graphMismatchMessage =
      "Prebuilt RAP graph does not contain symbolic RAP graph";
  const std::string graphMismatchDetails =
      localGraphMismatch
          ? graphMismatchMessage + ": " + localMismatchReason
          : graphMismatchMessage + " on another MPI rank.";
  TEUCHOS_TEST_FOR_EXCEPTION(
      globalGraphMismatch != 0,
      std::runtime_error,
      graphMismatchDetails);

  Teuchos::RCP<Matrix> difference;
  Xpetra::MatrixMatrix<SC, LO, GO, NO>::TwoMatrixAdd(
      *structuredAc, false, TST::one(),
      *referenceAc, false, -TST::one(),
      difference, out);

  if (!difference->isFillComplete())
    difference->fillComplete(structuredAc->getDomainMap(), structuredAc->getRangeMap());

  const real_type differenceNorm      = difference->getFrobeniusNorm();
  const real_type referenceNorm       = referenceAc->getFrobeniusNorm();
  const real_type scale               = std::max(referenceNorm, Teuchos::ScalarTraits<real_type>::one());
  const real_type relativeError       = differenceNorm / scale;
  const real_type comparisonTolerance = 1000.0 * Teuchos::ScalarTraits<real_type>::eps();

  TEUCHOS_TEST_FOR_EXCEPTION(
      !(relativeError <= comparisonTolerance),
      std::runtime_error,
      "Final RAP matrices differ: relative Frobenius error = "
          << relativeError << ", tolerance = " << comparisonTolerance);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void runStructuredRAPComparison(const std::string& matrixType,
                                const GlobalOrdinal nx,
                                const GlobalOrdinal ny,
                                const GlobalOrdinal nz,
                                const GlobalOrdinal mx,
                                const GlobalOrdinal my,
                                const GlobalOrdinal mz,
                                const int interpolationOrder,
                                const std::string& coarseningRate,
                                Teuchos::FancyOStream& out) {
  using SC = Scalar;
  using LO = LocalOrdinal;
  using GO = GlobalOrdinal;
  using NO = Node;

  StructuredProblemData<SC, LO, GO, NO> problem =
      buildStructuredProblem<SC, LO, GO, NO>(matrixType, nx, ny, nz, mx, my, mz);

  checkDetectedFineStencil<SC, LO, GO, NO>(problem, matrixType);

  StructuredTransferData<SC, LO, GO, NO> transferData =
      buildStructuredTransferData<SC, LO, GO, NO>(problem, interpolationOrder, coarseningRate);

  Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO> > structuredAc =
      buildCoarseMatrix<SC, LO, GO, NO>(problem, transferData, true);
  Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO> > referenceAc =
      buildCoarseMatrix<SC, LO, GO, NO>(problem, transferData, false);

  compareRAPMatrices<SC, LO, GO, NO>(structuredAc, referenceAc, out);
}

}  // namespace

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredRAPFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<StructuredRAPFactory> rapFactory = rcp(new StructuredRAPFactory);
  TEST_EQUALITY(rapFactory != Teuchos::null, true);

  out << *rapFactory << std::endl;
}  // Constructor test

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredRAPFactory, ConstantLaplace1D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  runStructuredRAPComparison<SC, LO, GO, NO>("Laplace1D", 10 * comm->getSize(), -1, -1,
                                             comm->getSize(), -1, -1, 0, "{3}", out);
}  // ConstantLaplace1D test

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredRAPFactory, ConstantLaplace2D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  runStructuredRAPComparison<SC, LO, GO, NO>("Laplace2D", 12, 12, -1,
                                             -1, -1, -1, 0, "{3}", out);
}  // ConstantLaplace2D test

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredRAPFactory, LinearLaplace2D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  runStructuredRAPComparison<SC, LO, GO, NO>("Laplace2D", 10, 10, -1,
                                             -1, -1, -1, 1, "{2}", out);
}  // LinearLaplace2D test

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredRAPFactory, ConstantLaplace3D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  runStructuredRAPComparison<SC, LO, GO, NO>("Laplace3D", 10, 10, 10,
                                             -1, -1, -1, 0, "{3}", out);
}  // ConstantLaplace3D test

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredRAPFactory, LinearLaplace3D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  runStructuredRAPComparison<SC, LO, GO, NO>("Laplace3D", 10, 10, 10,
                                             -1, -1, -1, 1, "{2}", out);
}  // LinearLaplace3D test

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredRAPFactory, ConstantElasticity2D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  runStructuredRAPComparison<SC, LO, GO, NO>("Elasticity2D", 12, 12, -1,
                                             -1, -1, -1, 0, "{3}", out);
}  // ConstantElasticity2D test

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredRAPFactory, ConstantElasticity3D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  runStructuredRAPComparison<SC, LO, GO, NO>("Elasticity3D", 10, 10, 10,
                                             -1, -1, -1, 0, "{3}", out);
}  // ConstantElasticity3D test

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredRAPFactory, LinearElasticity3D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  runStructuredRAPComparison<SC, LO, GO, NO>("Elasticity3D", 10, 10, 10,
                                             -1, -1, -1, 1, "{2}", out);
}  // LinearElasticity3D test

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredRAPFactory, Constructor, Scalar, LO, GO, Node)          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredRAPFactory, ConstantLaplace1D, Scalar, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredRAPFactory, ConstantLaplace2D, Scalar, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredRAPFactory, LinearLaplace2D, Scalar, LO, GO, Node)      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredRAPFactory, ConstantLaplace3D, Scalar, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredRAPFactory, LinearLaplace3D, Scalar, LO, GO, Node)      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredRAPFactory, ConstantElasticity2D, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredRAPFactory, ConstantElasticity3D, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredRAPFactory, LinearElasticity3D, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
