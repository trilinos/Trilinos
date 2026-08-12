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
#include <sstream>
#include <vector>

#include "MueLu_config.hpp"

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Galeri_XpetraMaps.hpp>

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_CreateXpetraPreconditioner.hpp"
#include "MueLu_StructuredRAPFactory.hpp"

namespace MueLuTests {

namespace {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
struct StructuredProblemData {
  using Matrix                = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using RealValuedMultiVector = Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node>;
  using Vector                = Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  Teuchos::RCP<Matrix> A;
  Teuchos::RCP<RealValuedMultiVector> coordinates;
  Teuchos::RCP<Vector> Mdiag;
  Teuchos::Array<LocalOrdinal> lNodesPerDim;
  std::string matrixType;
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
  using Vector                = typename StructuredProblemData<SC, LO, GO, NO>::Vector;

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
  problem.matrixType   = matrixType;
  problem.lNodesPerDim = Teuchos::Array<LO>(3, Teuchos::as<LO>(1));
  problem.dofsPerNode  = 1;

  Teuchos::RCP<const Map> map;
  if (matrixType == "Laplace1D") {
    problem.numDimensions = 1;
    map = Galeri::Xpetra::CreateMap<LO, GO, NO>(TestHelpers::Parameters::getLib(), "Cartesian1D", comm, galeriList);
    problem.coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<typename RealValuedMultiVector::scalar_type, LO, GO, Map, RealValuedMultiVector>("1D", map, galeriList);
    problem.lNodesPerDim[0] = galeriList.get<LO>("lnx");
  } else if (matrixType == "Laplace2D" || matrixType == "Elasticity2D") {
    problem.numDimensions = 2;
    map = Galeri::Xpetra::CreateMap<LO, GO, NO>(TestHelpers::Parameters::getLib(), "Cartesian2D", comm, galeriList);
    problem.coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<typename RealValuedMultiVector::scalar_type, LO, GO, Map, RealValuedMultiVector>("2D", map, galeriList);
    problem.lNodesPerDim[0] = galeriList.get<LO>("lnx");
    problem.lNodesPerDim[1] = galeriList.get<LO>("lny");
  } else if (matrixType == "Laplace3D" || matrixType == "Elasticity3D") {
    problem.numDimensions = 3;
    map = Galeri::Xpetra::CreateMap<LO, GO, NO>(TestHelpers::Parameters::getLib(), "Cartesian3D", comm, galeriList);
    problem.coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<typename RealValuedMultiVector::scalar_type, LO, GO, Map, RealValuedMultiVector>("3D", map, galeriList);
    problem.lNodesPerDim[0] = galeriList.get<LO>("lnx");
    problem.lNodesPerDim[1] = galeriList.get<LO>("lny");
    problem.lNodesPerDim[2] = galeriList.get<LO>("lnz");
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Unsupported StructuredRAPFactory test matrix type: " << matrixType);
  }

  if (matrixType == "Elasticity2D") {
    problem.dofsPerNode = 2;
    map = Xpetra::MapFactory<LO, GO, NO>::Build(map, problem.dofsPerNode);
  } else if (matrixType == "Elasticity3D") {
    problem.dofsPerNode = 3;
    map = Xpetra::MapFactory<LO, GO, NO>::Build(map, problem.dofsPerNode);
  }

  Teuchos::RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > galeriProblem =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixType, map, galeriList);
  problem.A = galeriProblem->BuildMatrix();
  if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D")
    problem.A->SetFixedBlockSize(problem.dofsPerNode);

  problem.Mdiag = Xpetra::VectorFactory<SC, LO, GO, NO>::Build(problem.A->getRowMap(), false);
  problem.A->getLocalDiagCopy(*problem.Mdiag);

  return problem;
}

inline Teuchos::ParameterList buildStructuredRAPParameterList(const int dofsPerNode,
                                                             const int interpolationOrder,
                                                             const std::string& coarseningRate,
                                                             const bool prebuildCoarseGraph) {
  Teuchos::ParameterList paramList;
  paramList.sublist("Matrix").set("PDE equations", dofsPerNode);

  Teuchos::ParameterList& factories = paramList.sublist("Factories");
  factories.sublist("myAmalgamationFact")
      .set("factory", "AmalgamationFactory");

  Teuchos::ParameterList& coalesceDrop = factories.sublist("myCoalesceDropFact");
  coalesceDrop.set("factory", "CoalesceDropFactory");
  coalesceDrop.set("lightweight wrap", true);
  coalesceDrop.set("aggregation: drop tol", 0.0);
  coalesceDrop.set("UnAmalgamationInfo", "myAmalgamationFact");

  Teuchos::ParameterList& aggregation = factories.sublist("myAggregationFact");
  aggregation.set("factory", "StructuredAggregationFactory");
  aggregation.set("aggregation: coupling", "uncoupled");
  aggregation.set("aggregation: output type", "CrsGraph");
  aggregation.set("aggregation: coarsening order", interpolationOrder);
  aggregation.set("aggregation: coarsening rate", coarseningRate);
  aggregation.set("Graph", "myCoalesceDropFact");

  Teuchos::ParameterList& coarseMap = factories.sublist("myCoarseMapFact");
  coarseMap.set("factory", "CoarseMapFactory");
  coarseMap.set("Aggregates", "myAggregationFact");

  Teuchos::ParameterList& prolongator = factories.sublist("myProlongatorFact");
  prolongator.set("factory", "GeometricInterpolationPFactory");
  prolongator.set("interp: build coarse coordinates", true);
  prolongator.set("structuredInterpolationOrder", "myAggregationFact");
  prolongator.set("prolongatorGraph", "myAggregationFact");
  prolongator.set("indexManager", "myAggregationFact");
  prolongator.set("coarseCoordinatesFineMap", "myAggregationFact");
  prolongator.set("coarseCoordinatesMap", "myAggregationFact");

  Teuchos::ParameterList& coordinatesTransfer = factories.sublist("myCoordTransferFact");
  coordinatesTransfer.set("factory", "CoordinatesTransferFactory");
  coordinatesTransfer.set("structured aggregation", true);
  coordinatesTransfer.set("numDimensions", "myAggregationFact");
  coordinatesTransfer.set("lCoarseNodesPerDim", "myAggregationFact");

  Teuchos::ParameterList& nullspace = factories.sublist("myNullspaceFact");
  nullspace.set("factory", "NullspaceFactory");
  nullspace.set("Nullspace", "myProlongatorFact");
  nullspace.set("CoarseMap", "myCoarseMapFact");

  Teuchos::ParameterList& rap = factories.sublist("myRAPFact");
  rap.set("factory", "StructuredRAPFactory");
  rap.set("P", "myProlongatorFact");
  rap.set("lCoarseNodesPerDim", "myAggregationFact");
  rap.set("structuredInterpolationOrder", "myAggregationFact");
  rap.set("rap: triple product", true);
  rap.set("rap: prebuild coarse graph", prebuildCoarseGraph);
  rap.set("transpose: use implicit", true);
  rap.sublist("TransferFactories").set("CoordinateTransfer", "myCoordTransferFact");

  Teuchos::ParameterList& smoother = factories.sublist("myMTSGS");
  smoother.set("factory", "TrilinosSmoother");
  smoother.set("type", "RELAXATION");
  Teuchos::ParameterList& smootherParams = smoother.sublist("ParameterList");
  smootherParams.set("relaxation: type", "MT Symmetric Gauss-Seidel");
  smootherParams.set("relaxation: symmetric matrix structure", true);
  smootherParams.set("relaxation: sweeps", 2);
  smootherParams.set("relaxation: use l1", true);

  Teuchos::ParameterList& hierarchy = paramList.sublist("Hierarchy");
  hierarchy.set("max levels", 2);
  hierarchy.set("cycle type", "V");
  hierarchy.set("coarse: max size", 3);
  hierarchy.set("verbosity", "None");

  Teuchos::ParameterList& all = hierarchy.sublist("All");
  all.set("PreSmoother", "myMTSGS");
  all.set("PostSmoother", "myMTSGS");
  all.set("Graph", "myCoalesceDropFact");
  all.set("Nullspace", "myNullspaceFact");
  all.set("Aggregates", "myAggregationFact");
  all.set("lCoarseNodesPerDim", "myAggregationFact");
  all.set("P", "myProlongatorFact");
  all.set("A", "myRAPFact");
  all.set("CoarseSolver", "myMTSGS");
  all.set("Coordinates", "myProlongatorFact");
  all.set("lNodesPerDim", "myCoordTransferFact");
  all.set("numDimensions", "myCoordTransferFact");

  return paramList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
buildFinalCoarseMatrix(const StructuredProblemData<Scalar, LocalOrdinal, GlobalOrdinal, Node>& problem,
                       Teuchos::ParameterList paramList) {
  using SC                    = Scalar;
  using LO                    = LocalOrdinal;
  using GO                    = GlobalOrdinal;
  using NO                    = Node;
  using Matrix                = Xpetra::Matrix<SC, LO, GO, NO>;
  using RealValuedMultiVector = typename StructuredProblemData<SC, LO, GO, NO>::RealValuedMultiVector;

  Teuchos::ParameterList& userParamList = paramList.sublist("user data");
  userParamList.set<int>("int numDimensions", problem.numDimensions);
  userParamList.set<Teuchos::Array<LO> >("Array<LO> lNodesPerDim", problem.lNodesPerDim);
  userParamList.set<Teuchos::RCP<RealValuedMultiVector> >("Coordinates", problem.coordinates);
  userParamList.set<double>("double cfl", 1.0);
  userParamList.set<double>("double deltaT", 1.0);
  userParamList.set("Mdiag", problem.Mdiag);
  userParamList.set<std::string>("string matrixType", problem.matrixType);

  problem.A->SetMaxEigenvalueEstimate(-Teuchos::ScalarTraits<SC>::one());
  Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO> > hierarchy =
      MueLu::CreateXpetraPreconditioner<SC, LO, GO, NO>(problem.A, paramList);

  Teuchos::RCP<Matrix> Ac;
  hierarchy->GetLevel(hierarchy->GetNumLevels() - 1)->Get("A", Ac);
  return Ac;
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
  const double graphOverhead = referenceNnz == 0
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

  const Teuchos::RCP<const Map> structuredRowMap = structuredAc->getRowMap();
  const Teuchos::RCP<const Map> structuredColMap = structuredAc->getColMap();
  const Teuchos::RCP<const Map> referenceColMap  = referenceAc->getColMap();
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
  TEUCHOS_TEST_FOR_EXCEPTION(
      globalGraphMismatch != 0,
      std::runtime_error,
      localGraphMismatch
          ? graphMismatchMessage + ": " + localMismatchReason
          : graphMismatchMessage + " on another MPI rank.");

  Teuchos::RCP<Matrix> difference;
  Xpetra::MatrixMatrix<SC, LO, GO, NO>::TwoMatrixAdd(
      *structuredAc, false, TST::one(),
      *referenceAc, false, -TST::one(),
      difference, out);

  if (!difference->isFillComplete())
    difference->fillComplete(structuredAc->getDomainMap(), structuredAc->getRangeMap());

  const real_type differenceNorm = difference->getFrobeniusNorm();
  const real_type referenceNorm  = referenceAc->getFrobeniusNorm();
  const real_type scale = std::max(referenceNorm, Teuchos::ScalarTraits<real_type>::one());
  const real_type relativeError = differenceNorm / scale;
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

  Teuchos::ParameterList structuredParamList =
      buildStructuredRAPParameterList(problem.dofsPerNode, interpolationOrder, coarseningRate, true);
  Teuchos::ParameterList referenceParamList = structuredParamList;
  referenceParamList.sublist("Factories").sublist("myRAPFact").set("rap: prebuild coarse graph", false);

  Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO> > structuredAc =
      buildFinalCoarseMatrix<SC, LO, GO, NO>(problem, structuredParamList);
  Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO> > referenceAc =
      buildFinalCoarseMatrix<SC, LO, GO, NO>(problem, referenceParamList);

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

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredRAPFactory, Constructor, Scalar, LO, GO, Node)   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredRAPFactory, ConstantLaplace1D, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredRAPFactory, ConstantLaplace2D, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredRAPFactory, LinearLaplace2D, Scalar, LO, GO, Node)   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredRAPFactory, ConstantLaplace3D, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredRAPFactory, LinearLaplace3D, Scalar, LO, GO, Node)   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredRAPFactory, ConstantElasticity2D, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredRAPFactory, ConstantElasticity3D, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredRAPFactory, LinearElasticity3D, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
