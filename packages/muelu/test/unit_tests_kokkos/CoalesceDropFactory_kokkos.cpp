// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <iomanip>
#include <ostream>

#include "MueLu_TestHelpers_kokkos.hpp"
#include "MueLu_Version.hpp"

#include <MueLu_InitialBlockNumberFactory.hpp>
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_FilteredAFactory.hpp"
#include "MueLu_CoalesceDropFactory_kokkos.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "Tpetra_Access.hpp"
#include "Teuchos_Assert.hpp"

#include <Galeri_XpetraParameters.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<CoalesceDropFactory_kokkos> coalesceDropFact = rcp(new CoalesceDropFactory_kokkos());
  TEST_EQUALITY(coalesceDropFact != Teuchos::null, true);

  out << *coalesceDropFact << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class CoalesceDropTestingHelper {
 public:
  using magnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;

  CoalesceDropTestingHelper(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &A_,
                            RCP<Xpetra::MultiVector<magnitudeType, LocalOrdinal, GlobalOrdinal, Node>> &coordinates_,
                            RCP<ParameterList> &params_)
    : A(A_)
    , coordinates(coordinates_)
    , params(params_) {}

  CoalesceDropTestingHelper(GlobalOrdinal nx,
                            RCP<ParameterList> &params_) {
    using Map                   = Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
    using RealValuedMultiVector = Xpetra::MultiVector<magnitudeType, LocalOrdinal, GlobalOrdinal, Node>;

    A = TestHelpers_kokkos::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build1DPoisson(nx);

    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<Scalar, LocalOrdinal, GlobalOrdinal, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);

    params = params_;
  }

  void replaceEntriesInA(std::vector<std::tuple<GlobalOrdinal, GlobalOrdinal, Scalar>> replace) {
    const auto INVALID = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();

    auto crsA = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(A, true)->getCrsMatrix();
    crsA->resumeFill();

    auto rowmap = A->getRowMap();
    auto colmap = A->getColMap();
    auto lclA   = A->getLocalMatrixHost();

    size_t found = 0;
    for (auto [row_gid, col_gid, value] : replace) {
      auto row_lid = rowmap->getLocalElement(row_gid);
      auto col_lid = colmap->getLocalElement(col_gid);
      if ((row_lid != INVALID) && (col_lid != INVALID)) {
        auto row = lclA.row(row_lid);
        for (LocalOrdinal k = 0; k < row.length; ++k) {
          if (row.colidx(k) == col_lid) {
            row.value(k) = value;
            ++found;
          }
        }
      }
    }
    crsA->fillComplete();

    auto comm = rowmap->getComm();
    size_t globalFound;
    MueLu_sumAll(comm, found, globalFound);
    TEUCHOS_ASSERT_EQUALITY(globalFound, replace.size());
  }

  void print(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &A_, Teuchos::FancyOStream &out) {
    auto rowmap          = A_->getRowMap();
    auto comm            = rowmap->getComm();
    auto single_proc_map = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(rowmap->lib(),
                                                                                        rowmap->getGlobalNumElements(),
                                                                                        comm->getRank() == 0 ? rowmap->getGlobalNumElements() : 0,
                                                                                        rowmap->getIndexBase(),
                                                                                        comm);
    auto importer        = Xpetra::ImportFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(rowmap, single_proc_map);
    auto single_proc_A   = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(A_, *importer);

    if (comm->getRank() == 0) {
      auto lclA = single_proc_A->getLocalMatrixHost();
      std::stringstream ss;
      ss << std::fixed << std::setprecision(2);

      ss << std::setw(5) << " ";
      for (LocalOrdinal col_id = 0; col_id < lclA.numCols(); ++col_id) {
        ss << " " << std::setw(6) << A->getColMap()->getIndexBase() + col_id;
      }
      ss << std::endl;
      ss << std::setw(5) << " " << std::setw(7 * lclA.numCols()) << std::setfill('-') << "";
      ss << std::endl;
      ss << std::setfill(' ');
      for (LocalOrdinal row_id = 0; row_id < lclA.numRows(); ++row_id) {
        auto row = lclA.rowConst(row_id);
        ss << std::setw(3) << rowmap->getIndexBase() + row_id << " |";
        LocalOrdinal K = 0;
        for (LocalOrdinal k = 0; k < lclA.numRows(); ++k) {
          if (row.colidx(K) == k) {
            ss << " " << std::setw(6) << row.value(K);
            ++K;
          } else {
            ss << " " << std::setw(6) << "";
          }
        }
        ss << std::endl;
      }
      out << ss.str();
    }
  }

  void printA(Teuchos::FancyOStream &out) {
    print(A, out);
  }

  void printFilteredA(Teuchos::FancyOStream &out) {
    TEUCHOS_ASSERT(!filteredA.is_null());
    print(filteredA, out);
  }

  void run() {
#include <MueLu_UseShortNames.hpp>
    Level fineLevel;
    TestHelpers_kokkos::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(fineLevel);
    fineLevel.Set("A", A);
    fineLevel.Set("Coordinates", coordinates);

    RCP<MueLu::SingleLevelFactoryBase> dropFact;
    RCP<MueLu::SingleLevelFactoryBase> filteredAFact;
    dropFact = rcp(new CoalesceDropFactory_kokkos());
    dropFact->SetParameterList(*params);
    filteredAFact = dropFact;

    fineLevel.Request("A", filteredAFact.get());
    fineLevel.Request("Graph", dropFact.get());
    fineLevel.Request("DofsPerNode", dropFact.get());
    dropFact->Build(fineLevel);

    filteredA = fineLevel.Get<RCP<Matrix>>("A", dropFact.get());

    dofsPerNode = fineLevel.Get<LocalOrdinal>("DofsPerNode", dropFact.get());

    auto graph_d = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", dropFact.get());
    graph        = graph_d->copyToHost();

    auto boundaryNodes_d = graph_d->GetBoundaryNodeMap();
    boundaryNodes        = Kokkos::View<bool *, Kokkos::HostSpace>("boundaryNodes_host", boundaryNodes_d.extent(0));
    Kokkos::deep_copy(boundaryNodes, boundaryNodes_d);
  }

  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A;
  RCP<Xpetra::MultiVector<magnitudeType, LocalOrdinal, GlobalOrdinal, Node>> coordinates;
  RCP<Teuchos::ParameterList> params;

  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> filteredA;
  RCP<MueLu::LWGraph<LocalOrdinal, GlobalOrdinal, Node>> graph;
  LocalOrdinal dofsPerNode;
  Kokkos::View<bool *, Kokkos::HostSpace> boundaryNodes;
};

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);
}  // Build

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, DistanceLaplacian, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = TestHelpers_kokkos::Parameters::getDefaultComm();
  auto lib                           = TestHelpers_kokkos::Parameters::getLib();

  // 1 dof per node
  {
    Level fineLevel;
    TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

    RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
    fineLevel.Set("A", A);

    Teuchos::ParameterList galeriList;
    galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
    RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
    fineLevel.Set("Coordinates", coordinates);

    CoalesceDropFactory_kokkos coalesceDropFact;
    coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
    // We're dropping all the interior off-diagonal entries.
    // dx = 1/36
    // L_ij = -36
    // L_ii = 72
    // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
    coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.51));
    coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
    fineLevel.Request("Graph", &coalesceDropFact);
    fineLevel.Request("DofsPerNode", &coalesceDropFact);

    coalesceDropFact.Build(fineLevel);

    RCP<LWGraph_kokkos> graph_d = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
    auto graph                  = graph_d->copyToHost();
    LO myDofsPerNode            = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
    TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

    const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

    TEST_EQUALITY(graph->GetGlobalNumEdges(), 40);
  }

  // 3 dof per node
  {
    Level fineLevel;
    TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

    GO nx = 11, ny = 10, nz = 10;
    double stretchx = 9.0, stretchy = 1.0, stretchz = 1.0;

    Teuchos::ParameterList galeriList;
    galeriList.set("nx", Teuchos::as<GlobalOrdinal>(nx));
    galeriList.set("ny", Teuchos::as<GlobalOrdinal>(ny));
    galeriList.set("nz", Teuchos::as<GlobalOrdinal>(nz));

    auto map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian3D", comm, galeriList);
    auto coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<real_type, LO, GO, Map, RealValuedMultiVector>("3D", map, galeriList);

    map = Xpetra::MapFactory<LO, GO, Node>::Build(map, 3);

    galeriList.set("stretchx", stretchx);
    galeriList.set("stretchy", stretchy);
    galeriList.set("stretchz", stretchz);

    galeriList.set("left boundary", "Dirichlet");
    galeriList.set("right boundary", "Neumann");
    galeriList.set("bottom boundary", "Neumann");
    galeriList.set("top boundary", "Neumann");
    galeriList.set("front boundary", "Neumann");
    galeriList.set("back boundary", "Neumann");
    galeriList.set("keepBCs", true);

    RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector>> Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Elasticity3D", map, galeriList);

    RCP<Matrix> A = Pr->BuildMatrix();
    A->SetFixedBlockSize(3);

    RCP<RealValuedMultiVector> newcoordinates = Pr->BuildCoords();

    // Galeri makes multiple copies of coordinates to deal with
    // some issues when Ndofs != Nmeshnodes
    for (size_t kkk = 0; kkk < coordinates->getNumVectors(); kkk++) {
      Teuchos::ArrayRCP<real_type> old     = coordinates->getDataNonConst(kkk);
      Teuchos::ArrayRCP<real_type> newvals = newcoordinates->getDataNonConst(kkk);
      int numCopies                        = newvals.size() / old.size();
      for (int jj = 0; jj < old.size(); jj++) old[jj] = newvals[numCopies * (Teuchos_Ordinal)jj];
    }

    fineLevel.Set("A", A);
    fineLevel.Set("Coordinates", coordinates);

    CoalesceDropFactory_kokkos coalesceDropFact;
    coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
    // We're dropping all the edges in weak ny direction and keep the ones in nx
    // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
    coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.04));
    coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
    fineLevel.Request("Graph", &coalesceDropFact);
    fineLevel.Request("DofsPerNode", &coalesceDropFact);

    coalesceDropFact.Build(fineLevel);

    LO myDofsPerNode            = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
    RCP<LWGraph_kokkos> graph_d = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
    auto graph                  = graph_d->copyToHost();
    auto boundaryNodeMap        = graph->GetBoundaryNodeMap();
    auto numGlobalEdges         = graph->GetGlobalNumEdges();

    GO numLocalBoundaryNodes  = 0;
    GO numGlobalBoundaryNodes = 0;
    for (size_t i = 0; i < boundaryNodeMap.size(); ++i)
      if (boundaryNodeMap(i)) numLocalBoundaryNodes++;
    MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);

    TEST_EQUALITY(myDofsPerNode, 3);
    TEST_EQUALITY(numGlobalBoundaryNodes, 100);
    TEST_EQUALITY(numGlobalEdges, 7940);
  }
}  // DistanceLaplacian

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, DistanceLaplacianScaledCut, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  // Now we doctor the coordinates so that the off-diagonal pair row 0 will want to keep (0,1) and row 1 will want to drop (1,0)
  if (comm->getRank() == 0) {
    auto vals = coordinates->getDataNonConst(0);
    vals[0]   = vals[0] - 2000 * 36;
  }

  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0 / 8.0));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: distance laplacian algo", Teuchos::ParameterEntry(std::string("scaled cut")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph_d = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  auto graph                  = graph_d->copyToHost();
  LO myDofsPerNode            = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 105);

}  // DistanceLaplacianScaledCut

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, DistanceLaplacianUnscaledCut, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);

  // Now we doctor the coordinates so that the off-diagonal pair row 0 will want to keep (0,1) and row 1 will want to drop (1,0)
  if (!comm->getRank()) {
    auto vals = coordinates->getDataNonConst(0);
    vals[0]   = vals[0] - 2000 * 36;
  }

  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0 / 8.0));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: distance laplacian algo", Teuchos::ParameterEntry(std::string("unscaled cut")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph_d = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  auto graph                  = graph_d->copyToHost();
  LO myDofsPerNode            = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 105);

}  // DistanceLaplacianUnscaleCut

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, DistanceLaplacianCutSym, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);

  // Now we doctor the coordinates so that the off-diagonal pair row 0 will want to keep (0,1) and row 1 will want to drop (1,0)
  if (!comm->getRank()) {
    auto vals = coordinates->getDataNonConst(0);
    vals[0]   = vals[0] - 2000 * 36;
  }

  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);

  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.5));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: distance laplacian algo", Teuchos::ParameterEntry(std::string("scaled cut symmetric")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph_d = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  auto graph                  = graph_d->copyToHost();
  LO myDofsPerNode            = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 106);

}  // DistanceLaplacianCutScaled

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, DistanceLaplacianScalarMaterial, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  auto material = MultiVectorFactory::Build(A->getRowMap(), 1);
  material->putScalar(1.0);
  fineLevel.Set("Material", material);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: distance laplacian metric", Teuchos::ParameterEntry(std::string("material")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph_kokkos = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  auto graph                       = graph_kokkos->copyToHost();

  LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 40);

}  // DistanceLaplacianScalarMaterial

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
using materialTestCase = std::tuple<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>,
                                    RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>,
                                    RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node>>,
                                    RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
materialTestCase<Scalar, LocalOrdinal, GlobalOrdinal, Node> constructVariableMaterialMatrix(RCP<const Teuchos::Comm<int>> &comm) {
#include <MueLu_UseShortNames.hpp>
#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<Scalar>;
#else
  using ATS     = Kokkos::ArithTraits<Scalar>;
#endif
  using impl_scalar_type = typename ATS::val_type;
#if KOKKOS_VERSION >= 40799
  using implATS = KokkosKernels::ArithTraits<impl_scalar_type>;
#else
  using implATS = Kokkos::ArithTraits<impl_scalar_type>;
#endif
  using magnitudeType = typename implATS::magnitudeType;

  RCP<const Map> map = MapFactory::Build(Xpetra::UseTpetra, 27 * comm->getSize(), 0, comm);

  RCP<Matrix> A;
  {
    using local_matrix_type              = typename Matrix::local_matrix_device_type;
    std::vector<LocalOrdinal> rowptr     = {0, 8, 20, 37, 49, 61, 79, 106, 124, 132, 143, 155, 173, 185, 193, 211, 223, 231, 243, 251, 263, 281, 293, 301, 313, 325, 333, 341};
    std::vector<LocalOrdinal> indices    = {0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 18, 19, 20, 21, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15, 18, 19, 20, 21, 24, 25, 1, 2, 5, 6, 8, 9, 10, 11, 1, 5, 6, 8, 9, 10, 11, 12, 14, 16, 17, 1, 2, 5, 6, 8, 9, 10, 11, 19, 20, 22, 23, 1, 2, 5, 6, 8, 9, 10, 11, 12, 14, 16, 17, 19, 20, 22, 23, 24, 26, 2, 3, 6, 7, 9, 11, 12, 13, 14, 15, 16, 17, 2, 3, 6, 7, 12, 13, 14, 15, 2, 3, 6, 7, 9, 11, 12, 13, 14, 15, 16, 17, 20, 21, 23, 24, 25, 26, 2, 3, 6, 7, 12, 13, 14, 15, 20, 21, 24, 25, 2, 6, 9, 11, 12, 14, 16, 17, 2, 6, 9, 11, 12, 14, 16, 17, 20, 23, 24, 26, 4, 5, 6, 7, 18, 19, 20, 21, 4, 5, 6, 7, 10, 11, 18, 19, 20, 21, 22, 23, 4, 5, 6, 7, 10, 11, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 4, 5, 6, 7, 14, 15, 18, 19, 20, 21, 24, 25, 5, 6, 10, 11, 19, 20, 22, 23, 5, 6, 10, 11, 14, 17, 19, 20, 22, 23, 24, 26, 6, 7, 11, 14, 15, 17, 20, 21, 23, 24, 25, 26, 6, 7, 14, 15, 20, 21, 24, 25, 6, 11, 14, 17, 20, 23, 24, 26};
    std::vector<impl_scalar_type> values = {42.6666666666667, -4.88498130835069e-15, -10.6666666666667, 4.44089209850063e-16, -7.105427357601e-15, -10.6666666666667, -10.6666666666667, -10.6666666666667, -4.88498130835069e-15, 85.3333333333333, -1.77635683940025e-15, -10.6666666666667, -10.6666666666667, -1.31006316905768e-14, -21.3333333333333, -10.6666666666667, -8.88178419700125e-16, -10.6666666666667, -10.6666666666667, -10.6666666666667, -10.6666666666667, -1.77635683940025e-15, 170.666666666667, -7.105427357601e-15, -10.6666666666667, -21.3333333333333, -2.44526621173691e-14, -21.3333333333333, -10.6666666666667, -10.6666666666667, -21.3333333333333, -1.02140518265514e-14, -10.6666666666667, -21.3333333333333, -10.6666666666667, -10.6666666666667, -10.6666666666667, 4.44089209850063e-16, -10.6666666666667, -7.105427357601e-15, 85.3333333333333, -10.6666666666667, -10.6666666666667, -21.3333333333333, -1.06581410364015e-14, -10.6666666666667, -2.66453525910038e-15, -10.6666666666667, -10.6666666666667, -7.105427357601e-15, -10.6666666666667, -10.6666666666667, -10.6666666666667, 43, 6.10622663543836e-16, -10.75, 8.81239525796218e-16, -2.42861286636753e-17, -0.0833333333333333, -0.0833333333333333, -0.0833333333333333, -10.6666666666667, -1.31006316905768e-14, -21.3333333333333, -10.6666666666667, 6.10622663543836e-16, 86, -2.66626998257635e-15, -10.75, -10.6666666666667, -10.6666666666667, -9.08995101411847e-16, -10.75, -0.0833333333333333, -3.46944695195361e-17, -0.166666666666667, -0.0833333333333333, -0.0833333333333333, -0.0833333333333333, -10.6666666666667, -21.3333333333333, -2.44526621173691e-14, -21.3333333333333, -10.75, -2.66626998257635e-15, 172, -6.65092980689508e-15, -10.6666666666667, -21.3333333333333, -10.75, -7.63278329429795e-17, -21.3333333333333, -10.6666666666667, -1.03146657881581e-14, -10.75, -10.6666666666667, -10.75, -0.0833333333333333, -0.166666666666667, -4.68375338513738e-17, -0.166666666666667, -0.0833333333333333, -0.166666666666667, -0.166666666666667, -0.0833333333333333, -0.0833333333333333, -10.6666666666667, -10.6666666666667, -21.3333333333333, -1.06581410364015e-14, 8.81239525796218e-16, -10.75, -6.65092980689508e-15, 86, -10.6666666666667, -10.6666666666667, -10.75, -3.59434704222394e-15, -0.0833333333333333, -0.0833333333333333, -0.166666666666667, -1.21430643318376e-17, -0.0833333333333333, -0.0833333333333333, -8.88178419700125e-16, -10.6666666666667, -10.6666666666667, -10.6666666666667, 42.6666666666667, -1.77635683940025e-15, -7.105427357601e-15, -10.6666666666667, -10.6666666666667, -10.6666666666667, -21.3333333333333, -1.77635683940025e-15, 85.3333333333333, -10.6666666666667, -1.33781874467331e-14, -10.6666666666667, -10.6666666666667, -6.21724893790088e-15, -10.6666666666667, -10.6666666666667, -10.6666666666667, -9.08995101411847e-16, -10.75, -7.105427357601e-15, -10.6666666666667, 43, -1.8249290967276e-15, -0.0833333333333333, -0.0833333333333333, -2.08166817117217e-17, -0.0833333333333333, -10.6666666666667, -21.3333333333333, -10.75, -7.63278329429795e-17, -10.6666666666667, -1.33781874467331e-14, -1.8249290967276e-15, 86, -10.6666666666667, -10.75, -10.6666666666667, -4.9960036108132e-15, -0.0833333333333333, -0.166666666666667, -0.0833333333333333, -3.46944695195361e-17, -0.0833333333333333, -0.0833333333333333, -1.02140518265514e-14, -10.6666666666667, -21.3333333333333, -10.6666666666667, -10.6666666666667, -10.6666666666667, 85.3333333333333, -3.5527136788005e-15, -1.06581410364015e-14, -10.6666666666667, 8.88178419700125e-16, -10.6666666666667, -10.6666666666667, -2.66453525910038e-15, -10.6666666666667, -10.6666666666667, -3.5527136788005e-15, 42.6666666666667, -10.6666666666667, -3.33066907387547e-15, -21.3333333333333, -10.6666666666667, -1.03146657881581e-14, -10.75, -10.6666666666667, -10.75, -1.06581410364015e-14, -10.6666666666667, 86, -4.48252546192407e-15, -10.6666666666667, 4.37150315946155e-16, -0.166666666666667, -0.0833333333333333, -0.0833333333333333, -1.7130394325271e-17, -0.0833333333333333, -0.0833333333333333, -10.6666666666667, -10.6666666666667, -10.75, -3.59434704222394e-15, -10.6666666666667, -3.33066907387547e-15, -4.48252546192407e-15, 43, -0.0833333333333333, -0.0833333333333333, -0.0833333333333333, 6.93889390390723e-18, -10.6666666666667, -10.6666666666667, -6.21724893790088e-15, -10.6666666666667, 8.88178419700125e-16, -10.6666666666667, 42.6666666666667, -6.27276008913213e-15, -10.6666666666667, -10.75, -10.6666666666667, -4.9960036108132e-15, -10.6666666666667, 4.37150315946155e-16, -6.27276008913213e-15, 43, -0.0833333333333333, -0.0833333333333333, -0.0833333333333333, -1.38777878078145e-17, -2.42861286636753e-17, -0.0833333333333333, -0.0833333333333333, -0.0833333333333333, 0.333333333333333, -2.08166817117217e-17, -0.0833333333333333, -1.38777878078145e-17, -0.0833333333333333, -3.46944695195361e-17, -0.166666666666667, -0.0833333333333333, -0.0833333333333333, -0.0833333333333333, -2.08166817117217e-17, 0.666666666666667, -3.46944695195361e-17, -0.0833333333333333, -2.08166817117217e-17, -0.0833333333333333, -0.0833333333333333, -0.166666666666667, -4.68375338513738e-17, -0.166666666666667, -0.0833333333333333, -0.166666666666667, -0.166666666666667, -0.0833333333333333, -0.0833333333333333, -0.0833333333333333, -3.46944695195361e-17, 1.33333333333333, -5.98479599211998e-17, -0.0833333333333333, -3.16587034365767e-17, -1.04517089427603e-16, -0.0833333333333333, -0.0833333333333333, -0.0833333333333333, -0.0833333333333333, -0.166666666666667, -1.21430643318376e-17, -0.0833333333333333, -0.0833333333333333, -1.38777878078145e-17, -0.0833333333333333, -5.98479599211998e-17, 0.666666666666667, -0.0833333333333333, -4.85722573273506e-17, -0.0833333333333333, -0.0833333333333333, -2.08166817117217e-17, -0.0833333333333333, -2.08166817117217e-17, -0.0833333333333333, 0.333333333333333, -2.42861286636753e-17, -0.0833333333333333, -0.166666666666667, -0.0833333333333333, -3.46944695195361e-17, -0.0833333333333333, -0.0833333333333333, -0.0833333333333333, -3.16587034365767e-17, -2.42861286636753e-17, 0.666666666666667, -0.0833333333333333, -5.24753851482984e-17, -0.166666666666667, -0.0833333333333333, -0.0833333333333333, -1.7130394325271e-17, -0.0833333333333333, -0.0833333333333333, -1.04517089427603e-16, -0.0833333333333333, -0.0833333333333333, 0.666666666666667, -5.29090660172926e-17, -3.90312782094782e-18, -0.0833333333333333, -0.0833333333333333, -0.0833333333333333, 6.93889390390723e-18, -0.0833333333333333, -4.85722573273506e-17, -5.29090660172926e-17, 0.333333333333333, -0.0833333333333333, -0.0833333333333333, -0.0833333333333333, -1.38777878078145e-17, -0.0833333333333333, -5.24753851482984e-17, -3.90312782094782e-18, 0.333333333333333};
    local_matrix_type lclMatrix("lclA", 27, 27, 341, values.data(), rowptr.data(), indices.data());
    A = MatrixFactory::Build(map, map, lclMatrix);
  }

  RCP<Matrix> droppedA;
  {
    using local_matrix_type              = typename Matrix::local_matrix_device_type;
    std::vector<LocalOrdinal> rowptr     = {0, 8, 20, 37, 49, 58, 70, 88, 100, 108, 119, 128, 140, 152, 160, 172, 181, 189, 198, 203, 209, 218, 224, 229, 235, 241, 246, 251};
    std::vector<LocalOrdinal> indices    = {0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 18, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15, 1, 2, 5, 6, 8, 9, 10, 11, 1, 5, 6, 8, 9, 10, 11, 12, 14, 16, 17, 1, 2, 5, 6, 8, 9, 10, 11, 22, 1, 2, 5, 6, 8, 9, 10, 11, 12, 14, 16, 17, 2, 3, 6, 7, 9, 11, 12, 13, 14, 15, 16, 17, 2, 3, 6, 7, 12, 13, 14, 15, 2, 3, 6, 7, 9, 11, 12, 13, 14, 15, 16, 17, 2, 3, 6, 7, 12, 13, 14, 15, 25, 2, 6, 9, 11, 12, 14, 16, 17, 2, 6, 9, 11, 12, 14, 16, 17, 26, 4, 18, 19, 20, 21, 18, 19, 20, 21, 22, 23, 18, 19, 20, 21, 22, 23, 24, 25, 26, 18, 19, 20, 21, 24, 25, 10, 19, 20, 22, 23, 19, 20, 22, 23, 24, 26, 20, 21, 23, 24, 25, 26, 15, 20, 21, 24, 25, 17, 20, 23, 24, 26};
    std::vector<impl_scalar_type> values = {42.6666666666667, -4.88498130835069e-15, -10.6666666666667, 4.44089209850063e-16, -7.105427357601e-15, -10.6666666666667, -10.6666666666667, -10.6666666666667, -4.88498130835069e-15, 85.3333333333333, -1.77635683940025e-15, -10.6666666666667, -10.6666666666667, -1.31006316905768e-14, -21.3333333333333, -10.6666666666667, -8.88178419700125e-16, -10.6666666666667, -10.6666666666667, -10.6666666666667, -10.6666666666667, -1.77635683940025e-15, 170.666666666667, -7.105427357601e-15, -10.6666666666667, -21.3333333333333, -2.44526621173691e-14, -21.3333333333333, -10.6666666666667, -10.6666666666667, -21.3333333333333, -1.02140518265514e-14, -10.6666666666667, -21.3333333333333, -10.6666666666667, -10.6666666666667, -10.6666666666667, 4.44089209850063e-16, -10.6666666666667, -7.105427357601e-15, 85.3333333333333, -10.6666666666667, -10.6666666666667, -21.3333333333333, -1.06581410364015e-14, -10.6666666666667, -2.66453525910038e-15, -10.6666666666667, -10.6666666666667, -7.105427357601e-15, -10.6666666666667, -10.6666666666667, -10.6666666666667, 43, 6.10622663543836e-16, -10.75, 8.81239525796218e-16, -2.42861286636753e-17, -10.6666666666667, -1.31006316905768e-14, -21.3333333333333, -10.6666666666667, 6.10622663543836e-16, 86, -2.66626998257635e-15, -10.75, -10.6666666666667, -10.6666666666667, -9.08995101411847e-16, -10.75, -10.6666666666667, -21.3333333333333, -2.44526621173691e-14, -21.3333333333333, -10.75, -2.66626998257635e-15, 172, -6.65092980689508e-15, -10.6666666666667, -21.3333333333333, -10.75, -7.63278329429795e-17, -21.3333333333333, -10.6666666666667, -1.03146657881581e-14, -10.75, -10.6666666666667, -10.75, -10.6666666666667, -10.6666666666667, -21.3333333333333, -1.06581410364015e-14, 8.81239525796218e-16, -10.75, -6.65092980689508e-15, 86, -10.6666666666667, -10.6666666666667, -10.75, -3.59434704222394e-15, -8.88178419700125e-16, -10.6666666666667, -10.6666666666667, -10.6666666666667, 42.6666666666667, -1.77635683940025e-15, -7.105427357601e-15, -10.6666666666667, -10.6666666666667, -10.6666666666667, -21.3333333333333, -1.77635683940025e-15, 85.3333333333333, -10.6666666666667, -1.33781874467331e-14, -10.6666666666667, -10.6666666666667, -6.21724893790088e-15, -10.6666666666667, -10.6666666666667, -10.6666666666667, -9.08995101411847e-16, -10.75, -7.105427357601e-15, -10.6666666666667, 43, -1.8249290967276e-15, -2.08166817117217e-17, -10.6666666666667, -21.3333333333333, -10.75, -7.63278329429795e-17, -10.6666666666667, -1.33781874467331e-14, -1.8249290967276e-15, 86, -10.6666666666667, -10.75, -10.6666666666667, -4.9960036108132e-15, -1.02140518265514e-14, -10.6666666666667, -21.3333333333333, -10.6666666666667, -10.6666666666667, -10.6666666666667, 85.3333333333333, -3.5527136788005e-15, -1.06581410364015e-14, -10.6666666666667, 8.88178419700125e-16, -10.6666666666667, -10.6666666666667, -2.66453525910038e-15, -10.6666666666667, -10.6666666666667, -3.5527136788005e-15, 42.6666666666667, -10.6666666666667, -3.33066907387547e-15, -21.3333333333333, -10.6666666666667, -1.03146657881581e-14, -10.75, -10.6666666666667, -10.75, -1.06581410364015e-14, -10.6666666666667, 86, -4.48252546192407e-15, -10.6666666666667, 4.37150315946155e-16, -10.6666666666667, -10.6666666666667, -10.75, -3.59434704222394e-15, -10.6666666666667, -3.33066907387547e-15, -4.48252546192407e-15, 43, 6.93889390390723e-18, -10.6666666666667, -10.6666666666667, -6.21724893790088e-15, -10.6666666666667, 8.88178419700125e-16, -10.6666666666667, 42.6666666666667, -6.27276008913213e-15, -10.6666666666667, -10.75, -10.6666666666667, -4.9960036108132e-15, -10.6666666666667, 4.37150315946155e-16, -6.27276008913213e-15, 43, -1.38777878078145e-17, -2.42861286636753e-17, 0.333333333333333, -2.08166817117217e-17, -0.0833333333333333, -1.38777878078145e-17, -2.08166817117217e-17, 0.666666666666667, -3.46944695195361e-17, -0.0833333333333333, -2.08166817117217e-17, -0.0833333333333333, -0.0833333333333333, -3.46944695195361e-17, 1.33333333333333, -5.98479599211998e-17, -0.0833333333333333, -3.16587034365767e-17, -1.04517089427603e-16, -0.0833333333333333, -0.0833333333333333, -1.38777878078145e-17, -0.0833333333333333, -5.98479599211998e-17, 0.666666666666667, -0.0833333333333333, -4.85722573273506e-17, -2.08166817117217e-17, -2.08166817117217e-17, -0.0833333333333333, 0.333333333333333, -2.42861286636753e-17, -0.0833333333333333, -3.16587034365767e-17, -2.42861286636753e-17, 0.666666666666667, -0.0833333333333333, -5.24753851482984e-17, -1.04517089427603e-16, -0.0833333333333333, -0.0833333333333333, 0.666666666666667, -5.29090660172926e-17, -3.90312782094782e-18, 6.93889390390723e-18, -0.0833333333333333, -4.85722573273506e-17, -5.29090660172926e-17, 0.333333333333333, -1.38777878078145e-17, -0.0833333333333333, -5.24753851482984e-17, -3.90312782094782e-18, 0.333333333333333};
    local_matrix_type lclMatrix("lclA", 27, 27, 251, values.data(), rowptr.data(), indices.data());

    droppedA = MatrixFactory::Build(map, map, lclMatrix);
  }

  auto coords = Xpetra::MultiVectorFactory<magnitudeType, LocalOrdinal, GlobalOrdinal, Node>::Build(map, 3);
  {
    magnitudeType data[27][3] = {{0.0, 0.0, 2.0}, {0.0, 0.0, 1.0}, {0.0, 1.0, 1.0}, {0.0, 1.0, 2.0}, {1.0, 0.0, 2.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {1.0, 1.0, 2.0}, {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 2.0, 1.0}, {0.0, 2.0, 2.0}, {1.0, 2.0, 1.0}, {1.0, 2.0, 2.0}, {0.0, 2.0, 0.0}, {1.0, 2.0, 0.0}, {2.0, 0.0, 2.0}, {2.0, 0.0, 1.0}, {2.0, 1.0, 1.0}, {2.0, 1.0, 2.0}, {2.0, 0.0, 0.0}, {2.0, 1.0, 0.0}, {2.0, 2.0, 1.0}, {2.0, 2.0, 2.0}, {2.0, 2.0, 0.0}};
    Kokkos::View<magnitudeType **, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> kv(&data[0][0], 27, 3);
    auto lclMV = coords->getLocalViewHost(Tpetra::Access::OverwriteAll);
    Kokkos::deep_copy(lclMV, kv);
  }

  auto material = MultiVectorFactory::Build(map, 9);
  {
    impl_scalar_type data[27][9] = {{128.0, 0.0, 0.0, 0.0, 128.0, 0.0, 0.0, 0.0, 128.0}, {128.0, 0.0, 0.0, 0.0, 128.0, 0.0, 0.0, 0.0, 128.0}, {128.0, 0.0, 0.0, 0.0, 128.0, 0.0, 0.0, 0.0, 128.0}, {128.0, 0.0, 0.0, 0.0, 128.0, 0.0, 0.0, 0.0, 128.0}, {64.5, 0.0, 0.0, 0.0, 64.5, 0.0, 0.0, 0.0, 64.5}, {64.5, 0.0, 0.0, 0.0, 64.5, 0.0, 0.0, 0.0, 64.5}, {64.5, 0.0, 0.0, 0.0, 64.5, 0.0, 0.0, 0.0, 64.5}, {64.5, 0.0, 0.0, 0.0, 64.5, 0.0, 0.0, 0.0, 64.5}, {128.0, 0.0, 0.0, 0.0, 128.0, 0.0, 0.0, 0.0, 128.0}, {128.0, 0.0, 0.0, 0.0, 128.0, 0.0, 0.0, 0.0, 128.0}, {64.5, 0.0, 0.0, 0.0, 64.5, 0.0, 0.0, 0.0, 64.5}, {64.5, 0.0, 0.0, 0.0, 64.5, 0.0, 0.0, 0.0, 64.5}, {128.0, 0.0, 0.0, 0.0, 128.0, 0.0, 0.0, 0.0, 128.0}, {128.0, 0.0, 0.0, 0.0, 128.0, 0.0, 0.0, 0.0, 128.0}, {64.5, 0.0, 0.0, 0.0, 64.5, 0.0, 0.0, 0.0, 64.5}, {64.5, 0.0, 0.0, 0.0, 64.5, 0.0, 0.0, 0.0, 64.5}, {128.0, 0.0, 0.0, 0.0, 128.0, 0.0, 0.0, 0.0, 128.0}, {64.5, 0.0, 0.0, 0.0, 64.5, 0.0, 0.0, 0.0, 64.5}, {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}, {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}, {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}, {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}, {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}, {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}, {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}, {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}, {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}};
    Kokkos::View<impl_scalar_type **, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> kv(&data[0][0], 27, 9);
    auto lclMV = material->getLocalViewHost(Tpetra::Access::OverwriteAll);
    Kokkos::deep_copy(lclMV, kv);
  }
  return std::make_tuple(A, droppedA, coords, material);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, DistanceLaplacianTensorMaterial, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  RCP<Matrix> A, expectedFilteredA;
  RCP<RealValuedMultiVector> coordinates;
  RCP<MultiVector> material;
  {
    auto testCase     = constructVariableMaterialMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(comm);
    A                 = std::get<0>(testCase);
    expectedFilteredA = std::get<1>(testCase);
    coordinates       = std::get<2>(testCase);
    material          = std::get<3>(testCase);
  }

  RCP<Matrix> filteredA;
  {
    Level fineLevel;
    TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);
    fineLevel.Set("A", A);
    fineLevel.Set("Coordinates", coordinates);
    fineLevel.Set("Material", material);

    CoalesceDropFactory_kokkos coalesceDropFact;
    coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
    coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.02));
    coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
    coalesceDropFact.SetParameter("aggregation: distance laplacian metric", Teuchos::ParameterEntry(std::string("material")));
    coalesceDropFact.SetParameter("filtered matrix: reuse graph", Teuchos::ParameterEntry(false));
    coalesceDropFact.SetParameter("filtered matrix: lumping choice", Teuchos::ParameterEntry(std::string("no lumping")));
    fineLevel.Request("Graph", &coalesceDropFact);
    fineLevel.Request("A", &coalesceDropFact);
    fineLevel.Request("DofsPerNode", &coalesceDropFact);

    coalesceDropFact.Build(fineLevel);

    // RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
    filteredA = fineLevel.Get<RCP<Matrix>>("A", &coalesceDropFact);
  }

  using TF = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>;
  {
    auto lclA                 = A->getLocalMatrixHost();
    auto lclFilteredA         = filteredA->getLocalMatrixHost();
    auto lclExpectedFilteredA = expectedFilteredA->getLocalMatrixHost();

    out << "A:\n"
        << TF::localMatToString(lclA);

    out << "\nFiltered A:\n"
        << TF::localMatToString(lclFilteredA);

    out << "\nExpected A:\n"
        << TF::localMatToString(lclExpectedFilteredA) << std::endl;

    TEST_EQUALITY(lclFilteredA.graph.row_map.extent(0), lclExpectedFilteredA.graph.row_map.extent(0));
    TEST_EQUALITY(lclFilteredA.graph.entries.extent(0), lclExpectedFilteredA.graph.entries.extent(0));

    for (size_t row = 0; row < std::min(lclFilteredA.graph.row_map.extent(0) - 1,
                                        lclExpectedFilteredA.graph.row_map.extent(0) - 1);
         ++row) {
      out << std::endl;
      TEST_EQUALITY(lclFilteredA.graph.row_map(row), lclExpectedFilteredA.graph.row_map(row));
      TEST_EQUALITY(lclFilteredA.graph.row_map(row + 1) - lclFilteredA.graph.row_map(row), lclExpectedFilteredA.graph.row_map(row + 1) - lclExpectedFilteredA.graph.row_map(row));

      for (size_t entry = 0; entry < std::min(lclFilteredA.graph.row_map(row + 1) - lclFilteredA.graph.row_map(row), lclExpectedFilteredA.graph.row_map(row + 1) - lclExpectedFilteredA.graph.row_map(row)); ++entry) {
        auto offset         = lclFilteredA.graph.row_map(row);
        auto offsetExpected = lclExpectedFilteredA.graph.row_map(row);
        TEST_EQUALITY(lclFilteredA.graph.entries(offset + entry), lclExpectedFilteredA.graph.entries(offsetExpected + entry));
        TEST_EQUALITY(lclFilteredA.values(offset + entry), lclExpectedFilteredA.values(offsetExpected + entry));
      }
    }
  }

}  // DistanceLaplacianTensorMaterial

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, DistanceLaplacianSgndRugeStuebenDistribLump, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;
  using TMT = Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  RCP<Matrix> A;
  RCP<RealValuedMultiVector> coordinates;
  {
    Teuchos::ParameterList galeriList;
    galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
    A           = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  }

  RCP<Matrix> filteredA;
  {
    Level fineLevel;
    TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);
    fineLevel.Set("A", A);
    // Doctor coordinates with goal that the filtered matrix drops all entries in lower
    // triangular portion of the matrix (except for final row).
    auto lclCoords = coordinates->getLocalViewHost(Tpetra::Access::OverwriteAll);
    double delta   = 1.1;
    for (size_t i = 0; i < coordinates->getMap()->getLocalNumElements(); i++) lclCoords(i, 0) = pow(delta, coordinates->getMap()->getGlobalElement((LO)i)) / 35.;
    fineLevel.Set("Coordinates", coordinates);

    CoalesceDropFactory_kokkos coalesceDropFact;
    coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
    coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(.95));
    coalesceDropFact.SetParameter("filtered matrix: reuse graph", Teuchos::ParameterEntry(false));
    coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("point-wise")));
    coalesceDropFact.SetParameter("aggregation: strength-of-connection: measure", Teuchos::ParameterEntry(std::string("signed ruge-stueben")));
    coalesceDropFact.SetParameter("aggregation: strength-of-connection: matrix", Teuchos::ParameterEntry(std::string("distance laplacian")));
    coalesceDropFact.SetParameter("filtered matrix: lumping choice", Teuchos::ParameterEntry(std::string("distributed lumping")));
    coalesceDropFact.SetParameter("aggregation: distance laplacian metric", Teuchos::ParameterEntry(std::string("unweighted")));
    fineLevel.Request("Graph", &coalesceDropFact);
    fineLevel.Request("A", &coalesceDropFact);
    fineLevel.Request("DofsPerNode", &coalesceDropFact);

    coalesceDropFact.Build(fineLevel);

    filteredA = fineLevel.Get<RCP<Matrix>>("A", &coalesceDropFact);
  }

  // All interior rows have 2 nonzeros.
  // If we instead chose "no lumping", the diagonal would be 2 and the off-diagonal would be -1
  // If we instead chose "diag lumping", the diagonal would be 1 and the off-diagonal would be -1
  // Since we chose "distributed lumping", the diagonal would be 4/3 and the off-diagonal shoudl be -4/3

  {
    auto lclFilteredA = filteredA->getLocalMatrixHost();

    SC thevalue;
    bool fourThirdsExists    = false;
    bool negFourThirdsExists = false;
    bool hasValidValues      = true;
    // check that all filtered values are either 2, -1, 4/3, -4/3
    // check that at least one row has a 4/3 and that at least one row has -4/3 .
    for (size_t entry = 0; entry < lclFilteredA.graph.entries.extent(0); ++entry) {
      thevalue = lclFilteredA.values(entry);
      if ((Teuchos::ScalarTraits<SC>::magnitude(thevalue - as<Scalar>(4. / 3.)) < 100 * TMT::eps())) fourThirdsExists = true;
      if ((Teuchos::ScalarTraits<SC>::magnitude(thevalue + as<Scalar>(4. / 3.)) < 100 * TMT::eps())) negFourThirdsExists = true;
      if ((Teuchos::ScalarTraits<SC>::magnitude(thevalue - as<Scalar>(2.0)) > 100 * TMT::eps()) &&
          (Teuchos::ScalarTraits<SC>::magnitude(thevalue + as<Scalar>(1.0)) > 100 * TMT::eps()) &&
          (Teuchos::ScalarTraits<SC>::magnitude(thevalue - as<Scalar>(4. / 3.)) > 100 * TMT::eps()) &&
          (Teuchos::ScalarTraits<SC>::magnitude(thevalue + as<Scalar>(4. / 3.)) > 100 * TMT::eps())) hasValidValues = false;
    }
    TEST_EQUALITY(hasValidValues, true);
    // The if's below avoid a test failure just because a processor contains no interior rows
    if (lclFilteredA.graph.entries.extent(0) > 1) TEST_EQUALITY(fourThirdsExists, true);
    if (lclFilteredA.graph.entries.extent(0) > 1) TEST_EQUALITY(negFourThirdsExists, true);
  }

}  // DistanceLaplacianSgndRugeStuebenDistribLump

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, ClassicalScaledCut, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  // Change entry (1,0)
  auto crsA = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A, true)->getCrsMatrix();
  crsA->resumeFill();
  if (comm->getRank() == 0) {
    Teuchos::Array<GlobalOrdinal> cols(3);
    Teuchos::Array<Scalar> vals(3);
    size_t numEntries;
    crsA->getGlobalRowCopy(1, cols, vals, numEntries);
    vals[0] = 0.5;
    crsA->replaceGlobalValues(1, cols, vals);
  }
  crsA->fillComplete();
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("classical")));
  coalesceDropFact.SetParameter("aggregation: classical algo", Teuchos::ParameterEntry(std::string("scaled cut")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph_d = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  auto graph                  = graph_d->copyToHost();
  LO myDofsPerNode            = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 105);

}  // ClassicalScaledCut

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, ClassicalUnScaledCut, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  // Change entry (1,0)
  auto crsA = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A, true)->getCrsMatrix();
  crsA->resumeFill();
  if (comm->getRank() == 0) {
    Teuchos::Array<GlobalOrdinal> cols(3);
    Teuchos::Array<Scalar> vals(3);
    size_t numEntries;
    crsA->getGlobalRowCopy(1, cols, vals, numEntries);
    vals[0] = 0.5;
    crsA->replaceGlobalValues(1, cols, vals);
  }
  crsA->fillComplete();
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);

  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("classical")));
  coalesceDropFact.SetParameter("aggregation: classical algo", Teuchos::ParameterEntry(std::string("unscaled cut")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph_d = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  auto graph                  = graph_d->copyToHost();
  LO myDofsPerNode            = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 105);

}  // ClassicalUnScaledCut

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, ClassicalCutSym, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  // Change entry (1,0)
  auto crsA = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A, true)->getCrsMatrix();
  crsA->resumeFill();
  if (comm->getRank() == 0) {
    Teuchos::Array<GlobalOrdinal> cols(3);
    Teuchos::Array<Scalar> vals(3);
    size_t numEntries;
    crsA->getGlobalRowCopy(1, cols, vals, numEntries);
    vals[0] = 0.5;
    crsA->replaceGlobalValues(1, cols, vals);
  }
  crsA->fillComplete();
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);

  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("classical")));
  coalesceDropFact.SetParameter("aggregation: classical algo", Teuchos::ParameterEntry(std::string("scaled cut symmetric")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph_d = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  auto graph                  = graph_d->copyToHost();
  LO myDofsPerNode            = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 106);

}  // ClassicalCutSym

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, pw_A_SA, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  RCP<ParameterList> params = rcp(new ParameterList());
  params->set("aggregation: drop scheme", "point-wise");
  params->set("aggregation: strength-of-connection: matrix", "A");
  params->set("aggregation: strength-of-connection: measure", "smoothed aggregation");

  GlobalOrdinal nx = 10;
  auto cdth        = CoalesceDropTestingHelper<Scalar, LocalOrdinal, GlobalOrdinal, Node>(nx, params);

  // Make an off-diagonal entry in row 5 positive.
  cdth.replaceEntriesInA({{1, 2, -0.5}});
  cdth.printA(out);

  // We are dropping if
  //   |a_ij| / sqrt(|A_{ii}| |A_{jj}|) < dropTol
  // Entry (1,0) should be dropped for dropTol > 0.25

  cdth.params->set("aggregation: drop tol", 0.24);
  cdth.run();
  cdth.printFilteredA(out);
  TEST_EQUALITY(cdth.dofsPerNode, 1);
  TEST_EQUALITY(cdth.graph->GetGlobalNumEdges(), cdth.A->getGlobalNumEntries());

  cdth.params->set("aggregation: drop tol", 0.26);
  cdth.run();
  cdth.printFilteredA(out);
  TEST_EQUALITY(cdth.dofsPerNode, 1);
  TEST_EQUALITY(cdth.graph->GetGlobalNumEdges(), cdth.A->getGlobalNumEntries() - 1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, pw_A_RS, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  RCP<ParameterList> params = rcp(new ParameterList());
  params->set("aggregation: drop scheme", "point-wise");
  params->set("aggregation: strength-of-connection: matrix", "A");
  params->set("aggregation: strength-of-connection: measure", "signed ruge-stueben");

  GlobalOrdinal nx = 10;
  auto cdth        = CoalesceDropTestingHelper<Scalar, LocalOrdinal, GlobalOrdinal, Node>(nx, params);

  // Make an off-diagonal entry in row 5 positive.
  cdth.replaceEntriesInA({{1, 2, -2.0}});
  cdth.printA(out);

  // We are dropping if
  //   -Re(a_ij) / | max_j -A_{ij}| < dropTol
  // Entry (1,0) should be dropped for dropTol > 0.5

  cdth.params->set("aggregation: drop tol", 0.49);
  cdth.run();
  cdth.printFilteredA(out);
  TEST_EQUALITY(cdth.dofsPerNode, 1);
  TEST_EQUALITY(cdth.graph->GetGlobalNumEdges(), cdth.A->getGlobalNumEntries());

  // dropTol = 0.51
  cdth.params->set("aggregation: drop tol", 0.51);
  cdth.run();
  cdth.printFilteredA(out);
  TEST_EQUALITY(cdth.dofsPerNode, 1);
  TEST_EQUALITY(cdth.graph->GetGlobalNumEdges(), cdth.A->getGlobalNumEntries() - 1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, pw_A_signedSA, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  RCP<ParameterList> params = rcp(new ParameterList());
  params->set("aggregation: drop scheme", "point-wise");
  params->set("aggregation: strength-of-connection: matrix", "A");
  params->set("aggregation: strength-of-connection: measure", "signed smoothed aggregation");

  GlobalOrdinal nx = 10;
  auto cdth        = CoalesceDropTestingHelper<Scalar, LocalOrdinal, GlobalOrdinal, Node>(nx, params);

  // Make an off-diagonal entry in row 5 positive.
  cdth.replaceEntriesInA({{5, 6, 1.0}});
  cdth.printA(out);

  // We are dropping if
  //   -sign(a_ij) |a_ij| / sqrt{|a_ii| |a_jj|} < dropTol

  // In row 5 the positive off-diagonal entry should always be dropped.
  // All other negative off-diagonal entries should be dropped for dropTol > 0.5

  cdth.params->set("aggregation: drop tol", 0.49);
  cdth.run();
  cdth.printFilteredA(out);
  TEST_EQUALITY(cdth.dofsPerNode, 1);
  TEST_EQUALITY(cdth.graph->GetGlobalNumEdges(), cdth.A->getGlobalNumEntries() - 1);

  cdth.params->set("aggregation: drop tol", 0.51);
  cdth.run();
  cdth.printFilteredA(out);
  TEST_EQUALITY(cdth.dofsPerNode, 1);
  TEST_EQUALITY(cdth.graph->GetGlobalNumEdges(), (Xpetra::global_size_t)nx);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, cutdrop_A_SA, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  RCP<ParameterList> params = rcp(new ParameterList());
  params->set("aggregation: drop scheme", "cut-drop");
  params->set("aggregation: strength-of-connection: matrix", "A");
  params->set("aggregation: strength-of-connection: measure", "smoothed aggregation");

  GlobalOrdinal nx = 10;
  auto cdth        = CoalesceDropTestingHelper<Scalar, LocalOrdinal, GlobalOrdinal, Node>(nx, params);

  // Make an off-diagonal entry in row 2 positive and change an entry in row 5.
  cdth.replaceEntriesInA({
      {2, 3, -0.5},
  });
  cdth.printA(out);

  // We are using the SoC
  //   |a_ij| / sqrt(|A_{ii}| |A_{jj}|)

  // In row 2 the unchanged off-diagonal entry should be dropped for dropTol>0.5

  cdth.params->set("aggregation: drop tol", 0.49);
  cdth.run();
  cdth.printFilteredA(out);
  TEST_EQUALITY(cdth.dofsPerNode, 1);
  TEST_EQUALITY(cdth.graph->GetGlobalNumEdges(), cdth.A->getGlobalNumEntries());

  cdth.params->set("aggregation: drop tol", 0.51);
  cdth.run();
  cdth.printFilteredA(out);
  TEST_EQUALITY(cdth.dofsPerNode, 1);
  TEST_EQUALITY(cdth.graph->GetGlobalNumEdges(), cdth.A->getGlobalNumEntries() - 1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, cutdrop_A_RS, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  RCP<ParameterList> params = rcp(new ParameterList());
  params->set("aggregation: drop scheme", "cut-drop");
  params->set("aggregation: strength-of-connection: matrix", "A");
  params->set("aggregation: strength-of-connection: measure", "signed ruge-stueben");

  GlobalOrdinal nx = 10;
  auto cdth        = CoalesceDropTestingHelper<Scalar, LocalOrdinal, GlobalOrdinal, Node>(nx, params);

  // Make an off-diagonal entry in row 2 positive and change an entry in row 5.
  cdth.replaceEntriesInA({{2, 3, 1.0},
                          {5, 6, -1.5}});
  cdth.printA(out);

  // We are using the SoC
  //   -Re(a_ij) / | max_j -A_{ij}|

  // In row 2 the positive off-diagonal entry should always be dropped.
  // In row 7 the unchanged off-diagonal entry should be dropped for dropTol>0.66666

  cdth.params->set("aggregation: drop tol", 0.65);
  cdth.run();
  cdth.printFilteredA(out);
  TEST_EQUALITY(cdth.dofsPerNode, 1);
  TEST_EQUALITY(cdth.graph->GetGlobalNumEdges(), cdth.A->getGlobalNumEntries() - 1);

  cdth.params->set("aggregation: drop tol", 0.67);
  cdth.run();
  cdth.printFilteredA(out);
  TEST_EQUALITY(cdth.dofsPerNode, 1);
  TEST_EQUALITY(cdth.graph->GetGlobalNumEdges(), cdth.A->getGlobalNumEntries() - 2);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, cutdrop_A_signedSA, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  RCP<ParameterList> params = rcp(new ParameterList());
  params->set("aggregation: drop scheme", "cut-drop");
  params->set("aggregation: strength-of-connection: matrix", "A");
  params->set("aggregation: strength-of-connection: measure", "signed smoothed aggregation");

  GlobalOrdinal nx = 10;
  auto cdth        = CoalesceDropTestingHelper<Scalar, LocalOrdinal, GlobalOrdinal, Node>(nx, params);

  // Make an off-diagonal entry in row 2 positive and change an entry in row 5.
  cdth.replaceEntriesInA({{2, 3, 1.0},
                          {5, 6, -1.5}});
  cdth.printA(out);

  // We are using the SoC
  //   -sign(a_ij) |a_ij| / sqrt{|a_ii| |a_jj|}

  // In row 2 the positive off-diagonal entry should always be dropped.
  // In row 7 the unchanged off-diagonal entry should be dropped for dropTol>0.66666

  cdth.params->set("aggregation: drop tol", 0.66);
  cdth.run();
  cdth.printFilteredA(out);
  TEST_EQUALITY(cdth.dofsPerNode, 1);
  TEST_EQUALITY(cdth.graph->GetGlobalNumEdges(), cdth.A->getGlobalNumEntries() - 1);

  cdth.params->set("aggregation: drop tol", 0.67);
  cdth.run();
  cdth.printFilteredA(out);
  TEST_EQUALITY(cdth.dofsPerNode, 1);
  TEST_EQUALITY(cdth.graph->GetGlobalNumEdges(), cdth.A->getGlobalNumEntries() - 2);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, SignedScaledCutClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("signed classical")));
  coalesceDropFact.SetParameter("aggregation: classical algo", Teuchos::ParameterEntry(std::string("scaled cut")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  RCP<LWGraph_kokkos> graph_d = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  auto graph                  = graph_d->copyToHost();
  LO myDofsPerNode            = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 106);

}  // SignedScaledCutClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, SignedUnscaledCutClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("signed classical")));
  coalesceDropFact.SetParameter("aggregation: classical algo", Teuchos::ParameterEntry(std::string("unscaled cut")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  RCP<LWGraph_kokkos> graph_d = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  auto graph                  = graph_d->copyToHost();
  LO myDofsPerNode            = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 106);

}  // SignedUnScaledCutClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, BlockDiagonalColoredSignedClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0 / 0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal colored signed classical")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph_d = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  auto graph                  = graph_d->copyToHost();
  LO myDofsPerNode            = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 36);

}  // BlockDiagonalColoredSignedClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, BlockDiagonalNoColoredSignedClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  // this test is only compatible with rank higher than 1
  if (comm->getSize() == 1) {
    return;
  }

  // Default is Laplace1D with nx = 8748.
  // It's a nice size for 1D and perfect aggregation. (6561 = 3^8)
  // Nice size for 1D and perfect aggregation on small numbers of processors. (8748 = 4*3^7)
  Teuchos::CommandLineProcessor clp(false);
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748);  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);

  RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  //    getCrsGraph()->getImporter()
  RCP<const Import> importer = ImportFactory::Build(A->getRowMap(), map);
  fineLevel.Set("Importer", importer);
  auto importerTest = A->getCrsGraph()->getImporter();  // NULL
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0 / 0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal colored signed classical")));
  coalesceDropFact.SetParameter("aggregation: coloring: localize color graph", Teuchos::ParameterEntry(false));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  // Need an importer
  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph_d = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  auto graph                  = graph_d->copyToHost();
  LO myDofsPerNode            = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 36);

}  // BlockDiagonalNoColoredSignedClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, BlockDiagonalSignedClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0 / 0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal signed classical")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph_d = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  auto graph                  = graph_d->copyToHost();
  LO myDofsPerNode            = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 36);

}  // BlockDiagonalSignedClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, BlockDiagonal, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;
  typedef Teuchos::ScalarTraits<LO> STL;
  LO zero = STL::zero(), one = STL::one();

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  // interleaved block diagonal
  {
    GO nx = 10 * comm->getSize();
    Teuchos::ParameterList matrixList;
    matrixList.set("nx", nx);
    RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(matrixList, lib);
    A->SetFixedBlockSize(1);  // So we can block diagonalize
    Level fineLevel;
    fineLevel.Set("A", A);

    RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), matrixList);
    fineLevel.Set("Coordinates", coordinates);

    RCP<InitialBlockNumberFactory> ibFact = rcp(new InitialBlockNumberFactory());
    Teuchos::ParameterList ibList;
    ibList.set("aggregation: block diagonal: interleaved blocksize", 3);
    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
    CoalesceDropFactory_kokkos coalesceDropFact;
    coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
    coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    coalesceDropFact.SetFactory("BlockNumber", ibFact);
    coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0 / 8.0));
    coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal")));
    coalesceDropFact.SetParameter("aggregation: block diagonal: interleaved blocksize", Teuchos::ParameterEntry(3));
    coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);

    fineLevel.Request("Graph", &coalesceDropFact);
    fineLevel.Request("DofsPerNode", &coalesceDropFact);

    coalesceDropFact.Build(fineLevel);

    RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
    LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
    TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
  }

  // block diagonal
  {
    GO n = 10 * comm->getSize();

    RCP<Matrix> A             = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build2DPoisson(n);
    RCP<LOVector> blocknumber = Xpetra::VectorFactory<LO, LO, GO, NO>::Build(A->getDomainMap());

    for (size_t row = 0; row < A->getLocalNumRows(); row++) {
      GO global_row = A->getRowMap()->getGlobalElement(row);

      if (global_row < 0.5 * n * n) {
        // lower part of domain get's 0
        blocknumber->replaceLocalValue(row, zero);
      } else {
        // upper part of domain get's 1
        blocknumber->replaceLocalValue(row, one);
      }
    }

    Level fineLevel;
    fineLevel.Set("A", A);
    fineLevel.Set("BlockNumber", blocknumber);

    Teuchos::ParameterList ibList;
    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
    CoalesceDropFactory_kokkos coalesceDropFact;
    coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal")));
    coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.0001));
    coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);

    fineLevel.Request("Graph", &coalesceDropFact);

    coalesceDropFact.Build(fineLevel);

    RCP<LWGraph_kokkos> graph_d = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
    auto graph                  = graph_d->copyToHost();
    auto numGlobalEdges         = graph->GetGlobalNumEdges();
    auto fullNumGlobalEdges     = A->getCrsGraph()->getGlobalNumEntries();

    // we drop exactly two off-diagonal blocks each of entry size n
    TEST_EQUALITY(Teuchos::as<GO>(fullNumGlobalEdges - numGlobalEdges) == 2 * n, true);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, BlockDiagonalVector, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;
  typedef Teuchos::ScalarTraits<LO> STL;
  LO zero = STL::zero(), one = STL::one();

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  // block diagonal - dof per node > 1
  {
    auto runAndCheck = [&](const RCP<Matrix> &A, const RCP<LOVector> &blocknumber, LO ndofn, const std::string &algo, const RCP<RealValuedMultiVector> &coords = Teuchos::null) {
      Level fineLevel;
      fineLevel.Set("A", A);
      fineLevel.Set("BlockNumber", blocknumber);

      RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
      CoalesceDropFactory_kokkos coalesceDropFact;
      coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
      coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(algo));
      coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.0001));
      if (coords != Teuchos::null)
        fineLevel.Set("Coordinates", coords);

      coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);

      fineLevel.Request("Graph", &coalesceDropFact);
      fineLevel.Request("A", &coalesceDropFact);
      fineLevel.Request("DofsPerNode", &coalesceDropFact);

      fineLevel.Request("UnAmalgamationInfo", amalgFact.get());
      RCP<AmalgamationInfo> unAmalgInfo = fineLevel.Get<RCP<AmalgamationInfo>>("UnAmalgamationInfo", amalgFact.get());

      LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
      TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == ndofn, true);

      RCP<LWGraph_kokkos> graph_d = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);

      auto graph_h  = graph_d->copyToHost();
      auto crsGraph = graph_h->GetCrsGraph();

      RCP<LOVector> blocknumberGhosted;
      const RCP<const Map> uniqueMap    = unAmalgInfo->getNodeRowMap();
      const RCP<const Map> nonUniqueMap = unAmalgInfo->getNodeColMap();
      auto importer                     = Xpetra::ImportFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(blocknumber->getMap(), nonUniqueMap);
      if (!importer.is_null()) {
        blocknumberGhosted = Xpetra::VectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>::Build(importer->getTargetMap());
        blocknumberGhosted->doImport(*blocknumber, *importer, Xpetra::INSERT);
      } else {
        blocknumberGhosted = blocknumber;
      }

      auto blockNums        = blocknumber->getData(0);
      auto blockNumsGhosted = blocknumberGhosted->getData(0);

      // test if there are any edges between nodes which do not belong to the same block
      bool testEmptyBlockDiagonals = true;
      for (LO i = 0; i < Teuchos::as<LO>(crsGraph->getLocalNumRows()); ++i) {
        Teuchos::ArrayView<const LO> rowEntries;
        crsGraph->getLocalRowView(i, rowEntries);

        for (LO k = 0; k < rowEntries.size(); ++k) {
          LO j = rowEntries[k];

          if (blockNums[i] != blockNumsGhosted[j])
            testEmptyBlockDiagonals = false;
        }
      }
      TEST_EQUALITY(testEmptyBlockDiagonals, true);
    };

    // Elasticity3D, contiguous blocknumber
    {
      constexpr LO ndofn = 3;

      // unrelated to 3 dofs per node; just the smallest matrix that isn't nonzero everywhere in serial
      GO n = 3 * comm->getSize();

      Teuchos::ParameterList galeriList;
      galeriList.set("nx", Teuchos::as<GlobalOrdinal>(n));
      galeriList.set("ny", Teuchos::as<GlobalOrdinal>(n));
      galeriList.set("nz", Teuchos::as<GlobalOrdinal>(n));

      // node map
      auto map = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian3D", comm, galeriList);

      RCP<LOVector> blocknumber = Xpetra::VectorFactory<LO, LO, GO, NO>::Build(map);
      // domain is split into 2 blocks (upper and lower)
      for (size_t row = 0; row < blocknumber->getLocalLength(); row++) {
        GO global_row = map->getGlobalElement(row);

        if (global_row < 0.5 * n * n * n) {
          blocknumber->replaceLocalValue(row, zero);
        } else {
          blocknumber->replaceLocalValue(row, one);
        }
      }

      // dof map
      map = Xpetra::MapFactory<LO, GO, Node>::Build(map, ndofn);

      RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector>> Pr =
          Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Elasticity3D", map, galeriList);

      RCP<Matrix> A = Pr->BuildMatrix();
      A->SetFixedBlockSize(ndofn);

      // run and check/test all algorithms involving block diagonalization
      runAndCheck(A, blocknumber, ndofn, "block diagonal");
      runAndCheck(A, blocknumber, ndofn, "block diagonal signed classical");
      runAndCheck(A, blocknumber, ndofn, "block diagonal classical");
      runAndCheck(A, blocknumber, ndofn, "block diagonal colored signed classical");

      auto coordinates                          = Galeri::Xpetra::Utils::CreateCartesianCoordinates<real_type, LO, GO, Map, RealValuedMultiVector>("3D", map, galeriList);
      RCP<RealValuedMultiVector> newcoordinates = Pr->BuildCoords();

      for (size_t kkk = 0; kkk < coordinates->getNumVectors(); kkk++) {
        Teuchos::ArrayRCP<real_type> old     = coordinates->getDataNonConst(kkk);
        Teuchos::ArrayRCP<real_type> newvals = newcoordinates->getDataNonConst(kkk);
        int numCopies                        = newvals.size() / old.size();
        for (int jj = 0; jj < old.size(); jj++) old[jj] = newvals[numCopies * (Teuchos_Ordinal)jj];
      }
      runAndCheck(A, blocknumber, ndofn, "block diagonal distance laplacian", coordinates);
    }

    // Elasticity2D, non-contiguous blocknumber
    {
      constexpr LO ndofn = 2;

      GO n = 5 * comm->getSize();

      Teuchos::ParameterList galeriList;
      galeriList.set("nx", Teuchos::as<GlobalOrdinal>(n));
      galeriList.set("ny", Teuchos::as<GlobalOrdinal>(n));

      // node map
      auto map = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian2D", comm, galeriList);

      RCP<LOVector> blocknumber = Xpetra::VectorFactory<LO, LO, GO, NO>::Build(map);
      // 3 non-contiguous blocks
      for (size_t row = 0; row < blocknumber->getLocalLength(); row++) {
        GO global_row = map->getGlobalElement(row);
        blocknumber->replaceLocalValue(row, Teuchos::as<LO>(global_row / 4 % 3));
      }

      // dof map
      map = Xpetra::MapFactory<LO, GO, Node>::Build(map, ndofn);

      RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector>> Pr =
          Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Elasticity2D", map, galeriList);
      RCP<Matrix> A = Pr->BuildMatrix();
      A->SetFixedBlockSize(ndofn);

      // run and check/test all algorithms involving block diagonalization
      runAndCheck(A, blocknumber, ndofn, "block diagonal");
      // runAndCheck(A, blocknumber, ndofn, "block diagonal signed classical");  //CAG: This bugs out on >1 threads with OpenMP. No idea why.
      runAndCheck(A, blocknumber, ndofn, "block diagonal classical");
      // runAndCheck(A, blocknumber, ndofn, "block diagonal colored signed classical");  //CAG: This bugs out on >1 threads with OpenMP. No idea why.

      auto coordinates                          = Galeri::Xpetra::Utils::CreateCartesianCoordinates<real_type, LO, GO, Map, RealValuedMultiVector>("2D", map, galeriList);
      RCP<RealValuedMultiVector> newcoordinates = Pr->BuildCoords();

      for (size_t kkk = 0; kkk < coordinates->getNumVectors(); kkk++) {
        Teuchos::ArrayRCP<real_type> old     = coordinates->getDataNonConst(kkk);
        Teuchos::ArrayRCP<real_type> newvals = newcoordinates->getDataNonConst(kkk);
        int numCopies                        = newvals.size() / old.size();
        for (int jj = 0; jj < old.size(); jj++) old[jj] = newvals[numCopies * (Teuchos_Ordinal)jj];
      }
      runAndCheck(A, blocknumber, ndofn, "block diagonal distance laplacian", coordinates);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, BlockDiagonalClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(matrixList, lib);
  A->SetFixedBlockSize(1);  // So we can block diagonalize
  Level fineLevel;
  fineLevel.Set("A", A);

  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), matrixList);
  fineLevel.Set("Coordinates", coordinates);

  RCP<InitialBlockNumberFactory> ibFact = rcp(new InitialBlockNumberFactory());
  Teuchos::ParameterList ibList;
  ibList.set("aggregation: block diagonal: interleaved blocksize", 3);
  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetFactory("BlockNumber", ibFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0 / 8.0));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal classical")));
  coalesceDropFact.SetParameter("aggregation: block diagonal: interleaved blocksize", Teuchos::ParameterEntry(3));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, BlockDiagonalDistanceLaplacian, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(matrixList, lib);
  A->SetFixedBlockSize(1);  // So we can block diagonalize
  Level fineLevel;
  fineLevel.Set("A", A);

  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), matrixList);
  fineLevel.Set("Coordinates", coordinates);

  RCP<InitialBlockNumberFactory> ibFact = rcp(new InitialBlockNumberFactory());
  Teuchos::ParameterList ibList;
  ibList.set("aggregation: block diagonal: interleaved blocksize", 3);
  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetFactory("BlockNumber", ibFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.025));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: block diagonal: interleaved blocksize", Teuchos::ParameterEntry(3));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, BlockDiagonalDistanceDifferentCoordinatesLaplacian, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(matrixList, lib);
  A->SetFixedBlockSize(1);  // So we can block diagonalize
  Level fineLevel;
  fineLevel.Set("A", A);

  GO bnx = 15 * comm->getSize();
  Teuchos::ParameterList bMatrixList;
  matrixList.set("bnx", bnx);
  RCP<Matrix> B = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(bMatrixList, lib);

  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", B->getRowMap(), matrixList);
  fineLevel.Set("Coordinates", coordinates);

  RCP<InitialBlockNumberFactory> ibFact = rcp(new InitialBlockNumberFactory());
  Teuchos::ParameterList ibList;
  ibList.set("aggregation: block diagonal: interleaved blocksize", 3);
  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetFactory("BlockNumber", ibFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.025));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: block diagonal: interleaved blocksize", Teuchos::ParameterEntry(3));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, BlockDiagonalDistanceLaplacianWeighted, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  matrixList.set("ny", (GO)10);
  matrixList.set("nz", (GO)10);
  matrixList.set("matrixType", "Laplace3D");
  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(matrixList, lib);
  A->SetFixedBlockSize(1);  // So we can block diagonalize
  Level fineLevel;
  fineLevel.Set("A", A);

  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("3D", A->getRowMap(), matrixList);
  fineLevel.Set("Coordinates", coordinates);

  RCP<InitialBlockNumberFactory> ibFact = rcp(new InitialBlockNumberFactory());
  Teuchos::ParameterList ibList;
  ibList.set("aggregation: block diagonal: interleaved blocksize", 3);
  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetFactory("BlockNumber", ibFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.025));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: block diagonal: interleaved blocksize", Teuchos::ParameterEntry(3));
  std::vector<double> weights_v{100.0, 1.0, 1.0, 1.0, 100, 1.0, 1.0, 1.0, 100.0};
  Teuchos::Array<double> weights(weights_v);
  coalesceDropFact.SetParameter("aggregation: distance laplacian directional weights", Teuchos::ParameterEntry(weights));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, DistanceLaplacianWeighted, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  matrixList.set("ny", (GO)10);
  matrixList.set("nz", (GO)10);
  matrixList.set("matrixType", "Laplace3D");
  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildMatrix(matrixList, lib);

  Level fineLevel;
  fineLevel.Set("A", A);

  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("3D", A->getRowMap(), matrixList);
  fineLevel.Set("Coordinates", coordinates);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.025));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
  std::vector<double> weights_v{100.0, 1.0, 1.0};
  Teuchos::Array<double> weights(weights_v);
  coalesceDropFact.SetParameter("aggregation: distance laplacian directional weights", Teuchos::ParameterEntry(weights));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, ClassicScalarWithoutFiltering, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory);
  CoalesceDropFactory_kokkos dropFact;
  dropFact.SetFactory("UnAmalgamationInfo", amalgFact);

  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  auto graph         = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &dropFact);
  auto myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);

  TEST_EQUALITY(as<int>(myDofsPerNode) == 1, true);

  bool bCorrectGraph = false;
  int reduction_val  = 0;
  int comm_size = comm->getSize(), comm_rank = comm->getRank();
  auto lclLWGraph = *graph;
  Kokkos::parallel_reduce(
      "MueLu:TentativePF:Build:compute_agg_sizes", Kokkos::RangePolicy<typename NO::execution_space, size_t>(0, 1),
      KOKKOS_LAMBDA(const LO i, int &correct) {
        if (comm_size == 1) {
          auto v0 = lclLWGraph.getNeighborVertices(0);
          auto v1 = lclLWGraph.getNeighborVertices(1);
          auto v2 = lclLWGraph.getNeighborVertices(2);
          if (v0.length == 2 && ((v0(0) == 0 && v0(1) == 1) || (v0(0) == 1 && v0(1) == 0)) &&
              v1.length == 3 && v2.length == 3)
            correct = true;
        } else {
          if (comm_rank == 0) {
            if (lclLWGraph.getNeighborVertices(0).length == 2)
              correct = true;

          } else {
            if (lclLWGraph.getNeighborVertices(0).length == 3)
              correct = true;
          }
        }
      },
      reduction_val);
  bCorrectGraph = reduction_val;
  TEST_EQUALITY(bCorrectGraph, true);

  auto myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  auto myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, ClassicScalarWithFiltering, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  auto dofMap = MapFactory::Build(lib, 3 * comm->getSize(), 0, comm);
  auto mtx    = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 1.0, -1.0, -0.0001);

  mtx->SetFixedBlockSize(1, 0);
  fineLevel.Set("A", mtx);

  RCP<AmalgamationFactory> amalgFact  = rcp(new AmalgamationFactory);
  CoalesceDropFactory_kokkos dropFact = CoalesceDropFactory_kokkos();
  dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  dropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.5));

  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  auto graph         = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &dropFact);
  auto myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);

  TEST_EQUALITY(as<int>(myDofsPerNode) == 1, true);
  TEST_EQUALITY(as<int>(graph->GetDomainMap()->getGlobalNumElements()) == 3 * comm->getSize(), true);

  bool bCorrectGraph = false;
  int reduction_val  = 0;
  int comm_size = comm->getSize(), comm_rank = comm->getRank();
  auto lclLWGraph = *graph;
  Kokkos::parallel_reduce(
      "MueLu:TentativePF:Build:compute_agg_sizes", Kokkos::RangePolicy<typename NO::execution_space, size_t>(0, 1),
      KOKKOS_LAMBDA(const LO i, int &correct) {
        if (comm_size == 1) {
          auto v0 = lclLWGraph.getNeighborVertices(0);
          auto v1 = lclLWGraph.getNeighborVertices(1);
          auto v2 = lclLWGraph.getNeighborVertices(2);
          if (v0.length == 1 && v0(0) == 0 &&
              v1.length == 2 && ((v1(0) == 0 && v1(1) == 1) || (v1(0) == 1 && v1(1) == 0)) &&
              v2.length == 2 && ((v2(0) == 1 && v2(1) == 2) || (v2(0) == 2 && v2(1) == 1)))
            correct = true;
        } else {
          if (comm_rank == 0) {
            if (lclLWGraph.getNeighborVertices(0).length == 1)
              correct = true;

          } else {
            if (lclLWGraph.getNeighborVertices(0).length == 2)
              correct = true;
          }
        }
      },
      reduction_val);
  bCorrectGraph = reduction_val;
  TEST_EQUALITY(bCorrectGraph, true);

  auto myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  auto myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 3 * comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), as<size_t>(3 * comm->getSize() + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 3 * comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), as<size_t>(3 * comm->getSize()));
  TEST_EQUALITY(myDomainMap->getLocalNumElements(), 3);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, ClassicBlockWithoutFiltering, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  int blockSize = 3;

  auto dofMap = MapFactory::Build(lib, blockSize * comm->getSize(), 0, comm);
  auto mtx    = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 2.0, -1.0, -1.0);
  mtx->SetFixedBlockSize(blockSize, 0);
  fineLevel.Set("A", mtx);

  RCP<AmalgamationFactory> amalgFact  = rcp(new AmalgamationFactory);
  CoalesceDropFactory_kokkos dropFact = CoalesceDropFactory_kokkos();
  dropFact.SetFactory("UnAmalgamationInfo", amalgFact);

  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  auto graph         = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &dropFact);
  auto myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);

  TEST_EQUALITY(as<int>(myDofsPerNode) == blockSize, true);
  TEST_EQUALITY(as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);

  bool bCorrectGraph = false;
  int reduction_val  = 0;
  int comm_size = comm->getSize(), comm_rank = comm->getRank();
  auto lclLWGraph = *graph;
  Kokkos::parallel_reduce(
      "MueLu:TentativePF:Build:compute_agg_sizes", Kokkos::RangePolicy<typename NO::execution_space, size_t>(0, 1),
      KOKKOS_LAMBDA(const LO i, int &correct) {
        if (comm_size == 1 && lclLWGraph.getNeighborVertices(0).length == 1) {
          correct = true;
        } else {
          if (comm_rank == 0 || comm_rank == comm_size - 1) {
            if (lclLWGraph.getNeighborVertices(0).length == 2)
              correct = true;

          } else {
            if (static_cast<int>(lclLWGraph.getNeighborVertices(0).length) == blockSize)
              correct = true;
          }
        }
      },
      reduction_val);
  bCorrectGraph = reduction_val;
  TEST_EQUALITY(bCorrectGraph, true);

  auto myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  auto myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), as<size_t>(comm->getSize() + 2 * (comm->getSize() - 1)));

  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t numLocalImportElts = myImportMap->getLocalNumElements();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(numLocalImportElts, numLocalRowMapElts + 1);
    } else {
      TEST_EQUALITY(numLocalImportElts, numLocalRowMapElts + 2);
    }
  }
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t maxLocalIndex      = myImportMap->getMaxLocalIndex();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(maxLocalIndex, numLocalRowMapElts * blockSize - 2);
    } else {
      TEST_EQUALITY(maxLocalIndex, numLocalRowMapElts * blockSize - 1);
    }
  }

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMaxLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), as<size_t>(comm->getSize()));
  TEST_EQUALITY(myDomainMap->getLocalNumElements(), 1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, ClassicBlockWithFiltering, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  auto dofMap = MapFactory::Build(lib, 3 * comm->getSize(), 0, comm);
  auto mtx    = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 2.0, -1.0, 0.00001);

  mtx->SetFixedBlockSize(3, 0);
  fineLevel.Set("A", mtx);

  RCP<AmalgamationFactory> amalgFact  = rcp(new AmalgamationFactory);
  CoalesceDropFactory_kokkos dropFact = CoalesceDropFactory_kokkos();
  dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  dropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0));

  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  auto graph         = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &dropFact);
  auto myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);

  TEST_EQUALITY(as<int>(myDofsPerNode) == 3, true);
  TEST_EQUALITY(as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);

  TEST_EQUALITY(graph->getNeighborVertices(0).size(), 1);

  auto myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  auto myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), as<size_t>(comm->getSize() + 2 * (comm->getSize() - 1)));
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t numLocalImportElts = myImportMap->getLocalNumElements();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(numLocalImportElts, numLocalRowMapElts + 1);
    } else {
      TEST_EQUALITY(numLocalImportElts, numLocalRowMapElts + 2);
    }
  }

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMaxLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), as<size_t>(comm->getSize()));
  TEST_EQUALITY(myDomainMap->getLocalNumElements(), 1);
}

#if 0
  TEUCHOS_UNIT_TEST(CoalesceDropFactory_kokkos, LaplacianScalarWithoutFiltering)
  {
#include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    Level fineLevel;
    TestHelpers_kokkos::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC,LO,GO,NO>::Build1DPoisson(36);
    fineLevel.Set("A", A);

    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory);
    CoalesceDropFactory_kokkos dropFact;
    dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    dropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));

    fineLevel.Request("Graph",       &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    auto graph_d         = fineLevel.Get<RCP<LWGraph_kokkos> >("Graph",       &dropFact);
    auto graph = graph_d->copyToHost();
    auto myDofsPerNode = fineLevel.Get<LO>                  ("DofsPerNode", &dropFact);
    TEST_EQUALITY(as<int>(myDofsPerNode) == 1, true);

    bool bCorrectGraph = false;
    if (comm->getSize() == 1) {
      auto v0 = graph->getNeighborVertices(0);
      auto v1 = graph->getNeighborVertices(1);
      auto v2 = graph->getNeighborVertices(2);
      if (v0.size() == 2 && ((v0(0) == 0 && v0(1) == 1) || (v0(0) == 1 && v0(1) == 0)) &&
          v1.size() == 3 && v2.size() == 3)
        bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 ) {
        if (graph->getNeighborVertices(0).size() == 2)
          bCorrectGraph = true;

      } else {
        if (graph->getNeighborVertices(0).size() == 3)
          bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    auto myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    auto myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(),  35);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(),  0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),      0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),  as<size_t>(36 + (comm->getSize()-1)*2));

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(),  35);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(),  0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),      0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),  36);
  }
#endif

#if 0
  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationStridedLW)
  {
#include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    // unit test for block size 3 using a strided map

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    int blockSize=3;

    GO nx = blockSize*comm->getSize();
    RCP<Matrix> A = TestHelpers::TestFactory<SC,LO,GO,NO>::Build1DPoisson(nx);

    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(as<size_t>(blockSize));
    LocalOrdinal stridedBlockId = -1;

    RCP<const Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node> > stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                                  A->getRangeMap(),
                                                  stridingInfo,
                                                  stridedBlockId,
                                                  0 /*offset*/
                                                  );
    RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                            A->getDomainMap(),
                                            stridingInfo,
                                            stridedBlockId,
                                            0 /*offset*/
                                            );

    if(A->IsView("stridedMaps") == true) A->RemoveView("stridedMaps");
    A->CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);

    fineLevel.Set("A", A);
    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory);
    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);
    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(as<int>(myDofsPerNode) == blockSize, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        if (graph->getNeighborVertices(0).size() == 2) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == blockSize) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    if (comm->getSize()>1) {
      size_t numLocalRowMapElts = graph->GetNodeNumVertices();
      size_t numLocalImportElts = myImportMap->getLocalNumElements();
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        TEST_EQUALITY(as<bool>(numLocalImportElts==numLocalRowMapElts+1), true);
      } else {
        TEST_EQUALITY(as<bool>(numLocalImportElts==numLocalRowMapElts+2), true);
      }
    }
    if (comm->getSize()>1) {
      size_t numLocalRowMapElts = graph->GetNodeNumVertices();
      size_t maxLocalIndex = myImportMap->getMaxLocalIndex();
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        TEST_EQUALITY(as<bool>(maxLocalIndex==numLocalRowMapElts*blockSize-2), true);
      } else {
        TEST_EQUALITY(as<bool>(maxLocalIndex==numLocalRowMapElts*blockSize-1), true);
      }
    }

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getMaxLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),as<size_t>(comm->getSize()));
    TEST_EQUALITY(as<bool>(myDomainMap->getLocalNumElements()==1), true);
  } // AmalgamationStridedLW

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationStrided2LW)
  {
#include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    // unit test for block size 3 = (2,1). wrap block 0

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    // create strided map information
    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(as<size_t>(2));
    stridingInfo.push_back(as<size_t>(1));
    LocalOrdinal stridedBlockId = 0;

    int blockSize=3;

    RCP<const StridedMap> dofMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, blockSize*comm->getSize(), 0,
                                  stridingInfo, comm,
                                  stridedBlockId /*blockId*/, 0 /*offset*/);

    /////////////////////////////////////////////////////

    Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC,LO,GO,NO>::BuildTridiag(dofMap, 2.0, -1.0, -1.0);

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    RCP<const Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node> > stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                                  mtx->getRangeMap(),
                                                  stridingInfo,
                                                  stridedBlockId,
                                                  0 /*offset*/
                                                  );
    RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                            mtx->getDomainMap(),
                                            stridingInfo,
                                            stridedBlockId,
                                            0 /*offset*/
                                            );
    if(mtx->IsView("stridedMaps") == true) mtx->RemoveView("stridedMaps");
    mtx->CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);

    fineLevel.Set("A", mtx);
    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory);
    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);

    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(as<int>(myDofsPerNode) == blockSize, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        if (graph->getNeighborVertices(0).size() == 2) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == blockSize) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    if (comm->getSize()>1) {
      size_t numLocalRowMapElts = graph->GetNodeNumVertices();
      size_t numLocalImportElts = myImportMap->getLocalNumElements();
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        TEST_EQUALITY(as<bool>(numLocalImportElts==numLocalRowMapElts+1), true);
      } else {
        TEST_EQUALITY(as<bool>(numLocalImportElts==numLocalRowMapElts+2), true);
      }
    }
    if (comm->getSize()>1) {
      size_t numLocalRowMapElts = graph->GetNodeNumVertices();
      size_t maxLocalIndex = myImportMap->getMaxLocalIndex();
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        TEST_EQUALITY(as<bool>(maxLocalIndex==numLocalRowMapElts*blockSize-2), true);
      } else {
        TEST_EQUALITY(as<bool>(maxLocalIndex==numLocalRowMapElts*blockSize-1), true);
      }
    }

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getMaxLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),as<size_t>(comm->getSize()));
    TEST_EQUALITY(as<bool>(myDomainMap->getLocalNumElements()==1), true);
  } // AmalgamationStrided2LW


  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationStridedOffsetDropping2LW)
  {
    // unit test for block size 9 = (2,3,4). wrap block 1.
    // drop small entries
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    // create strided map information
    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(as<size_t>(2));
    stridingInfo.push_back(as<size_t>(3));
    stridingInfo.push_back(as<size_t>(4));
    LocalOrdinal stridedBlockId = 1;
    GlobalOrdinal offset = 19;

    RCP<const StridedMap> dofMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 9*comm->getSize(), 0,
                                  stridingInfo, comm,
                                  stridedBlockId, offset);

    /////////////////////////////////////////////////////

    Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC,LO,GO,NO>::BuildTridiag(dofMap, 2.0, 1.0, 0.0001);

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    RCP<const Map> stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                                  mtx->getRangeMap(),
                                                  stridingInfo,
                                                  stridedBlockId,
                                                  offset
                                                  );
    RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                            mtx->getDomainMap(),
                                            stridingInfo,
                                            stridedBlockId,
                                            offset
                                            );

    if(mtx->IsView("stridedMaps") == true) mtx->RemoveView("stridedMaps");
    mtx->CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);

    fineLevel.Set("A", mtx);
    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory);
    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    dropFact.SetParameter("aggregation: drop tol",Teuchos::ParameterEntry(0.3));

    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);

    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(as<int>(myDofsPerNode) == 9, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0) {
        if (graph->getNeighborVertices(0).size() == 1) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == 2) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    if (comm->getSize()>1) {
      size_t numLocalRowMapElts = graph->GetNodeNumVertices();
      size_t numLocalImportElts = myImportMap->getLocalNumElements();
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        TEST_EQUALITY(as<bool>(numLocalImportElts==numLocalRowMapElts+1), true);
      } else {
        TEST_EQUALITY(as<bool>(numLocalImportElts==numLocalRowMapElts+2), true);
      }
    }
    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),as<size_t>(comm->getSize()));
    TEST_EQUALITY(as<bool>(myDomainMap->getLocalNumElements()==1), true);
  } // AmalgamationStridedOffsetDropping2LW
#endif

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, AggresiveDroppingIsMarkedAsBoundary, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
// Test that when everything but the diagonal is dropped, the node is marked as boundary
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  RCP<const Map> dofMap    = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 12 * comm->getSize(), 0, comm);
  Teuchos::RCP<Matrix> mtx = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 2.0, -1.0, -1.0);

  {
    Level fineLevel;
    TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

    mtx->SetFixedBlockSize(1);
    fineLevel.Set("A", mtx);

    CoalesceDropFactory_kokkos dropFact = CoalesceDropFactory_kokkos();
    RCP<AmalgamationFactory> amalgFact  = rcp(new AmalgamationFactory());
    dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    dropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(4.1));
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &dropFact);

    auto boundaryNodes     = graph->GetBoundaryNodeMap();
    auto boundaryNodesHost = Kokkos::create_mirror_view(boundaryNodes);
    Kokkos::deep_copy(boundaryNodesHost, boundaryNodes);
    bool allNodesAreOnBoundary = true;
    for (LO i = 0; i < Teuchos::as<LO>(boundaryNodesHost.size()); i++)
      allNodesAreOnBoundary &= boundaryNodesHost(i);
    TEST_EQUALITY(allNodesAreOnBoundary, true);
  }

  // {
  //   Level fineLevel;
  //   TestHelpers_kokkos::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

  //   mtx->SetFixedBlockSize(2);
  //   fineLevel.Set("A", mtx);

  //   CoalesceDropFactory_kokkos dropFact = CoalesceDropFactory_kokkos();
  //   RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  //   dropFact.SetFactory("UnAmalgamationInfo",amalgFact);
  //   dropFact.SetParameter("aggregation: drop tol",Teuchos::ParameterEntry(4.1));
  //   fineLevel.Request("Graph", &dropFact);
  //   fineLevel.Request("DofsPerNode", &dropFact);

  //   dropFact.Build(fineLevel);

  //   RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos> >("Graph", &dropFact);

  //   auto boundaryNodes = graph->GetBoundaryNodeMap();
  //   bool allNodesAreOnBoundary = true;
  //   for (LO i = 0; i < Teuchos::as<LO>(boundaryNodes.size()); i++)
  //     allNodesAreOnBoundary &= boundaryNodes[i];
  //   TEST_EQUALITY(allNodesAreOnBoundary, true);
  // }

  // {
  //   Level fineLevel;
  //   TestHelpers_kokkos::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

  //   mtx->SetFixedBlockSize(3);
  //   fineLevel.Set("A", mtx);

  //   CoalesceDropFactory_kokkos dropFact = CoalesceDropFactory_kokkos();
  //   RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  //   dropFact.SetFactory("UnAmalgamationInfo",amalgFact);
  //   dropFact.SetParameter("aggregation: drop tol",Teuchos::ParameterEntry(4.1));
  //   fineLevel.Request("Graph", &dropFact);
  //   fineLevel.Request("DofsPerNode", &dropFact);

  //   dropFact.Build(fineLevel);

  //   RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos> >("Graph", &dropFact);

  //   auto boundaryNodes = graph->GetBoundaryNodeMap();
  //   bool allNodesAreOnBoundary = true;
  //   for (LO i = 0; i < Teuchos::as<LO>(boundaryNodes.size()); i++)
  //     allNodesAreOnBoundary &= boundaryNodes[i];
  //   TEST_EQUALITY(allNodesAreOnBoundary, true);
  // }

}  // AggresiveDroppingIsMarkedAsBoundary

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, 2x2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  // Set up a 2x2 matrix.
  Teuchos::RCP<Matrix> mtx = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::build2x2(lib,
                                                                                       2.0, -1.0,
                                                                                       -1.5, 2.0);
  mtx->SetFixedBlockSize(1);

  using magnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  Teuchos::RCP<Xpetra::MultiVector<magnitudeType, LocalOrdinal, GlobalOrdinal, Node>> coords;
  coords = Xpetra::MultiVectorFactory<magnitudeType, LocalOrdinal, GlobalOrdinal, Node>::Build(mtx->getRowMap(), 1);
  {
    auto lclCoords  = coords->getLocalViewHost(Tpetra::Access::OverwriteAll);
    auto rank       = comm->getRank();
    lclCoords(0, 0) = 2 * rank;
    lclCoords(1, 0) = 2 * rank + 1;
  }

  using TF                = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>;
  using local_matrix_type = typename Matrix::local_matrix_host_type;

  std::vector<Teuchos::ParameterList> params;
  std::vector<local_matrix_type> expectedFilteredMatrices;
  std::vector<std::vector<bool>> expectedBoundaryNodesVector;

  for (bool reuseGraph : {false, true}) {
    // test case 0
    Teuchos::ParameterList params0 = Teuchos::ParameterList();
    params0.set("aggregation: Dirichlet threshold", 0.);
    // dropFact.SetParameter("aggregation: row sum drop tol", Teuchos::ParameterEntry(0.2));

    params0.set("aggregation: drop scheme", "classical");
    params0.set("aggregation: use ml scaling of drop tol", false);
    params0.set("aggregation: drop tol", 0.);
    params0.set("aggregation: dropping may create Dirichlet", false);
    params0.set("filtered matrix: reuse graph", reuseGraph);
    params0.set("filtered matrix: use lumping", false);
    params.push_back(params0);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(2.0, -1.0,
                                                             -1.5, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({false, false});

    // test case 1
    Teuchos::ParameterList params1 = Teuchos::ParameterList(params0);
    params1.set("aggregation: Dirichlet threshold", 0.2);
    params.push_back(params1);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(2.0, -1.0,
                                                             -1.5, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({false, false});

    // test case 2
    Teuchos::ParameterList params2 = Teuchos::ParameterList(params0);
    params2.set("aggregation: Dirichlet threshold", 1.1);
    params.push_back(params2);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(2.0, -1.0,
                                                             -1.5, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({true, false});

    // test case 3
    Teuchos::ParameterList params3 = Teuchos::ParameterList(params0);
    params3.set("aggregation: drop tol", 0.51);
    params.push_back(params3);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(2.0, 0.0,
                                                             -1.5, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({false, false});

    // test case 4
    Teuchos::ParameterList params4 = Teuchos::ParameterList(params0);
    params4.set("aggregation: Dirichlet threshold", 1.1);
    params4.set("aggregation: drop tol", 0.51);
    params.push_back(params4);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(2.0, 0.0,
                                                             -1.5, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({true, false});

    // test case 5
    Teuchos::ParameterList params5 = Teuchos::ParameterList(params0);
    params5.set("aggregation: Dirichlet threshold", 1.1);
    params5.set("aggregation: dropping may create Dirichlet", true);
    params5.set("aggregation: drop tol", 0.51);
    params.push_back(params5);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(2.0, 0.0,
                                                             -1.5, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({true, false});

    // test case 6
    Teuchos::ParameterList params6 = Teuchos::ParameterList(params4);
    params6.set("aggregation: Dirichlet threshold", 1.1);
    params6.set("filtered matrix: use lumping", true);
    params.push_back(params6);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(1.0, 0.0,
                                                             -1.5, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({true, false});

    // test case 7
    Teuchos::ParameterList params7 = Teuchos::ParameterList(params0);
    params7.set("aggregation: drop scheme", "distance laplacian");
    params0.set("aggregation: dropping may create Dirichlet", true);
    params.push_back(params7);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(2.0, -1.0,
                                                             -1.5, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({false, false});

    // test case 8
    Teuchos::ParameterList params8 = Teuchos::ParameterList(params0);
    params8.set("aggregation: drop scheme", "distance laplacian");
    params8.set("aggregation: drop tol", 1.01);
    params.push_back(params8);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(2.0, 0.0,
                                                             0.0, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({true, true});

    // test case 9
    Teuchos::ParameterList params9 = Teuchos::ParameterList(params0);
    params9.set("aggregation: drop scheme", "classical");
    params9.set("aggregation: classical algo", "unscaled cut");
    params9.set("aggregation: drop tol", 1.0 / 3.6);
    params.push_back(params9);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(2.0, -1.0,
                                                             -1.5, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({false, false});
  }

  for (size_t testNo = 0; testNo < params.size(); ++testNo) {
    out << "\n\nRunning test number " << testNo << std::endl
        << std::endl;

    auto param                  = params[testNo];
    auto expectedFilteredMatrix = expectedFilteredMatrices[testNo];
    auto expectedBoundaryNodes  = expectedBoundaryNodesVector[testNo];

    for (bool useKokkos : {false, true}) {
      RCP<Matrix> filteredA;
      Kokkos::View<bool *, Kokkos::HostSpace> boundaryNodes;
      {
        Level fineLevel;
        TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);
        fineLevel.Set("A", mtx);
        fineLevel.Set("Coordinates", coords);

        RCP<MueLu::SingleLevelFactoryBase> dropFact;
        RCP<MueLu::SingleLevelFactoryBase> filteredAFact;
        if (useKokkos) {
          out << "\nKokkos code path\nparams:\n"
              << param << "\n";
          dropFact = rcp(new CoalesceDropFactory_kokkos());
          dropFact->SetParameterList(param);
          filteredAFact = dropFact;
        } else {
          out << "\nNon-Kokkos code path\nparams:\n"
              << param << "\n";
          dropFact            = rcp(new CoalesceDropFactory());
          filteredAFact       = rcp(new FilteredAFactory());
          auto paramCopy      = Teuchos::ParameterList(param);
          auto paramFilteredA = Teuchos::ParameterList();
          paramFilteredA.set("filtered matrix: use lumping", paramCopy.get<bool>("filtered matrix: use lumping"));
          paramCopy.remove("filtered matrix: use lumping");
          paramFilteredA.set("filtered matrix: reuse graph", paramCopy.get<bool>("filtered matrix: reuse graph"));
          paramCopy.remove("filtered matrix: reuse graph");
          paramFilteredA.set("filtered matrix: Dirichlet threshold", paramCopy.get<double>("aggregation: Dirichlet threshold"));
          dropFact->SetParameterList(paramCopy);
          filteredAFact->SetParameterList(paramFilteredA);
          filteredAFact->SetFactory("Graph", dropFact);
          filteredAFact->SetFactory("Filtering", dropFact);
        }
        fineLevel.Request("A", filteredAFact.get());
        fineLevel.Request("Graph", dropFact.get());
        fineLevel.Request("DofsPerNode", dropFact.get());
        dropFact->Build(fineLevel);

        filteredA = fineLevel.Get<RCP<Matrix>>("A", filteredAFact.get());

        if (useKokkos) {
          auto graph           = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", dropFact.get());
          auto boundaryNodes_d = graph->GetBoundaryNodeMap();
          boundaryNodes        = Kokkos::View<bool *, Kokkos::HostSpace>("boundaryNodes_host", boundaryNodes_d.extent(0));
          Kokkos::deep_copy(boundaryNodes, boundaryNodes_d);
        } else {
          auto graph    = fineLevel.Get<RCP<LWGraph>>("Graph", dropFact.get());
          boundaryNodes = graph->GetBoundaryNodeMap();
        }
      }

      auto lclFilteredA         = filteredA->getLocalMatrixHost();
      auto lclExpectedFilteredA = expectedFilteredMatrix;

      out << "Filtered A:\n"
          << TF::localMatToString(lclFilteredA);

      TEST_EQUALITY(lclFilteredA.graph.row_map.extent(0), lclExpectedFilteredA.graph.row_map.extent(0));
      for (size_t row = 0; row < std::min(lclFilteredA.graph.row_map.extent(0), lclExpectedFilteredA.graph.row_map.extent(0)); ++row)
        TEST_EQUALITY(lclFilteredA.graph.row_map(row), lclExpectedFilteredA.graph.row_map(row));

      TEST_EQUALITY(lclFilteredA.graph.entries.extent(0), lclExpectedFilteredA.graph.entries.extent(0));
      for (size_t entry = 0; entry < std::min(lclFilteredA.graph.entries.extent(0), lclExpectedFilteredA.graph.entries.extent(0)); ++entry) {
        TEST_EQUALITY(lclFilteredA.graph.entries(entry), lclExpectedFilteredA.graph.entries(entry));
        TEST_EQUALITY(lclFilteredA.values(entry), lclExpectedFilteredA.values(entry));
      }

      auto boundaryNodes_h = Kokkos::create_mirror_view(boundaryNodes);
      Kokkos::deep_copy(boundaryNodes_h, boundaryNodes);
      TEST_EQUALITY(boundaryNodes_h.extent(0), expectedBoundaryNodes.size());
      for (size_t row = 0; row < std::min(boundaryNodes_h.extent(0), expectedBoundaryNodes.size()); ++row)
        TEST_EQUALITY(boundaryNodes_h(row), expectedBoundaryNodes[row]);
    }
    out << "Done with test number " << testNo << std::endl;
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, SignedClassicalDistanceLaplacian, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = TestHelpers_kokkos::Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  GO nx = 20;
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  matrixList.set("matrixType", "Laplace1D");
  auto [A, coords, nullspace, dofsPerNode] = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildMatrixCoordsNullspace(matrixList, lib);
  // coords are uniform on [0, 1].
  if (comm->getRank() == 0) {  // move first dof to -1e3
    auto lclCoords  = coords->getLocalViewHost(Tpetra::Access::ReadWrite);
    lclCoords(0, 0) = -1e3;
  }

  // d(x_i, x_{i+-1}) ~ 1/nx except for
  // d(x_0, x_1) ~ 1e3
  //
  // -> We drop the entry (1, 0).
  //    The entry (0, 1) is kept due to the unsymmetric nature of the criterion.

  Level fineLevel;
  fineLevel.Set("A", A);
  fineLevel.Set("Coordinates", coords);
  fineLevel.Set("Nullspace", nullspace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.6));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("signed classical distance laplacian")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph_d = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  auto graph                  = graph_d->copyToHost();
  LO myDofsPerNode            = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(graph->GetGlobalNumEdges(), 57);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, SignedClassicalSADistanceLaplacian, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = TestHelpers_kokkos::Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  GO nx = 20;
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  matrixList.set("matrixType", "Laplace1D");
  auto [A, coords, nullspace, dofsPerNode] = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildMatrixCoordsNullspace(matrixList, lib);
  // coords are uniform on [0, 1].
  if (comm->getRank() == 0) {  // move first dof to -1e3
    auto lclCoords  = coords->getLocalViewHost(Tpetra::Access::ReadWrite);
    lclCoords(0, 0) = -1e3;
  }

  // d(x_i, x_{i+-1}) ~ 1/nx except for
  // d(x_0, x_1) ~ 1e3
  //
  // L_01 = L_10 = -1e-6
  // L_ij ~ -1/nx^2
  //
  // L_ii ~ 2/nx^2 for i=1, .., nx-1
  // L_00 = L_11 ~ 1/nx^2
  // L_ii ~ 1/nx^2 for i=nx
  //
  //
  // -> We drop the entries (0, 1), (1, 0) for dropTol = 0.4

  Level fineLevel;
  fineLevel.Set("A", A);
  fineLevel.Set("Coordinates", coords);
  fineLevel.Set("Nullspace", nullspace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.4));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("signed classical sa distance laplacian")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph_d = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  auto graph                  = graph_d->copyToHost();
  LO myDofsPerNode            = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(graph->GetGlobalNumEdges(), 56);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, CountNegativeDiagonals, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = TestHelpers_kokkos::Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  matrixList.set("ny", (GO)10);
  matrixList.set("nz", (GO)10);
  matrixList.set("matrixType", "Laplace3D");
  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildMatrix(matrixList, lib);

  Level fineLevel;
  fineLevel.Set("A", A);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.5));
  coalesceDropFact.SetParameter("filtered matrix: count negative diagonals", Teuchos::ParameterEntry(true));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

#define MUELU_ETI_GROUP(SC, LO, GO, NO)                                                                                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, Constructor, SC, LO, GO, NO)                                 \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, Build, SC, LO, GO, NO)                                       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, DistanceLaplacian, SC, LO, GO, NO)                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, DistanceLaplacianScaledCut, SC, LO, GO, NO)                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, DistanceLaplacianSgndRugeStuebenDistribLump, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, DistanceLaplacianUnscaledCut, SC, LO, GO, NO)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, DistanceLaplacianCutSym, SC, LO, GO, NO)                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, DistanceLaplacianScalarMaterial, SC, LO, GO, NO)             \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, DistanceLaplacianTensorMaterial, SC, LO, GO, NO)             \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, ClassicalScaledCut, SC, LO, GO, NO)                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, ClassicalUnScaledCut, SC, LO, GO, NO)                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, ClassicalCutSym, SC, LO, GO, NO)                             \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, pw_A_SA, SC, LO, GO, NO)                                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, pw_A_RS, SC, LO, GO, NO)                                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, pw_A_signedSA, SC, LO, GO, NO)                               \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, cutdrop_A_SA, SC, LO, GO, NO)                                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, cutdrop_A_RS, SC, LO, GO, NO)                                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, cutdrop_A_signedSA, SC, LO, GO, NO)                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, SignedScaledCutClassical, SC, LO, GO, NO)                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, SignedUnscaledCutClassical, SC, LO, GO, NO)                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, BlockDiagonalColoredSignedClassical, SC, LO, GO, NO)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, BlockDiagonalNoColoredSignedClassical, SC, LO, GO, NO)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, BlockDiagonalSignedClassical, SC, LO, GO, NO)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, BlockDiagonal, SC, LO, GO, NO)                               \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, BlockDiagonalVector, SC, LO, GO, NO)                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, BlockDiagonalDistanceLaplacian, SC, LO, GO, NO)              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, BlockDiagonalDistanceLaplacianWeighted, SC, LO, GO, NO)      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, DistanceLaplacianWeighted, SC, LO, GO, NO)                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, ClassicScalarWithoutFiltering, SC, LO, GO, NO)               \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, ClassicScalarWithFiltering, SC, LO, GO, NO)                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, ClassicBlockWithoutFiltering, SC, LO, GO, NO)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, AggresiveDroppingIsMarkedAsBoundary, SC, LO, GO, NO)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, 2x2, SC, LO, GO, NO)                                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, SignedClassicalDistanceLaplacian, SC, LO, GO, NO)            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, SignedClassicalSADistanceLaplacian, SC, LO, GO, NO)          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, CountNegativeDiagonals, SC, LO, GO, NO)

// TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, ClassicBlockWithFiltering,     SC, LO, GO, NO) // not implemented yet

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
