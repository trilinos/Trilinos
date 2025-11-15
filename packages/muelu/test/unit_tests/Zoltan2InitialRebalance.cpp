// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MueLu_config.hpp"

#ifdef HAVE_MUELU_ZOLTAN2

#include <vector>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_IO.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <MueLu_RepartitionHeuristicFactory.hpp>
#include <MueLu_RepartitionFactory.hpp>
#include <MueLu_ZoltanInterface.hpp>

#include <MueLu_SingleLevelFactoryBase.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include "MueLu_RepartitionFactory_decl.hpp"

#include "MueLu_Zoltan2Interface.hpp"

#include <vector>
#include <algorithm>
#include <utility>
#include <cmath>

namespace MueLuTests {

Teuchos::ParameterList
repartition_parameters() {
  Teuchos::ParameterList pl;
  pl.set("algorithm", "multijagged");
  pl.set("partitioning_approach", "partition");
  pl.set("mj_premigration_coordinate_count", 1000);
  pl.set("mj_premigration_option", 1);
  return pl;
}

struct Point {
  double x, y;
};

double cross_product(const Point& p1, const Point& p2, const Point& p3) {
  return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
}

std::vector<Point> compute_convex_hull(std::vector<Point>& points) {
  std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
  });

  std::vector<Point> lowerHull;
  for (const auto& p : points) {
    while (lowerHull.size() >= 2 && cross_product(lowerHull[lowerHull.size() - 2], lowerHull.back(), p) <= 0) {
      lowerHull.pop_back();
    }
    lowerHull.push_back(p);
  }

  std::vector<Point> upperHull;
  for (auto it = points.rbegin(); it != points.rend(); ++it) {
    while (upperHull.size() >= 2 && cross_product(upperHull[upperHull.size() - 2], upperHull.back(), *it) <= 0) {
      upperHull.pop_back();
    }
    upperHull.push_back(*it);
  }

  lowerHull.pop_back();
  upperHull.pop_back();

  lowerHull.insert(lowerHull.end(), upperHull.begin(), upperHull.end());
  return lowerHull;
}

bool segments_intersect(const Point& p1, const Point& p2, const Point& p3, const Point& p4) {
  double d1 = cross_product(p3, p4, p1);
  double d2 = cross_product(p3, p4, p2);
  double d3 = cross_product(p1, p2, p3);
  double d4 = cross_product(p1, p2, p4);

  // points on opposite sides of the segment
  if ((d1 * d2 < 0) && (d3 * d4 < 0)) {
    return true;
  }

  // colinear
  if (d1 == 0 && std::min(p3.x, p4.x) <= p1.x && p1.x <= std::max(p3.x, p4.x) &&
      std::min(p3.y, p4.y) <= p1.y && p1.y <= std::max(p3.y, p4.y)) {
    return true;
  }
  if (d2 == 0 && std::min(p3.x, p4.x) <= p2.x && p2.x <= std::max(p3.x, p4.x) &&
      std::min(p3.y, p4.y) <= p2.y && p2.y <= std::max(p3.y, p4.y)) {
    return true;
  }
  if (d3 == 0 && std::min(p1.x, p2.x) <= p3.x && p3.x <= std::max(p1.x, p2.x) &&
      std::min(p1.y, p2.y) <= p3.y && p3.y <= std::max(p1.y, p2.y)) {
    return true;
  }
  if (d4 == 0 && std::min(p1.x, p2.x) <= p4.x && p4.x <= std::max(p1.x, p2.x) &&
      std::min(p1.y, p2.y) <= p4.y && p4.y <= std::max(p1.y, p2.y)) {
    return true;
  }

  return false;
}

bool convex_hulls_intersect(const std::vector<Point>& hull1, const std::vector<Point>& hull2) {
  for (size_t i = 0; i < hull1.size(); ++i) {
    for (size_t j = 0; j < hull2.size(); ++j) {
      Point p1 = hull1[i];
      Point p2 = hull1[(i + 1) % hull1.size()];
      Point p3 = hull2[j];
      Point p4 = hull2[(j + 1) % hull2.size()];

      if (segments_intersect(p1, p2, p3, p4)) {
        return true;
      }
    }
  }
  return false;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Zoltan2Repartition, DeterminePartition, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  out << "Tests the algorithm for assigning partitions to PIDs." << std::endl;
  out << std::endl;

  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();
  int numProcs                       = comm->getSize();
  int myRank                         = comm->getRank();

  if (numProcs != 4) {
    std::cout << "\nThis test must be run on 4 processors!\n"
              << std::endl;
    return;
  }

  const GlobalOrdinal nx = 4, ny = 6;

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  matrixList.set("ny", ny);
  matrixList.set("keepBCs", false);

  // Describes the initial layout of matrix rows across processors.
  const GlobalOrdinal numGlobalElements = nx * ny;  // 24
  const size_t numMyElements            = myRank == 0 ? numGlobalElements : 0;
  const GlobalOrdinal indexBase         = 0;

  RCP<const Map> rank0Map = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, numMyElements, indexBase, comm);

  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector>> Pr =
      Galeri::Xpetra::BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, CrsMatrixWrap, MultiVector>("Laplace2D", rank0Map, matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();
  auto coords   = Pr->BuildCoords();

  auto params        = repartition_parameters();
  auto decomposition = MueLu::Zoltan2Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ComputeDecomposition(numProcs, A, coords, params);
  auto [GIDs, _]     = MueLu::RepartitionUtilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ConstructGIDs(decomposition);

  auto rowMap = A->getRowMap();
  auto newMap = MapFactory ::Build(TestHelpers::Parameters::getLib(), rowMap->getGlobalNumElements(), GIDs(), indexBase, comm);

  RCP<Import> importer = ImportFactory::Build(A->getRowMap(), newMap);
  auto targetMap       = importer->getTargetMap();

  Teuchos::ParameterList XpetraList;
  auto newMatrix  = MatrixFactory::Build(A, *importer, *importer, targetMap, targetMap, rcp(&XpetraList, false));
  auto distCoords = MultiVectorFactory::Build(newMap, coords->getNumVectors());
  distCoords->doImport(*coords, *importer, Xpetra::INSERT);

  TEST_EQUALITY(A->getGlobalNumRows(), numGlobalElements);
  TEST_EQUALITY(A->getLocalNumRows(), numMyElements);
  TEST_EQUALITY(coords->getGlobalLength(), numGlobalElements);
  TEST_EQUALITY(coords->getLocalLength(), numMyElements);

  const auto expectedNumLocalRows = numGlobalElements / numProcs;
  TEST_EQUALITY(newMatrix->getLocalNumRows(), expectedNumLocalRows);
  TEST_EQUALITY(distCoords->getLocalLength(), expectedNumLocalRows);
  TEST_EQUALITY(distCoords->getGlobalLength(), numGlobalElements);

  // gather coordinates on all procs so we can check the quality of the repartition
  std::vector<double> globalCoords_x(numGlobalElements);
  std::vector<double> globalCoords_y(numGlobalElements);
  {
    std::vector<double> localCoords_x(expectedNumLocalRows);
    std::vector<double> localCoords_y(expectedNumLocalRows);

    auto distCoordsView = distCoords->getLocalViewHost(Tpetra::Access::ReadOnly);
    for (auto i = 0; i < distCoordsView.extent(0); i++) {
      localCoords_x[i] = distCoordsView(i, 0);
      localCoords_y[i] = distCoordsView(i, 1);
    }

    Teuchos::gatherAll<int, double>(*comm, expectedNumLocalRows, localCoords_x.data(), numGlobalElements, globalCoords_x.data());
    Teuchos::gatherAll<int, double>(*comm, expectedNumLocalRows, localCoords_y.data(), numGlobalElements, globalCoords_y.data());
  }

  int numHullsIntersectLocal = 0;
  if (myRank == 0) {
    std::vector<std::vector<Point>> processorLocalCoords;
    processorLocalCoords.resize(numProcs);
    for (int rank = 0; rank < numProcs; ++rank) {
      processorLocalCoords[rank].resize(expectedNumLocalRows);
      for (auto i = 0; i < expectedNumLocalRows; i++) {
        const auto gRow               = expectedNumLocalRows * rank + i;
        processorLocalCoords[rank][i] = Point{
            .x = globalCoords_x[gRow],
            .y = globalCoords_y[gRow]};
      }
    }

    for (int rank_i = 0; rank_i < numProcs; ++rank_i) {
      auto hull_i = compute_convex_hull(processorLocalCoords[rank_i]);
      for (int rank_j = rank_i + 1; rank_j < numProcs; ++rank_j) {
        auto hull_j = compute_convex_hull(processorLocalCoords[rank_j]);

        const auto hullsIntersect = convex_hulls_intersect(hull_i, hull_j);
        if (hullsIntersect) numHullsIntersectLocal++;
      }
    }
  }

  int numHullsIntersectGlobal = -1;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &numHullsIntersectLocal, &numHullsIntersectGlobal);

  // If the convex hulls of the partition intersect each other, the partition is off.
  TEST_EQUALITY(numHullsIntersectGlobal, 0);
}

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Zoltan2Repartition, DeterminePartition, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests

#endif
