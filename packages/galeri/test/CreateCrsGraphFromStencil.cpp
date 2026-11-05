// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Galeri_XpetraMaps.hpp"
#include "Galeri_XpetraProblemFactory.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CreateCrsGraphFromStencil.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using scalar_type = double;
using map_type    = Tpetra::Map<>;
using matrix_type = Tpetra::CrsMatrix<>;
using vector_type = Tpetra::MultiVector<>;
using graph_type  = Tpetra::CrsGraph<>;
using LO          = typename map_type::local_ordinal_type;
using GO          = typename map_type::global_ordinal_type;
using offset_type = Tpetra::Experimental::StencilOffset;

struct TestCase {
  std::string problem;
  Teuchos::Array<GO> globalDimensions;
  Teuchos::Array<offset_type> offsets;
  LO dofsPerNode;
  int ranks;
};

Teuchos::Array<offset_type> cross(const int dimensions) {
  Teuchos::Array<offset_type> offsets{{0, 0, 0}, {-1, 0, 0}, {1, 0, 0}};
  if (dimensions > 1) {
    offsets.push_back({0, -1, 0});
    offsets.push_back({0, 1, 0});
  }
  if (dimensions > 2) {
    offsets.push_back({0, 0, -1});
    offsets.push_back({0, 0, 1});
  }
  return offsets;
}

Teuchos::Array<offset_type> box(const int dimensions) {
  Teuchos::Array<offset_type> offsets;
  for (int z = dimensions > 2 ? -1 : 0; z <= (dimensions > 2 ? 1 : 0); ++z)
    for (int y = dimensions > 1 ? -1 : 0; y <= (dimensions > 1 ? 1 : 0); ++y)
      for (int x = -1; x <= 1; ++x)
        offsets.push_back({x, y, z});
  return offsets;
}

TestCase getTestCase(const std::string& name) {
  if (name == "Serial1DThreePoint")
    return {"Laplace1D", {8}, cross(1), 1, 1};
  if (name == "Serial2DFivePoint")
    return {"Laplace2D", {5, 4}, cross(2), 1, 1};
  if (name == "Serial3DSevenPoint")
    return {"Laplace3D", {4, 3, 5}, cross(3), 1, 1};
  if (name == "Serial2DNinePointTwoDofs")
    return {"Elasticity2D", {5, 4}, box(2), 2, 1};
  if (name == "Mpi2DFivePoint")
    return {"Laplace2D", {8, 6}, cross(2), 1, 4};
  if (name == "Mpi3D27PointThreeDofs")
    return {"Elasticity3D", {6, 4, 4}, box(3), 3, 4};
  if (name == "Mpi2DBigStarRadiusTwo")
    return {"BigStar2D", {14, 12}, {{0, 0, 0}, {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {-1, -1, 0}, {1, -1, 0}, {-1, 1, 0}, {1, 1, 0}, {-2, 0, 0}, {2, 0, 0}, {0, -2, 0}, {0, 2, 0}}, 1, 4};
  throw std::invalid_argument("Unknown test case: " + name);
}

Teuchos::RCP<const map_type> expandMap(
    const Teuchos::RCP<const map_type>& nodeMap, const LO dofsPerNode) {
  Teuchos::Array<GO> gids(nodeMap->getLocalNumElements() * dofsPerNode);
  size_t output       = 0;
  const auto nodeGids = nodeMap->getMyGlobalIndices();
  for (size_t i = 0; i < nodeGids.extent(0); ++i) {
    const GO node = nodeGids(i);
    for (LO dof = 0; dof < dofsPerNode; ++dof)
      gids[output++] = node * dofsPerNode + dof;
  }
  return Teuchos::rcp(new map_type(
      nodeMap->getGlobalNumElements() * dofsPerNode, gids(), 0,
      nodeMap->getComm()));
}

std::vector<GO> globalColumns(const graph_type& graph, const LO row) {
  graph_type::local_inds_host_view_type columns;
  graph.getLocalRowView(row, columns);
  std::vector<GO> result(columns.size());
  for (size_t i = 0; i < result.size(); ++i)
    result[i] = graph.getColMap()->getGlobalElement(columns[i]);
  std::sort(result.begin(), result.end());
  return result;
}

void compareGraphs(const graph_type& actual, const graph_type& expected) {
  if (!actual.getRowMap()->isSameAs(*expected.getRowMap()) ||
      actual.getGlobalNumEntries() != expected.getGlobalNumEntries() ||
      actual.getLocalNumEntries() != expected.getLocalNumEntries())
    throw std::runtime_error("Graph maps or entry counts differ");
  for (size_t row = 0; row < actual.getLocalNumRows(); ++row)
    if (globalColumns(actual, static_cast<LO>(row)) !=
        globalColumns(expected, static_cast<LO>(row)))
      throw std::runtime_error("Columns differ in local row " +
                               std::to_string(row));
}

void runTest(const TestCase& test,
             const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  if (comm->getSize() != test.ranks)
    throw std::runtime_error("Wrong number of MPI ranks");

  Teuchos::ParameterList parameters;
  const int dimensions = test.globalDimensions.size();
  parameters.set("nx", test.globalDimensions[0]);
  if (dimensions > 1) parameters.set("ny", test.globalDimensions[1]);
  if (dimensions > 2) parameters.set("nz", test.globalDimensions[2]);
  const std::string mapType = dimensions == 1 ? "Cartesian1D" : dimensions == 2 ? "Cartesian2D"
                                                                                : "Cartesian3D";
  auto nodeMap              = Galeri::Xpetra::CreateMap<LO, GO, map_type>(
      mapType, comm, parameters);
  Teuchos::RCP<const map_type> rowMap = nodeMap;
  if (test.dofsPerNode != 1)
    rowMap = expandMap(nodeMap, test.dofsPerNode);
  auto problem = Galeri::Xpetra::BuildProblem<
      scalar_type, LO, GO, map_type, matrix_type, vector_type>(
      test.problem, rowMap, parameters);
  auto reference = problem->BuildMatrix()->getCrsGraph();

  Teuchos::Array<LO> localDimensions(dimensions);
  localDimensions[0] = parameters.get<LO>("lnx");
  if (dimensions > 1) localDimensions[1] = parameters.get<LO>("lny");
  if (dimensions > 2) localDimensions[2] = parameters.get<LO>("lnz");
  const Teuchos::Array<LO>& constLocalDimensions = localDimensions;
  auto actual                                    = Tpetra::Experimental::createCrsGraphFromStencil(
                                         rowMap, constLocalDimensions(), test.offsets(), test.dofsPerNode, true);
  compareGraphs(*actual, *reference);
}

}  // namespace

int main(int argc, char* argv[]) {
  Tpetra::ScopeGuard scope(&argc, &argv);
  const auto comm = Tpetra::getDefaultComm();
  int success     = argc == 2;
  try {
    if (success) runTest(getTestCase(argv[1]), comm);
  } catch (const std::exception& error) {
    std::cerr << "Rank " << comm->getRank() << ": " << error.what() << '\n';
    success = 0;
  }
  int globalSuccess = 0;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, success,
                     Teuchos::outArg(globalSuccess));
  if (comm->getRank() == 0)
    std::cout << "End Result: TEST "
              << (globalSuccess ? "PASSED" : "FAILED") << '\n';
  return globalSuccess ? EXIT_SUCCESS : EXIT_FAILURE;
}
