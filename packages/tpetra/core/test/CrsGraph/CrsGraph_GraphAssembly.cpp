// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file CrsGraph_GraphAssembly.cpp
///
/// Unit tests for Tpetra::Experimental::GraphAssembly.
///
/// The finite-element example directory
/// (core/example/Finite-Element-Assembly) only exercises GraphAssembly on a 2D
/// structured quad mesh.  GraphAssembly itself is written generically over the
/// number of nodes per element, so these tests check that it produces the
/// correct sparse graph for three different element geometries:
///
///   1. quads       (2D, 4 nodes per element)
///   2. triangles   (2D, 3 nodes per element)
///   3. tetrahedra  (3D, 4 nodes per element)
///
/// Each test builds a small global mesh, assigns
/// nodes and elements to MPI ranks, and then the inputs required by GraphAssembly:
/// - the owned map
/// - ownedElementToNode, a 2D view describing element-to-node connectivity of the owned elements
/// - the owned+shared map
/// It then uses GraphAssembly to generate the graph, and checks it has the expected connectivity
/// based on the mesh structure.
///
/// The tests are written to run correctly on either 1 or 4 MPI ranks (like most
/// Tpetra tests).

#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_GraphAssembly.hpp"
#include "Teuchos_CommHelpers.hpp"  // REDUCE_MIN, reduceAll
#include "Tpetra_TestingUtilities.hpp"

#include <algorithm>
#include <set>
#include <vector>

using std::endl;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using Tpetra::TestingUtilities::getDefaultComm;

// ---------------------------------------------------------------------------
// A small, self-contained description of a global finite-element mesh.  The
// entire global mesh is listed explicitly (this is a correctness test on tiny
// meshes, so that is fine).  Elements are described by their global node IDs.
// ---------------------------------------------------------------------------
template <class GO>
struct GlobalMesh {
  size_t numGlobalNodes = 0;
  int nodesPerElement   = 0;
  // element_to_node[e] is the list of global node IDs of global element e.
  std::vector<std::vector<GO>> element_to_node;

  size_t numGlobalElements() const {
    return element_to_node.size();
  }
};

// Build a 2D structured grid of quad elements, nx by ny elements.
// Node (i,j), 0 <= i <= nx, 0 <= j <= ny, has GID j*(nx+1) + i.
template <class GO>
GlobalMesh<GO> makeQuadMesh(int nx, int ny) {
  GlobalMesh<GO> mesh;
  mesh.nodesPerElement = 4;
  mesh.numGlobalNodes  = size_t(nx + 1) * size_t(ny + 1);
  auto nodeGID         = [nx](int i, int j) -> GO {
    return j * (nx + 1) + i;
  };
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      // Counter-clockwise ordering of the quad's four corners.
      mesh.element_to_node.push_back({nodeGID(i, j), nodeGID(i + 1, j),
                                      nodeGID(i + 1, j + 1), nodeGID(i, j + 1)});
    }
  }
  return mesh;
}

// Build a 2D structured grid of triangles: nx by ny cells, each split into two
// triangles.  Node numbering matches the quad grid.
template <class GO>
GlobalMesh<GO> makeTriMesh(int nx, int ny) {
  GlobalMesh<GO> mesh;
  mesh.nodesPerElement = 3;
  mesh.numGlobalNodes  = size_t(nx + 1) * size_t(ny + 1);
  auto nodeGID         = [nx](int i, int j) -> GO {
    return j * (nx + 1) + i;
  };
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      const GO n00 = nodeGID(i, j);
      const GO n10 = nodeGID(i + 1, j);
      const GO n11 = nodeGID(i + 1, j + 1);
      const GO n01 = nodeGID(i, j + 1);
      // Lower-left triangle and upper-right triangle of the cell.
      mesh.element_to_node.push_back({n00, n10, n11});
      mesh.element_to_node.push_back({n00, n11, n01});
    }
  }
  return mesh;
}

// Build a 3D structured grid of tetrahedra: nx by ny by nz hex cells, each hex
// split into 6 tets using the standard Kuhn / diagonal (0,7) decomposition.
// Node (i,j,k) has GID k*(nx+1)*(ny+1) + j*(nx+1) + i.
template <class GO>
GlobalMesh<GO> makeTetMesh(int nx, int ny, int nz) {
  GlobalMesh<GO> mesh;
  mesh.nodesPerElement = 4;
  mesh.numGlobalNodes  = size_t(nx + 1) * size_t(ny + 1) * size_t(nz + 1);
  auto nodeGID         = [nx, ny](int i, int j, int k) -> GO {
    return k * (nx + 1) * (ny + 1) + j * (nx + 1) + i;
  };
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        // The 8 corners of the hex cell.  Local corner index c has bit 0 = di,
        // bit 1 = dj, bit 2 = dk.
        GO c[8];
        for (int b = 0; b < 8; ++b) {
          const int di = b & 1;
          const int dj = (b >> 1) & 1;
          const int dk = (b >> 2) & 1;
          c[b]         = nodeGID(i + di, j + dj, k + dk);
        }
        // Kuhn decomposition of a cube into 6 tets sharing the main diagonal
        // 0 - 7.  Each tet is {0, 7, a, b} where a->b walks the cube edges.
        static const int tets[6][4] = {
            {0, 1, 3, 7}, {0, 3, 2, 7}, {0, 2, 6, 7}, {0, 6, 4, 7}, {0, 4, 5, 7}, {0, 5, 1, 7}};
        for (int t = 0; t < 6; ++t) {
          mesh.element_to_node.push_back(
              {c[tets[t][0]], c[tets[t][1]], c[tets[t][2]], c[tets[t][3]]});
        }
      }
    }
  }
  return mesh;
}

// Given the global mesh, compute the expected connectivity of every global
// node: node N's row contains exactly the union of the nodes of every element
// that N belongs to (including N itself).  This is the independent "ground
// truth" the assembled graph is checked against.
template <class GO>
std::vector<std::set<GO>> expectedConnectivity(const GlobalMesh<GO>& mesh) {
  std::vector<std::set<GO>> adj(mesh.numGlobalNodes);
  for (const auto& elem : mesh.element_to_node) {
    for (const GO a : elem) {
      for (const GO b : elem) {
        adj[a].insert(b);
      }
    }
  }
  return adj;
}

// The core of the test: distribute the given global mesh over the
// communicator, run GraphAssembly, and verify.
template <class LO, class GO, class NT>
void testMesh(const GlobalMesh<GO>& mesh, Teuchos::FancyOStream& out,
              bool& success, const std::string& label) {
  using map_type            = Tpetra::Map<LO, GO, NT>;
  using graph_assembly_type = Tpetra::Experimental::GraphAssembly<LO, GO, NT>;

  out << "=== GraphAssembly test: " << label << " ===" << endl;
  Teuchos::OSTab tab1(out);

  auto comm              = getDefaultComm();
  const int myRank       = comm->getRank();
  const int numProcs     = comm->getSize();
  constexpr GO indexBase = 0;

  const size_t numGlobalNodes    = mesh.numGlobalNodes;
  const size_t numGlobalElements = mesh.numGlobalElements();

  // ---- Decide which nodes this rank owns (simple contiguous 1-to-1 map). ----
  RCP<const map_type> ownedMap(
      new map_type(numGlobalNodes, indexBase, comm));
  const size_t numOwnedNodes = ownedMap->getLocalNumElements();

  // ---- Decide which elements this rank owns (contiguous block split). ----
  // Rank r owns the global elements in [elemBegin, elemEnd).
  const size_t base      = numGlobalElements / numProcs;
  const size_t remainder = numGlobalElements % numProcs;
  auto elemsOnRank       = [&](size_t r) -> size_t {
    return base + (r < remainder ? 1 : 0);
  };
  size_t elemBegin = 0;
  for (int r = 0; r < myRank; ++r) elemBegin += elemsOnRank(r);
  const size_t elemEnd          = elemBegin + elemsOnRank(myRank);
  const size_t numOwnedElements = elemEnd - elemBegin;

  // ---- Build the element-to-node connectivity of the owned elements. ----
  typename graph_assembly_type::element_to_node_type::non_const_type
      ownedElementToNode(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "ownedElementToNode"),
          numOwnedElements, mesh.nodesPerElement);
  auto ownedElementToNodeHost =
      Kokkos::create_mirror_view(ownedElementToNode);
  for (size_t e = 0; e < numOwnedElements; ++e) {
    const size_t globalElem = elemBegin + e;
    const auto& elemNodes   = mesh.element_to_node[globalElem];
    for (int n = 0; n < mesh.nodesPerElement; ++n) {
      ownedElementToNodeHost(e, n) = elemNodes[n];
    }
  }
  Kokkos::deep_copy(ownedElementToNode, ownedElementToNodeHost);

  // ---- Build the owned+shared map: every node adjacent to an owned element. --
  // Loop over all nodes of this rank's owned elements, collect the unique set.
  std::set<GO> ownedPlusSharedSet;
  // Seed with the owned nodes so they always come first / are always present.
  for (size_t i = 0; i < numOwnedNodes; ++i) {
    ownedPlusSharedSet.insert(ownedMap->getGlobalElement(i));
  }
  for (size_t e = 0; e < numOwnedElements; ++e) {
    const size_t globalElem = elemBegin + e;
    const auto& elemNodes   = mesh.element_to_node[globalElem];
    for (const GO node : elemNodes) ownedPlusSharedSet.insert(node);
  }
  std::vector<GO> ownedPlusSharedGIDs(ownedPlusSharedSet.begin(),
                                      ownedPlusSharedSet.end());
  const Tpetra::global_size_t INVALID =
      Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  RCP<const map_type> ownedPlusSharedMap(
      new map_type(INVALID, ownedPlusSharedGIDs.data(),
                   ownedPlusSharedGIDs.size(), indexBase,
                   comm));

  // ---- Run the assembly. ----
  typename graph_assembly_type::element_to_node_type ownedElementToNodeConst =
      ownedElementToNode;
  graph_assembly_type assembler(ownedMap, ownedPlusSharedMap,
                                ownedElementToNodeConst);
  assembler.build();
  auto graph = assembler.getGraph();
  TEST_ASSERT(!graph.is_null());
  if (graph.is_null()) return;
  TEST_ASSERT(graph->isFillComplete());

  // The assembled graph's row map must be the owned map.
  TEST_ASSERT(graph->getRowMap()->isSameAs(*ownedMap));

  // ---- Verify each owned row against the independently-computed ground truth.
  const std::vector<std::set<GO>> expected = expectedConnectivity(mesh);

  using nonconst_global_inds_host_view_type =
      typename Tpetra::CrsGraph<LO, GO, NT>::nonconst_global_inds_host_view_type;

  for (size_t i = 0; i < numOwnedNodes; ++i) {
    const GO gblRow                  = ownedMap->getGlobalElement(i);
    const std::set<GO>& expectedCols = expected[gblRow];

    const size_t expectedNumEntries = expectedCols.size();
    const size_t reportedNumEntries =
        graph->getNumEntriesInGlobalRow(gblRow);
    TEST_EQUALITY(reportedNumEntries, expectedNumEntries);

    nonconst_global_inds_host_view_type gblColInds("gblColInds",
                                                   reportedNumEntries);
    size_t numColInds = 0;
    graph->getGlobalRowCopy(gblRow, gblColInds, numColInds);
    TEST_EQUALITY(numColInds, expectedNumEntries);

    // Collect the reported columns into a set and compare against expected.
    std::set<GO> reportedCols;
    for (size_t k = 0; k < numColInds; ++k) {
      reportedCols.insert(gblColInds(k));
    }
    TEST_EQUALITY(reportedCols.size(), expectedNumEntries);
    const bool colsMatch = (reportedCols == expectedCols);
    TEST_ASSERT(colsMatch);
    if (!colsMatch) {
      out << "Row " << gblRow << " mismatch.\n  expected: {";
      for (const GO c : expectedCols) out << c << " ";
      out << "}\n  got:      {";
      for (const GO c : reportedCols) out << c << " ";
      out << "}" << endl;
    }
  }

  // Make the pass/fail collective across all ranks.
  int lclSuccess = success ? 1 : 0;
  int gblSuccess = 0;
  reduceAll<int, int>(*comm, REDUCE_MIN, lclSuccess, outArg(gblSuccess));
  TEST_EQUALITY_CONST(gblSuccess, 1);
}

//
// UNIT TESTS
//

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, GraphAssembly_Quad, LO, GO, NT) {
  // 4 x 4 grid of quads: 16 elements, 25 nodes.  16 / 4 = 4 elements per rank.
  const auto mesh = makeQuadMesh<GO>(4, 4);
  testMesh<LO, GO, NT>(mesh, out, success, "quads (2D)");
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, GraphAssembly_Tri, LO, GO, NT) {
  // 2 x 4 grid of cells split into triangles: 16 triangles, 15 nodes.
  // 16 / 4 = 4 elements per rank.
  const auto mesh = makeTriMesh<GO>(2, 4);
  testMesh<LO, GO, NT>(mesh, out, success, "triangles (2D)");
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, GraphAssembly_Tet, LO, GO, NT) {
  // 2 x 2 x 1 grid of hexes, each split into 6 tets: 24 tets, 18 nodes.
  // 24 / 4 = 6 elements per rank.
  const auto mesh = makeTetMesh<GO>(2, 2, 1);
  testMesh<LO, GO, NT>(mesh, out, success, "tetrahedra (3D)");
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP(LO, GO, NT)                                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_Quad, LO, GO, \
                                       NT)                                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_Tri, LO, GO,  \
                                       NT)                                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_Tet, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LGN(UNIT_TEST_GROUP)
