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
/// Unit tests for Tpetra::Details::GraphAssembly.
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
#include "Tpetra_Details_GraphAssembly.hpp"
#include "Tpetra_Geometry.hpp"
#include "Tpetra_FECrsGraph.hpp"
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

// ---------------------------------------------------------------------------
// A multi-block mesh: several GlobalMesh blocks that share a common global node
// numbering.  This models a mesh with more than one element type (e.g. some
// triangles and some quads) assembled into a single graph -- exactly what
// Tpetra::Geometry supports.
// ---------------------------------------------------------------------------
template <class GO>
struct MultiBlockMesh {
  size_t numGlobalNodes = 0;
  // One GlobalMesh per block; all share the same global node numbering.
  std::vector<GlobalMesh<GO>> blocks;
};

// Build a mixed mesh over an nx-by-ny grid of cells: the left half of the cells
// are split into two triangles each (block 0), and the right half are kept as
// quads (block 1).  All blocks share the standard structured node numbering, so
// nodes on the interface are shared between the two blocks.
template <class GO>
MultiBlockMesh<GO> makeMixedTriQuadMesh(int nx, int ny) {
  MultiBlockMesh<GO> mesh;
  mesh.numGlobalNodes = size_t(nx + 1) * size_t(ny + 1);
  auto nodeGID        = [nx](int i, int j) -> GO { return j * (nx + 1) + i; };

  GlobalMesh<GO> triBlock;
  triBlock.nodesPerElement = 3;
  triBlock.numGlobalNodes  = mesh.numGlobalNodes;
  GlobalMesh<GO> quadBlock;
  quadBlock.nodesPerElement = 4;
  quadBlock.numGlobalNodes  = mesh.numGlobalNodes;

  const int split = nx / 2;
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      const GO n00 = nodeGID(i, j);
      const GO n10 = nodeGID(i + 1, j);
      const GO n11 = nodeGID(i + 1, j + 1);
      const GO n01 = nodeGID(i, j + 1);
      if (i < split) {
        // Two triangles.
        triBlock.element_to_node.push_back({n00, n10, n11});
        triBlock.element_to_node.push_back({n00, n11, n01});
      } else {
        // One quad.
        quadBlock.element_to_node.push_back({n00, n10, n11, n01});
      }
    }
  }
  mesh.blocks.push_back(triBlock);
  mesh.blocks.push_back(quadBlock);
  return mesh;
}

// Ground-truth connectivity of a multi-block mesh: the union over all blocks.
template <class GO>
std::vector<std::set<GO>> expectedConnectivity(const MultiBlockMesh<GO>& mesh) {
  std::vector<std::set<GO>> adj(mesh.numGlobalNodes);
  for (const auto& block : mesh.blocks) {
    for (const auto& elem : block.element_to_node) {
      for (const GO a : elem) {
        for (const GO b : elem) {
          adj[a].insert(b);
        }
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
  using graph_assembly_type = Tpetra::Details::GraphAssembly<LO, GO, NT>;

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

// The core of the multi-block test: distribute a multi-block global mesh over
// the communicator, run the assembly (either via GraphAssembly directly with a
// Geometry, or via the free-standing Tpetra::assembleFECrsGraph overloads), and
// verify against the ground-truth connectivity.
//
// mode selects which assembly entry point is exercised:
//   0: Tpetra::Details::GraphAssembly(rowMap, ownedPlusSharedMap, geometry)
//   1: Tpetra::assembleFECrsGraph(geometry, ownedRowMap, ownedPlusSharedMap)
//   2: Tpetra::assembleFECrsGraph(geometry, ownedPlusSharedMap)  [simplified]
//   3: Tpetra::assembleFECrsGraph(geometry, ownedPlusSharedGIDs, comm) [simplified]
template <class LO, class GO, class NT>
void testMultiBlockMesh(const MultiBlockMesh<GO>& mesh, int mode,
                        Teuchos::FancyOStream& out, bool& success,
                        const std::string& label) {
  using map_type      = Tpetra::Map<LO, GO, NT>;
  using geometry_type = Tpetra::Geometry<GO, NT>;
  using e2n_type      = typename geometry_type::element_to_node_type;

  out << "=== Multi-block GraphAssembly test: " << label
      << " (mode " << mode << ") ===" << endl;
  Teuchos::OSTab tab1(out);

  auto comm              = getDefaultComm();
  const int myRank       = comm->getRank();
  const int numProcs     = comm->getSize();
  constexpr GO indexBase = 0;

  const size_t numGlobalNodes = mesh.numGlobalNodes;

  // Owned nodes: simple contiguous 1-to-1 map.
  RCP<const map_type> ownedMap(new map_type(numGlobalNodes, indexBase, comm));
  const size_t numOwnedNodes = ownedMap->getLocalNumElements();

  // Build, per block, the element-to-node connectivity of this rank's owned
  // elements (contiguous block split of each block's elements), and collect the
  // owned+shared node set.
  std::set<GO> ownedPlusSharedSet;
  for (size_t i = 0; i < numOwnedNodes; ++i)
    ownedPlusSharedSet.insert(ownedMap->getGlobalElement(i));

  geometry_type geometry;
  const int numBlocks = static_cast<int>(mesh.blocks.size());
  // Keep the device views alive for the duration of the test.
  std::vector<typename e2n_type::non_const_type> ownedE2N(numBlocks);

  for (int b = 0; b < numBlocks; ++b) {
    const auto& block                = mesh.blocks[b];
    const size_t numGlobalBlockElems = block.element_to_node.size();
    const size_t base                = numGlobalBlockElems / numProcs;
    const size_t remainder           = numGlobalBlockElems % numProcs;
    auto elemsOnRank = [&](size_t r) -> size_t { return base + (r < remainder ? 1 : 0); };
    size_t elemBegin = 0;
    for (int r = 0; r < myRank; ++r) elemBegin += elemsOnRank(r);
    const size_t numOwnedElements = elemsOnRank(myRank);

    typename e2n_type::non_const_type e2n(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "ownedE2N"),
        numOwnedElements, block.nodesPerElement);
    auto e2nHost = Kokkos::create_mirror_view(e2n);
    for (size_t e = 0; e < numOwnedElements; ++e) {
      const auto& elemNodes = block.element_to_node[elemBegin + e];
      for (int n = 0; n < block.nodesPerElement; ++n) {
        e2nHost(e, n) = elemNodes[n];
        ownedPlusSharedSet.insert(elemNodes[n]);
      }
    }
    Kokkos::deep_copy(e2n, e2nHost);
    ownedE2N[b] = e2n;
    geometry.addBlock(e2n_type(e2n));
  }

  TEST_EQUALITY_CONST(geometry.getNumBlocks(), numBlocks);

  // Build the owned+shared GID list so that it is LOCALLY FITTED to the owned
  // map: the owned GIDs come first (in the owned map's local order), followed by
  // the remaining (shared) GIDs.  This is required by GraphAssembly /
  // FECrsGraph (the owned rows must be the leading chunk of the owned+shared
  // rows), and a mixed-element mesh can have shared nodes with smaller GIDs than
  // some owned nodes, so we cannot rely on a plain sorted set here.
  std::vector<GO> ownedPlusSharedGIDs;
  ownedPlusSharedGIDs.reserve(ownedPlusSharedSet.size());
  for (size_t i = 0; i < numOwnedNodes; ++i)
    ownedPlusSharedGIDs.push_back(ownedMap->getGlobalElement(i));
  for (const GO gid : ownedPlusSharedSet)
    if (!ownedMap->isNodeGlobalElement(gid))
      ownedPlusSharedGIDs.push_back(gid);

  const Tpetra::global_size_t INVALID =
      Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  RCP<const map_type> ownedPlusSharedMap(
      new map_type(INVALID, ownedPlusSharedGIDs.data(),
                   ownedPlusSharedGIDs.size(), indexBase, comm));

  // Run the requested assembly path and obtain the owned graph.
  RCP<Tpetra::CrsGraph<LO, GO, NT>> graph;
  if (mode == 0) {
    Tpetra::Details::GraphAssembly<LO, GO, NT> assembler(ownedMap, ownedPlusSharedMap, geometry);
    assembler.build();
    graph = assembler.getGraph();
  } else if (mode == 1) {
    graph = Tpetra::assembleFECrsGraph<LO, GO, NT>(geometry, ownedMap, ownedPlusSharedMap);
  } else if (mode == 2) {
    graph = Tpetra::assembleFECrsGraph<LO, GO, NT>(geometry, ownedPlusSharedMap);
  } else {  // mode == 3
    Teuchos::ArrayView<const GO> gidView(ownedPlusSharedGIDs.data(),
                                         ownedPlusSharedGIDs.size());
    graph = Tpetra::assembleFECrsGraph<LO, GO, NT>(geometry, gidView, comm);
  }

  TEST_ASSERT(!graph.is_null());
  if (graph.is_null()) return;
  TEST_ASSERT(graph->isFillComplete());

  // For modes 0 and 1 the owned row map is the contiguous ownedMap.  For the
  // simplified modes 2 and 3 the owned map is built internally by
  // createOneToOne, so we don't assume it equals ownedMap; we instead verify
  // each owned row against the ground truth using the graph's own row map.
  const std::vector<std::set<GO>> expected = expectedConnectivity(mesh);

  using nonconst_global_inds_host_view_type =
      typename Tpetra::CrsGraph<LO, GO, NT>::nonconst_global_inds_host_view_type;

  auto rowMap = graph->getRowMap();
  for (size_t i = 0; i < rowMap->getLocalNumElements(); ++i) {
    const GO gblRow                  = rowMap->getGlobalElement(i);
    const std::set<GO>& expectedCols = expected[gblRow];

    const size_t expectedNumEntries = expectedCols.size();
    const size_t reportedNumEntries = graph->getNumEntriesInGlobalRow(gblRow);
    TEST_EQUALITY(reportedNumEntries, expectedNumEntries);

    nonconst_global_inds_host_view_type gblColInds("gblColInds", reportedNumEntries);
    size_t numColInds = 0;
    graph->getGlobalRowCopy(gblRow, gblColInds, numColInds);
    TEST_EQUALITY(numColInds, expectedNumEntries);

    std::set<GO> reportedCols;
    for (size_t k = 0; k < numColInds; ++k) reportedCols.insert(gblColInds(k));
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

// Multi-block (mixed triangle + quad) mesh, exercised through each assembly
// entry point (GraphAssembly with a Geometry, and the various free-standing
// assembleFECrsGraph overloads including the simplified Feature-2 ones).
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, GraphAssembly_MultiBlock, LO, GO, NT) {
  // 4 x 4 grid: left 2 columns of cells are triangles (2 per cell = 16 tris),
  // right 2 columns are quads (8 quads).  25 nodes total.
  const auto mesh = makeMixedTriQuadMesh<GO>(4, 4);
  for (int mode = 0; mode <= 3; ++mode) {
    testMultiBlockMesh<LO, GO, NT>(mesh, mode, out, success, "mixed tri+quad (2D)");
  }
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP(LO, GO, NT)                                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_Quad, LO, GO, \
                                       NT)                                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_Tri, LO, GO,  \
                                       NT)                                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_Tet, LO, GO,  \
                                       NT)                                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_MultiBlock,   \
                                       LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LGN(UNIT_TEST_GROUP)
