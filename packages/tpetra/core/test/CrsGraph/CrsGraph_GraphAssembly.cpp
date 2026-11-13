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
/// Unit tests for Tpetra::Details::GraphAssembly, Tpetra::Geometry and the
/// Tpetra::assembleFECrsGraph / Tpetra::assembleFECrsGraphWithMaps functions.
///
/// There are two kinds of tests:
///
///   - Ground truth: build a small global mesh (quads, triangles, tetrahedra,
///     or a mix of element types in several blocks), distribute it, assemble
///     the graph, and check every owned row against the connectivity computed
///     independently from the global mesh.
///
///   - Identical to insertion: assemble the graph, and also build the
///     FECrsGraph the traditional way (insertGlobalIndices between
///     beginAssembly() and endAssembly()) with the same maps.  The two must be
///     identical: same row, column, domain and range maps (the column map
///     compared index by index), and the same local row offsets and column
///     indices, in the same order, for both the owned and the owned+shared
///     graphs.  These tests cover node ownerships where an element may own
///     none of its nodes, rectangular blocks (different row and column
///     spaces), padded blocks (elements coupled to a neighbor's nodes), every
///     combination of optional maps, and Geometries with local indices.
///
/// The tests are written to run correctly on 1 to 4 MPI ranks.

#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Details_GraphAssembly.hpp"
#include "Tpetra_Geometry.hpp"
#include "Tpetra_FECrsGraph.hpp"
#include "Tpetra_Assembly_Helpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_CommHelpers.hpp"  // REDUCE_MIN, reduceAll
#include "Tpetra_TestingUtilities.hpp"

#include <algorithm>
#include <functional>
#include <map>
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

// Split each block of a multi-block mesh into numSplits smaller blocks of the
// same element type (assigning elements round-robin).  This produces a mesh
// with many blocks, some of which may have no elements on a given rank.
template <class GO>
MultiBlockMesh<GO> splitBlocks(const MultiBlockMesh<GO>& mesh, int numSplits) {
  MultiBlockMesh<GO> split;
  split.numGlobalNodes = mesh.numGlobalNodes;
  for (const auto& block : mesh.blocks) {
    std::vector<GlobalMesh<GO>> pieces(numSplits);
    for (auto& piece : pieces) {
      piece.nodesPerElement = block.nodesPerElement;
      piece.numGlobalNodes  = block.numGlobalNodes;
    }
    for (size_t e = 0; e < block.element_to_node.size(); ++e)
      pieces[e % numSplits].element_to_node.push_back(block.element_to_node[e]);
    for (auto& piece : pieces) split.blocks.push_back(piece);
  }
  return split;
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

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Copy a list of rows (each a list of indices) into a rank-2 device view, padding
// short rows with Tpetra's invalid index (which Geometry ignores).
template <class T, class DeviceType>
Kokkos::View<T**, DeviceType> toPaddedDeviceView(const std::vector<std::vector<T>>& rows,
                                                 const std::string& label) {
  size_t width = 0;
  for (const auto& r : rows) width = std::max(width, r.size());
  Kokkos::View<T**, DeviceType> view(label, rows.size(), width);
  auto host = Kokkos::create_mirror_view(view);
  for (size_t i = 0; i < rows.size(); ++i)
    for (size_t j = 0; j < width; ++j)
      host(i, j) = j < rows[i].size() ? rows[i][j] : Teuchos::OrdinalTraits<T>::invalid();
  Kokkos::deep_copy(view, host);
  return view;
}

// Make a map with the given global indices, in the given order, on this process.
template <class LO, class GO, class NT>
RCP<const Tpetra::Map<LO, GO, NT>> makeMap(const std::vector<GO>& gids) {
  const Tpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  return rcp(new Tpetra::Map<LO, GO, NT>(INVALID, gids.data(), gids.size(), GO(0), getDefaultComm()));
}

// Make the pass/fail state of a test collective across all processes.
inline void collectiveSuccess(bool& success, Teuchos::FancyOStream& out) {
  int lclSuccess = success ? 1 : 0;
  int gblSuccess = 0;
  reduceAll<int, int>(*getDefaultComm(), REDUCE_MIN, lclSuccess, outArg(gblSuccess));
  TEST_EQUALITY_CONST(gblSuccess, 1);
}

// ---------------------------------------------------------------------------
// Comparing graphs
// ---------------------------------------------------------------------------

// Whether two maps have the same global indices in the same order on this process.
template <class MapType>
bool sameLocalGIDs(const MapType& a, const MapType& b) {
  if (a.getLocalNumElements() != b.getLocalNumElements()) return false;
  auto aGIDs = a.getMyGlobalIndices();
  auto bGIDs = b.getMyGlobalIndices();
  for (size_t i = 0; i < aGIDs.extent(0); ++i)
    if (aGIDs(i) != bGIDs(i)) return false;
  return true;
}

// Whether two fill-complete graphs are identical on this process: same maps
// (index by index) and same local row offsets and column indices, in order.
template <class GraphType>
bool identicalLocalGraphs(const GraphType& ref, const GraphType& test, const std::string& label) {
  const int myRank = getDefaultComm()->getRank();
  auto report      = [&](const std::string& what) {
    std::cerr << "Proc " << myRank << ": " << label << ": " << what << std::endl;
    return false;
  };
  if (!ref.isFillComplete() || !test.isFillComplete()) return report("not fill complete");
  if (!sameLocalGIDs(*ref.getRowMap(), *test.getRowMap())) return report("row maps differ");
  if (!sameLocalGIDs(*ref.getColMap(), *test.getColMap())) return report("column maps differ");
  if (!sameLocalGIDs(*ref.getDomainMap(), *test.getDomainMap())) return report("domain maps differ");
  if (!sameLocalGIDs(*ref.getRangeMap(), *test.getRangeMap())) return report("range maps differ");
  auto refGraph  = ref.getLocalGraphHost();
  auto testGraph = test.getLocalGraphHost();
  if (refGraph.row_map.extent(0) != testGraph.row_map.extent(0)) return report("number of rows differs");
  for (size_t i = 0; i < refGraph.row_map.extent(0); ++i)
    if (refGraph.row_map(i) != testGraph.row_map(i))
      return report("row offsets differ at local row " + std::to_string(i));
  const size_t nnz = refGraph.row_map.extent(0) ? refGraph.row_map(refGraph.row_map.extent(0) - 1) : 0;
  for (size_t k = 0; k < nnz; ++k)
    if (refGraph.entries(k) != testGraph.entries(k))
      return report("column indices differ at local entry " + std::to_string(k));
  return true;
}

// Whether two fully-assembled FECrsGraphs have identical owned and owned+shared
// graphs on every process.  Collective.
template <class FEGraphType>
bool identicalFEGraphs(FEGraphType& ref, FEGraphType& test) {
  bool ok = identicalLocalGraphs(ref, test, "owned graph");
  ref.switchActiveCrsGraph();
  test.switchActiveCrsGraph();
  ok = identicalLocalGraphs(ref, test, "owned+shared graph") && ok;
  ref.switchActiveCrsGraph();
  test.switchActiveCrsGraph();
  int lclOK = ok ? 1 : 0;
  int gblOK = 0;
  reduceAll<int, int>(*getDefaultComm(), REDUCE_MIN, lclOK, outArg(gblOK));
  return gblOK == 1;
}

// ---------------------------------------------------------------------------
// Ground-truth tests
// ---------------------------------------------------------------------------

// Check every row of graph (in its row map) against the expected connectivity.
template <class GraphType, class GO>
void checkAgainstGroundTruth(const GraphType& graph, const std::vector<std::set<GO>>& expected,
                             Teuchos::FancyOStream& out, bool& success) {
  using global_inds_type = typename GraphType::nonconst_global_inds_host_view_type;
  auto rowMap            = graph.getRowMap();
  for (size_t i = 0; i < rowMap->getLocalNumElements(); ++i) {
    const GO gblRow                  = rowMap->getGlobalElement(i);
    const std::set<GO>& expectedCols = expected[gblRow];
    const size_t numEntries          = graph.getNumEntriesInGlobalRow(gblRow);
    TEST_EQUALITY(numEntries, expectedCols.size());
    global_inds_type gblColInds("gblColInds", numEntries);
    size_t numColInds = 0;
    graph.getGlobalRowCopy(gblRow, gblColInds, numColInds);
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
}

// The owned+shared GIDs of a process: the owned GIDs (in owned map order) then
// the other GIDs of nodeSet in ascending order -- locally fitted to the owned
// map.  If fitted is false, all GIDs in ascending order instead (generally not
// locally fitted).
template <class LO, class GO, class NT>
std::vector<GO> ownedPlusSharedGIDs(const Tpetra::Map<LO, GO, NT>& ownedMap, const std::set<GO>& nodeSet,
                                    const bool fitted) {
  std::vector<GO> gids;
  if (fitted) {
    for (size_t i = 0; i < ownedMap.getLocalNumElements(); ++i) gids.push_back(ownedMap.getGlobalElement(i));
    for (const GO gid : nodeSet)
      if (!ownedMap.isNodeGlobalElement(gid)) gids.push_back(gid);
  } else {
    std::set<GO> all(nodeSet);
    for (size_t i = 0; i < ownedMap.getLocalNumElements(); ++i) all.insert(ownedMap.getGlobalElement(i));
    gids.assign(all.begin(), all.end());
  }
  return gids;
}

// Distribute a single-block mesh, assemble it with GraphAssembly, and check it.
template <class LO, class GO, class NT>
void testMesh(const GlobalMesh<GO>& mesh, Teuchos::FancyOStream& out, bool& success, const std::string& label) {
  using map_type      = Tpetra::Map<LO, GO, NT>;
  using geometry_type = Tpetra::Geometry<LO, GO, NT>;

  out << "=== GraphAssembly test: " << label << " ===" << endl;
  Teuchos::OSTab tab1(out);

  auto comm          = getDefaultComm();
  const int myRank   = comm->getRank();
  const int numProcs = comm->getSize();

  // Owned nodes: a contiguous 1-to-1 map.  Owned elements: a contiguous split.
  RCP<const map_type> ownedMap(new map_type(mesh.numGlobalNodes, GO(0), comm));
  const size_t numGlobalElements = mesh.numGlobalElements();
  const size_t base              = numGlobalElements / numProcs;
  const size_t remainder         = numGlobalElements % numProcs;
  auto elemsOnRank               = [&](size_t r) -> size_t { return base + (r < remainder ? 1 : 0); };
  size_t elemBegin               = 0;
  for (int r = 0; r < myRank; ++r) elemBegin += elemsOnRank(r);
  std::vector<std::vector<GO>> ownedElements(mesh.element_to_node.begin() + elemBegin,
                                             mesh.element_to_node.begin() + elemBegin + elemsOnRank(myRank));
  std::set<GO> nodeSet;
  for (const auto& elem : ownedElements) nodeSet.insert(elem.begin(), elem.end());
  auto ownedPlusSharedMap = makeMap<LO, GO, NT>(ownedPlusSharedGIDs(*ownedMap, nodeSet, true));

  geometry_type geometry;
  geometry.addGlobalBlock(toPaddedDeviceView<GO, typename NT::device_type>(ownedElements, "e2n"));
  Tpetra::Details::GraphAssembly<LO, GO, NT> assembly(geometry, ownedMap, ownedPlusSharedMap);
  assembly.build();
  auto graph = assembly.getGraph();
  TEST_ASSERT(!graph.is_null());
  if (!graph.is_null()) {
    TEST_ASSERT(graph->isFillComplete());
    TEST_ASSERT(graph->getRowMap()->isSameAs(*ownedMap));
    checkAgainstGroundTruth(*graph, expectedConnectivity(mesh), out, success);
  }
  collectiveSuccess(success, out);
}

// Distribute a multi-block mesh, assemble it, and check it against the ground
// truth.  mode selects the entry point:
//   0: Tpetra::Details::GraphAssembly
//   1: Tpetra::assembleFECrsGraphWithMaps(geometry, ownedMap, ownedPlusSharedMap)
//   2: Tpetra::assembleFECrsGraph(geometry, ownedPlusSharedMap)
//   3: Tpetra::assembleFECrsGraph(geometry, ownedPlusSharedGIDs, comm)
// For modes 2 and 3, the owned+shared GIDs are in ascending order, which is
// generally not locally fitted to the owned map, so these modes also exercise
// the fitting done by the simplified interface.  If localIndices is true, the
// Geometry holds local indices of the owned+shared map passed in.
template <class LO, class GO, class NT>
void testMultiBlockMesh(const MultiBlockMesh<GO>& mesh, const int mode, const bool localIndices,
                        Teuchos::FancyOStream& out, bool& success, const std::string& label) {
  using map_type      = Tpetra::Map<LO, GO, NT>;
  using geometry_type = Tpetra::Geometry<LO, GO, NT>;
  using device_type   = typename NT::device_type;

  out << "=== Multi-block test: " << label << " (mode " << mode << ", "
      << (localIndices ? "local" : "global") << " indices) ===" << endl;
  Teuchos::OSTab tab1(out);

  auto comm          = getDefaultComm();
  const int myRank   = comm->getRank();
  const int numProcs = comm->getSize();

  RCP<const map_type> ownedMap(new map_type(mesh.numGlobalNodes, GO(0), comm));

  // Each block's elements are split contiguously over the processes.
  std::vector<std::vector<std::vector<GO>>> ownedBlocks;
  std::set<GO> nodeSet;
  for (const auto& block : mesh.blocks) {
    const size_t n         = block.element_to_node.size();
    const size_t base      = n / numProcs;
    const size_t remainder = n % numProcs;
    auto elemsOnRank       = [&](size_t r) -> size_t { return base + (r < remainder ? 1 : 0); };
    size_t elemBegin       = 0;
    for (int r = 0; r < myRank; ++r) elemBegin += elemsOnRank(r);
    ownedBlocks.emplace_back(block.element_to_node.begin() + elemBegin,
                             block.element_to_node.begin() + elemBegin + elemsOnRank(myRank));
    for (const auto& elem : ownedBlocks.back()) nodeSet.insert(elem.begin(), elem.end());
  }
  const std::vector<GO> opsGIDs = ownedPlusSharedGIDs(*ownedMap, nodeSet, mode <= 1);
  auto ownedPlusSharedMap       = makeMap<LO, GO, NT>(opsGIDs);

  geometry_type geometry;
  for (const auto& block : ownedBlocks) {
    if (localIndices) {
      std::vector<std::vector<LO>> lids;
      for (const auto& elem : block) {
        lids.emplace_back();
        for (const GO gid : elem) lids.back().push_back(ownedPlusSharedMap->getLocalElement(gid));
      }
      geometry.addLocalBlock(toPaddedDeviceView<LO, device_type>(lids, "e2n"));
    } else {
      geometry.addGlobalBlock(toPaddedDeviceView<GO, device_type>(block, "e2n"));
    }
  }
  TEST_EQUALITY(geometry.getNumBlocks(), static_cast<int>(mesh.blocks.size()));
  TEST_EQUALITY(geometry.hasLocalIndices(), localIndices);

  RCP<const Tpetra::CrsGraph<LO, GO, NT>> graph;
  if (mode == 0) {
    Tpetra::Details::GraphAssembly<LO, GO, NT> assembly(geometry, ownedMap, ownedPlusSharedMap);
    assembly.build();
    graph = assembly.getGraph();
  } else if (mode == 1) {
    graph = Tpetra::assembleFECrsGraphWithMaps(geometry, ownedMap, ownedPlusSharedMap);
  } else if (mode == 2) {
    graph = Tpetra::assembleFECrsGraph(geometry, ownedPlusSharedMap);
  } else {
    graph = Tpetra::assembleFECrsGraph(geometry, Teuchos::ArrayView<const GO>(opsGIDs), comm);
  }
  TEST_ASSERT(!graph.is_null());
  if (!graph.is_null()) {
    TEST_ASSERT(graph->isFillComplete());
    // The simplified modes build their own owned map (with createOneToOne), so
    // check each row of the graph's own row map.
    checkAgainstGroundTruth(*graph, expectedConnectivity(mesh), out, success);
  }
  collectiveSuccess(success, out);
}

// ---------------------------------------------------------------------------
// Identical-to-insertion tests
// ---------------------------------------------------------------------------

// A global finite-element problem: blocks of elements, each element with its
// row indices and (for rectangular problems) column indices; and the owning
// process of each element, row index, and column index.
template <class GO>
struct FEProblem {
  struct Element {
    std::vector<GO> rows;
    std::vector<GO> cols;  // unused unless rectangular
  };
  std::vector<std::vector<Element>> blocks;
  bool rectangular = false;
  std::function<int(size_t globalElement, int numProcs)> elementOwner;
  std::function<int(GO gid, int numProcs)> rowOwner;
  std::function<int(GO gid, int numProcs)> colOwner;  // rectangular only

  const std::vector<GO>& colsOf(const Element& e) const { return rectangular ? e.cols : e.rows; }
};

// Options for which optional maps are given to the assembly.
struct MapOptions {
  bool giveColMap       = false;  // the owned+shared column map
  bool giveDomainMap    = false;  // the owned+shared domain map
  bool localIndices     = false;  // Geometry with local indices
  std::string describe() const {
    return std::string(giveColMap ? "colMap" : "no colMap") + ", " +
           (giveDomainMap ? "ownedPlusSharedDomainMap" : "no ownedPlusSharedDomainMap") + ", " +
           (localIndices ? "local" : "global") + " indices";
  }
};

// Distribute the problem, assemble its graph with assembleFECrsGraphWithMaps,
// build the same graph by insertion, and check that the two are identical.
template <class LO, class GO, class NT>
void testIdenticalToInsertion(const FEProblem<GO>& problem, const MapOptions& opts,
                              Teuchos::FancyOStream& out, bool& success, const std::string& label) {
  using map_type      = Tpetra::Map<LO, GO, NT>;
  using import_type   = Tpetra::Import<LO, GO, NT>;
  using fe_graph_type = Tpetra::FECrsGraph<LO, GO, NT>;
  using geometry_type = Tpetra::Geometry<LO, GO, NT>;
  using device_type   = typename NT::device_type;
  using Element       = typename FEProblem<GO>::Element;

  out << "=== Identical to insertion: " << label << " (" << opts.describe() << ") ===" << endl;
  Teuchos::OSTab tab1(out);

  auto comm          = getDefaultComm();
  const int myRank   = comm->getRank();
  const int numProcs = comm->getSize();
  auto rowOwner      = problem.rowOwner;
  auto colOwner      = problem.rectangular ? problem.colOwner : problem.rowOwner;

  // This process' elements, and the global row and column index sets.
  std::vector<std::vector<Element>> myBlocks;
  std::set<GO> allRows, allCols, myRows, myCols, mergedCols;
  size_t globalElement = 0, numGlobalElements = 0, maxCols = 1;
  for (const auto& block : problem.blocks) {
    myBlocks.emplace_back();
    for (const auto& elem : block) {
      const auto& cols = problem.colsOf(elem);
      allRows.insert(elem.rows.begin(), elem.rows.end());
      allCols.insert(cols.begin(), cols.end());
      maxCols = std::max(maxCols, cols.size());
      if (problem.elementOwner(globalElement++, numProcs) == myRank) {
        myBlocks.back().push_back(elem);
        myRows.insert(elem.rows.begin(), elem.rows.end());
        myCols.insert(cols.begin(), cols.end());
      }
      // The columns of the owned rows, after merging all processes' contributions.
      for (const GO row : elem.rows)
        if (rowOwner(row, numProcs) == myRank) mergedCols.insert(cols.begin(), cols.end());
      ++numGlobalElements;
    }
  }

  // The maps.  Owned+shared maps: owned indices first, then the others, sorted.
  std::vector<GO> ownedRowGIDs, ownedColGIDs;
  for (const GO gid : allRows)
    if (rowOwner(gid, numProcs) == myRank) ownedRowGIDs.push_back(gid);
  for (const GO gid : allCols)
    if (colOwner(gid, numProcs) == myRank) ownedColGIDs.push_back(gid);
  auto withOthers = [](std::vector<GO> gids, const std::set<GO>& others) {
    const std::set<GO> first(gids.begin(), gids.end());
    for (const GO gid : others)
      if (!first.count(gid)) gids.push_back(gid);
    return gids;
  };
  const std::vector<GO> opsColGIDs = withOthers(ownedColGIDs, myCols);
  auto ownedRowMap                 = makeMap<LO, GO, NT>(ownedRowGIDs);
  auto ownedPlusSharedRowMap       = makeMap<LO, GO, NT>(withOthers(ownedRowGIDs, myRows));
  auto ownedColMap                 = makeMap<LO, GO, NT>(ownedColGIDs);
  auto ownedPlusSharedColGIDsMap   = makeMap<LO, GO, NT>(opsColGIDs);
  // A column map must also hold the columns that other processes contribute to
  // the owned rows.
  RCP<const map_type> colMap =
      opts.giveColMap ? makeMap<LO, GO, NT>(withOthers(opsColGIDs, mergedCols)) : Teuchos::null;
  RCP<const map_type> ownedPlusSharedDomainMap = opts.giveDomainMap ? ownedPlusSharedColGIDsMap : Teuchos::null;
  // The owned domain map: the owned columns (the default, ownedRowMap, if square).
  RCP<const map_type> ownedDomainMap = problem.rectangular ? ownedColMap : Teuchos::null;

  // The reference: an FECrsGraph built by insertion, with the same maps.
  RCP<fe_graph_type> ref;
  const size_t maxEntries = numGlobalElements * maxCols;
  const RCP<const import_type> noImporter;
  if (colMap.is_null()) {
    if (ownedPlusSharedDomainMap.is_null())
      ref = rcp(new fe_graph_type(ownedRowMap, ownedPlusSharedRowMap, maxEntries, noImporter, ownedDomainMap));
    else
      ref = rcp(new fe_graph_type(ownedRowMap, ownedPlusSharedRowMap, maxEntries, ownedPlusSharedDomainMap,
                                  noImporter, ownedDomainMap));
  } else {
    if (ownedPlusSharedDomainMap.is_null())
      ref = rcp(new fe_graph_type(ownedRowMap, ownedPlusSharedRowMap, colMap, maxEntries, noImporter, ownedDomainMap));
    else
      ref = rcp(new fe_graph_type(ownedRowMap, ownedPlusSharedRowMap, colMap, maxEntries, ownedPlusSharedDomainMap,
                                  noImporter, ownedDomainMap));
  }
  {
    // Elements may own none of their nodes.
    auto params = Teuchos::parameterList();
    params->set("Check Col GIDs In At Least One Owned Row", false);
    ref->setParameterList(params);
  }
  Tpetra::beginAssembly(*ref);
  for (const auto& block : myBlocks)
    for (const auto& elem : block) {
      const auto& cols = problem.colsOf(elem);
      if (colMap.is_null()) {
        for (const GO row : elem.rows)
          ref->insertGlobalIndices(row, Teuchos::ArrayView<const GO>(cols));
      } else {
        // With a column map, the graph is locally indexed.
        std::vector<LO> lclCols;
        for (const GO col : cols) lclCols.push_back(colMap->getLocalElement(col));
        for (const GO row : elem.rows)
          ref->insertLocalIndices(ownedPlusSharedRowMap->getLocalElement(row), Teuchos::ArrayView<const LO>(lclCols));
      }
    }
  Tpetra::endAssembly(*ref);

  // The Geometry.
  RCP<const map_type> asmColMap = !colMap.is_null()                     ? colMap
                                  : !ownedPlusSharedDomainMap.is_null() ? ownedPlusSharedDomainMap
                                                                        : ownedPlusSharedRowMap;
  geometry_type geometry;
  for (const auto& block : myBlocks) {
    std::vector<std::vector<GO>> rows, cols;
    for (const auto& elem : block) {
      rows.push_back(elem.rows);
      cols.push_back(elem.cols);
    }
    if (opts.localIndices) {
      auto toLocal = [](const std::vector<std::vector<GO>>& gids, const map_type& map) {
        std::vector<std::vector<LO>> lids;
        for (const auto& r : gids) {
          lids.emplace_back();
          for (const GO gid : r) lids.back().push_back(map.getLocalElement(gid));
        }
        return lids;
      };
      auto localRows = toPaddedDeviceView<LO, device_type>(toLocal(rows, *ownedPlusSharedRowMap), "rows");
      if (problem.rectangular) {
        // Pass the columns as a strided subview of a wider array, to exercise
        // LayoutStride.
        auto padded = toPaddedDeviceView<LO, device_type>(toLocal(cols, *asmColMap), "cols");
        Kokkos::View<LO**, Kokkos::LayoutRight, device_type> wide("wide", padded.extent(0), padded.extent(1) + 1);
        auto strided = Kokkos::subview(wide, Kokkos::ALL(), Kokkos::make_pair(size_t(0), size_t(padded.extent(1))));
        Kokkos::deep_copy(strided, padded);
        geometry.addLocalBlock(localRows, strided);
      } else {
        geometry.addLocalBlock(localRows);
      }
    } else {
      auto globalRows = toPaddedDeviceView<GO, device_type>(rows, "rows");
      if (problem.rectangular)
        geometry.addGlobalBlock(globalRows, toPaddedDeviceView<GO, device_type>(cols, "cols"));
      else
        geometry.addGlobalBlock(globalRows);
    }
  }

  auto assembled = Tpetra::assembleFECrsGraphWithMaps(geometry, ownedRowMap, ownedPlusSharedRowMap, colMap,
                                                      ownedPlusSharedDomainMap, Teuchos::null, ownedDomainMap);
  TEST_ASSERT(!assembled.is_null());
  if (!assembled.is_null()) {
    const bool identical = identicalFEGraphs(*ref, *assembled);
    TEST_ASSERT(identical);
  }
  collectiveSuccess(success, out);
}

// Run testIdenticalToInsertion with every valid combination of options.
template <class LO, class GO, class NT>
void testIdenticalAllOptions(const FEProblem<GO>& problem, Teuchos::FancyOStream& out, bool& success,
                             const std::string& label) {
  for (const bool giveColMap : {false, true})
    for (const bool giveDomainMap : {false, true})
      for (const bool localIndices : {false, true}) {
        // A rectangular problem's columns are not in the owned+shared row map.
        if (problem.rectangular && !giveColMap && !giveDomainMap) continue;
        MapOptions opts;
        opts.giveColMap    = giveColMap;
        opts.giveDomainMap = giveDomainMap;
        opts.localIndices  = localIndices;
        testIdenticalToInsertion<LO, GO, NT>(problem, opts, out, success, label);
      }
}

// Square problems from a multi-block mesh.  Nodes and elements are assigned to
// processes round-robin, so that elements often own none of their nodes, and
// processes own nodes that none of their elements touch.
template <class GO>
FEProblem<GO> squareProblem(const MultiBlockMesh<GO>& mesh, const bool roundRobin) {
  FEProblem<GO> problem;
  size_t numElements = 0;
  for (const auto& block : mesh.blocks) {
    problem.blocks.emplace_back();
    for (const auto& elem : block.element_to_node) problem.blocks.back().push_back({elem, {}});
    numElements += block.element_to_node.size();
  }
  const GO numNodes = mesh.numGlobalNodes;
  if (roundRobin) {
    problem.elementOwner = [](size_t e, int p) { return int(e % p); };
    problem.rowOwner     = [](GO gid, int p) { return int((gid * 7 + 3) % p); };
  } else {
    problem.elementOwner = [numElements](size_t e, int p) { return int(e * p / numElements); };
    problem.rowOwner     = [numNodes](GO gid, int p) { return int(gid * p / numNodes); };
  }
  return problem;
}

template <class GO>
MultiBlockMesh<GO> singleBlock(const GlobalMesh<GO>& mesh) {
  MultiBlockMesh<GO> multi;
  multi.numGlobalNodes = mesh.numGlobalNodes;
  multi.blocks.push_back(mesh);
  return multi;
}

// A rectangular "two-field" problem on a quad mesh: the rows are the nodes; the
// columns are the nodes and one cell-centered unknown per element (column
// numNodes + element), owned by the element's process.  Each process' owned
// rows receive cell-centered columns of other processes' elements.
template <class GO>
FEProblem<GO> twoFieldProblem(int nx, int ny) {
  const auto mesh = makeQuadMesh<GO>(nx, ny);
  const GO numNodes = mesh.numGlobalNodes;
  FEProblem<GO> problem;
  problem.rectangular = true;
  problem.blocks.emplace_back();
  for (size_t e = 0; e < mesh.element_to_node.size(); ++e) {
    auto cols = mesh.element_to_node[e];
    cols.push_back(numNodes + GO(e));
    problem.blocks.back().push_back({mesh.element_to_node[e], cols});
  }
  problem.elementOwner = [](size_t e, int p) { return int(e % p); };
  problem.rowOwner     = [](GO gid, int p) { return int((gid * 7 + 3) % p); };
  problem.colOwner     = [numNodes](GO gid, int p) {
    return gid < numNodes ? int((gid * 7 + 3) % p) : int((gid - numNodes) % p);
  };
  return problem;
}

// A neighbor-coupled problem on a quad mesh: each element's rows are its nodes;
// its columns are its nodes and the nodes of the element to its right, if any
// (so the column lists have different lengths, and the column block is padded).
// This models e.g. Panzer's "associated neighbors".
template <class GO>
FEProblem<GO> neighborProblem(int nx, int ny) {
  const auto mesh = makeQuadMesh<GO>(nx, ny);
  FEProblem<GO> problem;
  problem.rectangular = true;
  problem.blocks.emplace_back();
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const auto& nodes = mesh.element_to_node[j * nx + i];
      auto cols         = nodes;
      if (i + 1 < nx) {
        const auto& right = mesh.element_to_node[j * nx + i + 1];
        cols.insert(cols.end(), right.begin(), right.end());
      }
      problem.blocks.back().push_back({nodes, cols});
    }
  problem.elementOwner = [](size_t e, int p) { return int((e / 3) % p); };
  problem.rowOwner     = [](GO gid, int p) { return int((gid * 5 + 1) % p); };
  problem.colOwner     = problem.rowOwner;
  return problem;
}

// A square block whose elements have different numbers of nodes (triangles and
// quads in one block, padded to 4 nodes).
template <class GO>
MultiBlockMesh<GO> paddedTriQuadMesh(int nx, int ny) {
  const auto mixed = makeMixedTriQuadMesh<GO>(nx, ny);
  MultiBlockMesh<GO> mesh;
  mesh.numGlobalNodes = mixed.numGlobalNodes;
  GlobalMesh<GO> block;
  block.numGlobalNodes  = mixed.numGlobalNodes;
  block.nodesPerElement = 4;
  // Interleave the triangles and quads.
  const auto& tris  = mixed.blocks[0].element_to_node;
  const auto& quads = mixed.blocks[1].element_to_node;
  for (size_t k = 0; k < std::max(tris.size(), quads.size()); ++k) {
    if (k < tris.size()) block.element_to_node.push_back(tris[k]);
    if (k < quads.size()) block.element_to_node.push_back(quads[k]);
  }
  mesh.blocks.push_back(block);
  return mesh;
}

//
// UNIT TESTS
//

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, GraphAssembly_Quad, LO, GO, NT) {
  testMesh<LO, GO, NT>(makeQuadMesh<GO>(4, 4), out, success, "quads (2D)");
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, GraphAssembly_Tri, LO, GO, NT) {
  testMesh<LO, GO, NT>(makeTriMesh<GO>(2, 4), out, success, "triangles (2D)");
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, GraphAssembly_Tet, LO, GO, NT) {
  testMesh<LO, GO, NT>(makeTetMesh<GO>(2, 2, 1), out, success, "tetrahedra (3D)");
}

// Mixed triangle + quad mesh in two blocks, through every entry point, with
// global and local indices.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, GraphAssembly_MultiBlock, LO, GO, NT) {
  const auto mesh = makeMixedTriQuadMesh<GO>(4, 4);
  for (int mode = 0; mode <= 3; ++mode)
    for (const bool localIndices : {false, true})
      testMultiBlockMesh<LO, GO, NT>(mesh, mode, localIndices, out, success, "mixed tri+quad (2D)");
}

// Many blocks (some empty on some processes).
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, GraphAssembly_ManyBlocks, LO, GO, NT) {
  const auto mesh = splitBlocks(makeMixedTriQuadMesh<GO>(4, 4), 7);
  for (int mode = 0; mode <= 3; ++mode)
    for (const bool localIndices : {false, true})
      testMultiBlockMesh<LO, GO, NT>(mesh, mode, localIndices, out, success, "many-block tri+quad (2D)");
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, GraphAssembly_IdenticalQuad, LO, GO, NT) {
  for (const bool roundRobin : {false, true})
    testIdenticalAllOptions<LO, GO, NT>(squareProblem(singleBlock(makeQuadMesh<GO>(5, 4)), roundRobin), out, success,
                                        roundRobin ? "quads, round-robin" : "quads, contiguous");
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, GraphAssembly_IdenticalTet, LO, GO, NT) {
  for (const bool roundRobin : {false, true})
    testIdenticalAllOptions<LO, GO, NT>(squareProblem(singleBlock(makeTetMesh<GO>(2, 2, 2)), roundRobin), out, success,
                                        roundRobin ? "tets, round-robin" : "tets, contiguous");
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, GraphAssembly_IdenticalManyBlocks, LO, GO, NT) {
  testIdenticalAllOptions<LO, GO, NT>(squareProblem(splitBlocks(makeMixedTriQuadMesh<GO>(4, 4), 3), true), out,
                                      success, "many-block tri+quad, round-robin");
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, GraphAssembly_IdenticalPadded, LO, GO, NT) {
  testIdenticalAllOptions<LO, GO, NT>(squareProblem(paddedTriQuadMesh<GO>(4, 4), true), out, success,
                                      "tri+quad padded in one block, round-robin");
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, GraphAssembly_IdenticalRectangular, LO, GO, NT) {
  testIdenticalAllOptions<LO, GO, NT>(twoFieldProblem<GO>(5, 3), out, success, "two fields (rectangular)");
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, GraphAssembly_IdenticalNeighbors, LO, GO, NT) {
  testIdenticalAllOptions<LO, GO, NT>(neighborProblem<GO>(5, 3), out, success, "neighbor-coupled (padded columns)");
}

// The comparison used by the tests above must detect a different graph.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, GraphAssembly_ComparisonDetectsDifferences, LO, GO, NT) {
  using map_type      = Tpetra::Map<LO, GO, NT>;
  using geometry_type = Tpetra::Geometry<LO, GO, NT>;
  using device_type   = typename NT::device_type;
  auto comm           = getDefaultComm();
  // Assemble a quad mesh and a triangle mesh on the same nodes, with the same maps.
  const auto quads = makeQuadMesh<GO>(2, 4);
  const auto tris  = makeTriMesh<GO>(2, 4);
  RCP<const map_type> ownedMap(new map_type(quads.numGlobalNodes, GO(0), comm));
  auto assemble = [&](const GlobalMesh<GO>& mesh) {
    // Every process has every element, so the owned+shared map has every node.
    std::set<GO> nodes;
    for (GO n = 0; n < GO(mesh.numGlobalNodes); ++n) nodes.insert(n);
    auto ops = makeMap<LO, GO, NT>(ownedPlusSharedGIDs(*ownedMap, nodes, true));
    geometry_type geometry;
    geometry.addGlobalBlock(toPaddedDeviceView<GO, device_type>(mesh.element_to_node, "e2n"));
    return Tpetra::assembleFECrsGraphWithMaps(geometry, ownedMap, ops);
  };
  auto quadGraph = assemble(quads);
  auto triGraph  = assemble(tris);
  auto quadAgain = assemble(quads);
  TEST_ASSERT(identicalFEGraphs(*quadGraph, *quadAgain));
  TEST_ASSERT(!identicalFEGraphs(*quadGraph, *triGraph));
}

// Errors: mixing kinds of index, mismatched rectangular views, indices not in
// the maps, and an owned+shared row map that is not locally fitted.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, GraphAssembly_Errors, LO, GO, NT) {
  using map_type      = Tpetra::Map<LO, GO, NT>;
  using geometry_type = Tpetra::Geometry<LO, GO, NT>;
  using device_type   = typename NT::device_type;
  auto comm           = getDefaultComm();

  Kokkos::View<GO**, device_type> globalBlock("g", 2, 3);
  Kokkos::View<GO**, device_type> globalBlock3("g3", 3, 3);
  Kokkos::View<LO**, device_type> localBlock("l", 2, 3);
  {
    geometry_type geometry;
    geometry.addGlobalBlock(globalBlock);
    TEST_THROW(geometry.addLocalBlock(localBlock), std::invalid_argument);
    TEST_THROW(geometry.addGlobalBlock(globalBlock, globalBlock3), std::invalid_argument);
  }

  // Every process has one triangle of nodes {0, 1, 2}; the owned+shared map is
  // {0, 1, 2} (owned first).
  const GO numNodes = 3;
  RCP<const map_type> ownedMap(new map_type(numNodes, GO(0), comm));
  std::set<GO> nodes{0, 1, 2};
  auto ops = makeMap<LO, GO, NT>(ownedPlusSharedGIDs(*ownedMap, nodes, true));
  {
    // A node that is not in the owned+shared map.
    geometry_type geometry;
    geometry.addGlobalBlock(toPaddedDeviceView<GO, device_type>({{0, 1, numNodes + 5}}, "e2n"));
    TEST_THROW(Tpetra::assembleFECrsGraphWithMaps(geometry, ownedMap, ops), std::invalid_argument);
  }
  {
    // A local index out of range.
    geometry_type geometry;
    geometry.addLocalBlock(toPaddedDeviceView<LO, device_type>({{0, 1, LO(numNodes + 5)}}, "e2n"));
    TEST_THROW(Tpetra::assembleFECrsGraphWithMaps(geometry, ownedMap, ops), std::invalid_argument);
  }
  if (comm->getSize() > 1) {
    // An owned+shared row map that is not locally fitted on some process: the
    // shared indices first.
    std::vector<GO> reversed = ownedPlusSharedGIDs(*ownedMap, nodes, true);
    std::reverse(reversed.begin(), reversed.end());
    auto notFitted = makeMap<LO, GO, NT>(reversed);
    geometry_type geometry;
    geometry.addGlobalBlock(toPaddedDeviceView<GO, device_type>({{0, 1, 2}}, "e2n"));
    TEST_THROW(Tpetra::assembleFECrsGraphWithMaps(geometry, ownedMap, notFitted), std::invalid_argument);
  }
  collectiveSuccess(success, out);
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP(LO, GO, NT)                                                                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_Quad, LO, GO, NT)                       \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_Tri, LO, GO, NT)                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_Tet, LO, GO, NT)                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_MultiBlock, LO, GO, NT)                 \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_ManyBlocks, LO, GO, NT)                 \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_IdenticalQuad, LO, GO, NT)              \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_IdenticalTet, LO, GO, NT)               \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_IdenticalManyBlocks, LO, GO, NT)        \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_IdenticalPadded, LO, GO, NT)            \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_IdenticalRectangular, LO, GO, NT)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_IdenticalNeighbors, LO, GO, NT)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_ComparisonDetectsDifferences, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, GraphAssembly_Errors, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LGN(UNIT_TEST_GROUP)
