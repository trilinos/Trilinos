// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_UnitTestHarness.hpp>
#include <MueLu_UnitTestHelpers.hpp>  // provides TEUCHOS_UNIT_TEST_TEMPLATE_5_INSTANT, etc.
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "MueLu_TestHelpers_HO.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MatrixMatrix.hpp>

#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_CreateXpetraPreconditioner.hpp"

#ifdef HAVE_MUELU_INTREPID2
#include "MueLu_IntrepidPCoarsenFactory.hpp"
#include "MueLu_IntrepidPCoarsenFactory_def.hpp"  // Why does ETI suddenly decide to hate right here?
#include "Intrepid2_Types.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"
//#include "Intrepid2_HGRAD_TRI_Cn_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_Cn_FEM.hpp"
//#include "Intrepid2_HGRAD_TET_Cn_FEM.hpp"
#include "Kokkos_DynRankView.hpp"

namespace MueLuTests {

/**** some helper methods and classes by Nate ****/
#ifndef TEST_MORE_COMBINATIONS
static const int MAX_LINE_DEGREE = Intrepid2::Parameters::MaxOrder;
static const int MAX_QUAD_DEGREE = Intrepid2::Parameters::MaxOrder;
static const int MAX_HEX_DEGREE  = (Intrepid2::Parameters::MaxOrder < 4) ? Intrepid2::Parameters::MaxOrder : 4;
static const int MAX_RANK_COUNT  = 4;
#else
static const int MAX_LINE_DEGREE = Intrepid2::Parameters::MaxOrder;
static const int MAX_QUAD_DEGREE = Intrepid2::Parameters::MaxOrder;
static const int MAX_HEX_DEGREE  = Intrepid2::Parameters::MaxOrder;
static const int MAX_RANK_COUNT  = 16;
#endif

using namespace std;
// pair is subcell dim, subcell ordinal in cellTopo.  Includes (spaceDim-1, sideOrdinal), where spaceDim is the extent of the cellTopo
std::vector<std::pair<int, int>> subcellEntitiesForSide(const shards::CellTopology &cellTopo, int sideOrdinal) {
  using namespace std;
  vector<pair<int, int>> subcellEntities;
  set<int> nodesForSide;
  int spaceDim = (int)cellTopo.getDimension();
  if (spaceDim == 0) return {{}};
  int sideDim   = spaceDim - 1;
  int nodeCount = (int)cellTopo.getNodeCount(sideDim, sideOrdinal);
  // first, collect all the nodes that match the side
  for (int nodeOrdinal = 0; nodeOrdinal < nodeCount; nodeOrdinal++) {
    int node = (int)cellTopo.getNodeMap(sideDim, sideOrdinal, nodeOrdinal);
    nodesForSide.insert(node);
  }
  // now, iterate over extents.
  // Any subcells that only have nodes that match nodesForSide should be included.
  for (int d = 0; d <= sideDim; d++) {
    int subcellCount = cellTopo.getSubcellCount(d);
    for (int subcord = 0; subcord < subcellCount; subcord++) {
      bool allNodesMatch = true;
      if (d == 0) {
        // subcord is a node; just check whether that node is in nodesForSide
        allNodesMatch = (nodesForSide.find(subcord) != nodesForSide.end());
      } else {
        int scNodeCount = cellTopo.getNodeCount(d, subcord);
        for (int scNodeOrdinal = 0; scNodeOrdinal < scNodeCount; scNodeOrdinal++) {
          int scNode = (int)cellTopo.getNodeMap(d, subcord, scNodeOrdinal);
          if (nodesForSide.find(scNode) == nodesForSide.end()) {
            allNodesMatch = false;
            break;
          }
        }
      }
      if (allNodesMatch) {
        subcellEntities.push_back({d, subcord});
      }
    }
  }
  return subcellEntities;
}

//! Returns ordinals in the basis that have nodes on the specified side, as a sorted vector<int>
template <class Basis>
std::vector<int> localDofOrdinalsForSide(RCP<Basis> basis, int sideOrdinal) {
  using namespace std;
  // to use dof tags for this, we first need to determine the subcells of the domain that
  // are part of the specified side
  auto subcellEntities = subcellEntitiesForSide(basis->getBaseCellTopology(), sideOrdinal);

  auto dofOrdinalData = basis->getAllDofOrdinal();

  // determine size of first two parts of dofOrdinalData container
  // for extents > 0, there may no entries at all (for lower-order bases)
  int maxDim            = dofOrdinalData.extent(0);
  int maxSubcellOrdinal = dofOrdinalData.extent(1);

  vector<int> localDofOrdinals;

  for (auto subcellEntity : subcellEntities) {
    int subcellDim     = subcellEntity.first;
    int subcellOrdinal = subcellEntity.second;

    if (subcellDim >= maxDim) continue;  // no entries

    int dofContainerSize = dofOrdinalData.extent(2);  // 3rd extent: max dof count

    if (subcellOrdinal >= maxSubcellOrdinal) continue;  // no entries

    for (int entryOrdinal = 0; entryOrdinal < dofContainerSize; entryOrdinal++) {
      int localDofOrdinal = dofOrdinalData(subcellDim, subcellOrdinal, entryOrdinal);

      if (localDofOrdinal >= 0)
        localDofOrdinals.push_back(localDofOrdinal);
      else
        break;
    }
  }

  std::sort(localDofOrdinals.begin(), localDofOrdinals.end());

  return localDofOrdinals;
}

class Symmetries {
  // the symmetries on a given topology (e.g. cube) are a subset of the
  // permutations of the nodes: namely, the ones for which nodal connectivities
  // (edges) are preserved.
  // This class enumerates the symmetries.  The 0 symmetry is the identity permutation.
  vector<vector<int>> _symmetries;
  vector<vector<int>> _inverses;
  int _N;

 public:
  Symmetries(const vector<vector<int>> &symmetriesList) {
    _symmetries = symmetriesList;
    // sanity checks on the input:
    // - require that all entries are of equal length N
    // - require that all entries contain each of 0, ..., N-1

    if (symmetriesList.size() == 0) {
      _N = 0;
    } else {
      _N = int(symmetriesList[0].size());
      for (auto symmetry : symmetriesList) {
        TEUCHOS_TEST_FOR_EXCEPTION(int(symmetry.size()) != _N, std::invalid_argument, "Each symmetry must have the same length as all the others.");
        vector<int> inverse(_N, -1);
        for (int i = 0; i < _N; i++) {
          int i_mapped = symmetry[i];
          TEUCHOS_TEST_FOR_EXCEPTION((i_mapped < 0) || (i_mapped > _N), std::invalid_argument,
                                     "Each symmetry entry must be between 0 and N-1, inclusive.");
          inverse[i_mapped] = i;
        }
        for (int j = 0; j < _N; j++) {
          TEUCHOS_TEST_FOR_EXCEPTION(j == -1, std::invalid_argument, "Each symmetry must include every integer between 0 and N-1, inclusive.");
        }
        _inverses.push_back(inverse);
      }
    }
  }

  const vector<int> &getPermutation(int permutationOrdinal) {
    TEUCHOS_TEST_FOR_EXCEPTION(permutationOrdinal < 0, std::invalid_argument, "permutationOrdinal must be positive");
    TEUCHOS_TEST_FOR_EXCEPTION(permutationOrdinal >= int(_symmetries.size()), std::invalid_argument, "permutationOrdinal out of range");
    return _symmetries[permutationOrdinal];
  }

  int getPermutationCount() {
    return _symmetries.size();
  }

  void printAll() {
    for (auto permutation : _symmetries) {
      cout << permutationString(permutation) << endl;
    }
  }

  static string permutationString(vector<int> &permutation) {
    ostringstream permString;
    permString << "{ ";
    for (int entry : permutation) {
      permString << entry << " ";
    }
    permString << "}";
    return permString.str();
  }

  typedef vector<vector<int>> EdgeContainer;

  // assumes that the inner vectors are sorted:
  static bool hasEdge(const EdgeContainer &edges, int node1, int node2) {
    return std::find(edges[node1].begin(), edges[node1].end(), node2) != edges[node1].end();
  }

  /*
   For any topology, we can determine the symmetries from
   the graph of the edges.  This method takes as argument
   a representation of said graph.  The edges container has
   size equal to the node count, and the entries for a given
   node are exactly the nodes that it is connected to by an edge.
   */
  static Symmetries symmetries(const EdgeContainer &unsortedEdges) {
    EdgeContainer edges = unsortedEdges;
    // now sort
    int nodeCount = (int)edges.size();
    for (int node = 0; node < nodeCount; node++) {
      std::sort(edges[node].begin(), edges[node].end());
    }

    vector<vector<int>> permutations;

    // recursive lambda function to determine valid permutations that start with a specified set of choices:
    function<void(const vector<int> &)> addPermutationsThatMatch;
    addPermutationsThatMatch = [nodeCount, edges, &addPermutationsThatMatch, &permutations](const vector<int> &permutationStart) -> void {
      int nextNodeToMap = int(permutationStart.size());
      // iterate through the possible choices for nodes to map to:
      for (int possibleMappedNode = 0; possibleMappedNode < nodeCount; possibleMappedNode++) {
        // is this node already mapped?
        bool isAlreadyMapped = false;
        for (int mappedNode : permutationStart) {
          if (mappedNode == possibleMappedNode) isAlreadyMapped = true;
        }
        if (isAlreadyMapped) continue;
        // not already mapped: let's check whether the edge relationships agree
        bool mappingAgrees = true;
        int numNodesMapped = (int)permutationStart.size();
        for (int node = 0; node < numNodesMapped; node++) {
          int mappedNode       = permutationStart[node];
          bool originalHasEdge = hasEdge(edges, node, nextNodeToMap);
          bool mappedHasEdge   = hasEdge(edges, mappedNode, possibleMappedNode);
          if (originalHasEdge != mappedHasEdge) {
            mappingAgrees = false;
            break;
          }
        }
        if (mappingAgrees) {
          // no conflict detected: try this mapping choice
          vector<int> permutationStartCopy = permutationStart;
          permutationStartCopy.push_back(possibleMappedNode);
          if (int(permutationStartCopy.size()) == nodeCount) {
            permutations.push_back(permutationStartCopy);
          } else {
            addPermutationsThatMatch(permutationStartCopy);
          }
        }
      }
    };

    vector<int> permutation;
    addPermutationsThatMatch(permutation);

    return Symmetries(permutations);
  }

  static Symmetries shardsSymmetries(shards::CellTopology &cellTopo) {
    int edgeDim   = 1;
    int nodeCount = cellTopo.getNodeCount();
    TEUCHOS_TEST_FOR_EXCEPTION(nodeCount != int(cellTopo.getVertexCount()), std::invalid_argument, "Higher-order topologies are not supported");
    int edgeCount = cellTopo.getSubcellCount(edgeDim);
    EdgeContainer edges(nodeCount);  // set --> sorted, which is handy
    for (int edgeOrdinal = 0; edgeOrdinal < edgeCount; edgeOrdinal++) {
      int node0 = cellTopo.getNodeMap(edgeDim, edgeOrdinal, 0);
      int node1 = cellTopo.getNodeMap(edgeDim, edgeOrdinal, 1);
      edges[node0].push_back(node1);
      edges[node1].push_back(node0);
    }
    return symmetries(edges);
  }
};

class UniqueNumbering {
  vector<vector<double>> _knownCoords;  // x,y,z, depending on spatial extent; inner vector is sorted
  double _tol;                          // what counts as a match
  map<vector<double>, int> _numbering;  // maps from tuple selected from the sets in _knownCoords to unique identifier for tuple
  void getSanitizedCoords(const vector<double> &coords, vector<double> &sanitizedCoords) {
    sanitizedCoords.resize(coords.size());
    for (int d = 0; d < int(coords.size()); d++) {
      double val = coords[d];
      // look to lower_bound to see if we already have something within _tol
      auto lowerBoundIt = lower_bound(_knownCoords[d].begin(), _knownCoords[d].end(), val);
      // lower_bound returns iterator to first value that does not compare less than val (i.e. lower_bound >= val)
      if ((lowerBoundIt != _knownCoords[d].end()) && ((*lowerBoundIt - val) < _tol)) {
        sanitizedCoords[d] = *lowerBoundIt;
      } else {
        bool shouldInsert = true;  // unless we find that the previous entry is a match
        if (lowerBoundIt != _knownCoords[d].begin()) {
          // then check the prior guy to see if he's within _tol
          double previousEntry = *(lowerBoundIt - 1);
          if ((val - previousEntry) < _tol) {
            sanitizedCoords[d] = previousEntry;
            shouldInsert       = false;
          }
        }
        if (shouldInsert) {
          // neither prior nor following entry is within _tol
          // add to our list:
          _knownCoords[d].insert(lowerBoundIt, val);
          sanitizedCoords[d] = val;
        }
      }
    }
  }

 public:
  UniqueNumbering(int spaceDim, double tol) {
    _tol         = tol;
    _knownCoords = vector<vector<double>>(spaceDim);
  }
  template <class ArrayScalar, class ArrayOrdinal>
  void getIDs(const ArrayScalar &points, ArrayOrdinal &globalIDs) {
    auto points_host    = Kokkos::create_mirror_view(points);
    auto globalIDs_host = Kokkos::create_mirror_view(globalIDs);
    Kokkos::deep_copy(points_host, points);
    int spaceDim = _knownCoords.size();
    TEUCHOS_TEST_FOR_EXCEPTION(spaceDim != (int)points.extent(points.rank() - 1), std::invalid_argument, "final extent of points container must equal spaceDim");
    if (points.rank() == 2) {
      int numPoints = points.extent(0);
      for (int pointOrdinal = 0; pointOrdinal < numPoints; pointOrdinal++) {
        vector<double> coords(spaceDim);
        for (int d = 0; d < spaceDim; d++) {
          coords[d] = points_host(pointOrdinal, d);
        }
        globalIDs_host(pointOrdinal) = getGlobalID(coords);
      }
    } else if (points.rank() == 3) {
      int numCells  = points.extent(0);
      int numPoints = points.extent(1);
      for (int cellOrdinal = 0; cellOrdinal < numCells; cellOrdinal++) {
        for (int pointOrdinal = 0; pointOrdinal < numPoints; pointOrdinal++) {
          vector<double> coords(spaceDim);
          for (int d = 0; d < spaceDim; d++) {
            coords[d] = points_host(cellOrdinal, pointOrdinal, d);
          }
          globalIDs_host(cellOrdinal, pointOrdinal) = getGlobalID(coords);
        }
      }
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "points must be a rank 2 or rank 3 container");
    }
    Kokkos::deep_copy(globalIDs, globalIDs_host);
  }
  int getGlobalID(const vector<double> &coords) {
    TEUCHOS_TEST_FOR_EXCEPTION(coords.size() != _knownCoords.size(), std::invalid_argument,
                               "coords must be the same size as spaceDim!");
    // we first filter ("sanitize") the coords to match _knownCoords within _tol
    // (adding entries to _knownCoords as needed)
    vector<double> coordsKey;
    getSanitizedCoords(coords, coordsKey);

    // we look up using the sanitized coordsKey
    if (_numbering.find(coordsKey) == _numbering.end()) {
      int newID             = _numbering.size();
      _numbering[coordsKey] = newID;
    }
    return _numbering[coordsKey];
  }
  int totalCount() {
    return _numbering.size();
  }
};

using namespace std;
using namespace Teuchos;

template <class Basis, class ExecutionSpace, class ArrayScalar, class ArrayOrdinal>
void buildSampleElementToNodeMapThreeElements(RCP<Basis> basis, ArrayOrdinal &elemToNodeMap,
                                              vector<vector<int>> &ordinalsForSubcellDimension,
                                              vector<int> &subcellCountForDimension) {
  using namespace Kokkos::Experimental;
  typedef Intrepid2::CellTools<ExecutionSpace> CellTools;
  double pointTol = 1e-10;

  ArrayScalar cellWorkset;
  shards::CellTopology cellTopo = basis->getBaseCellTopology();
  int numCells                  = 3;
  int vertexCount               = (int)cellTopo.getVertexCount();
  int spaceDim                  = (int)cellTopo.getDimension();

  ArrayScalar refDofCoords, cellDofCoords;
  ArrayOrdinal cellDofIDs;  // for each cell, stores the global ordinal assigned to local ordinal
  resize(cellWorkset, numCells, vertexCount, spaceDim);
  resize(refDofCoords, basis->getCardinality(), spaceDim);
  resize(cellDofCoords, numCells, basis->getCardinality(), spaceDim);
  resize(cellDofIDs, numCells, basis->getCardinality());
  resize(elemToNodeMap, numCells, basis->getCardinality());

  basis->getDofCoords(refDofCoords);

  subcellCountForDimension.resize(spaceDim + 1);
  subcellCountForDimension[spaceDim] = numCells;
  auto cellWorkset_host              = Kokkos::create_mirror_view(cellWorkset);
  if (spaceDim == 1) {
    // line
    cellWorkset_host(0, 0, 0)   = 0.0;
    cellWorkset_host(0, 1, 0)   = 1.0;
    cellWorkset_host(1, 0, 0)   = 1.0;
    cellWorkset_host(1, 1, 0)   = 2.0;
    cellWorkset_host(2, 0, 0)   = 2.0;
    cellWorkset_host(2, 1, 0)   = 3.0;
    subcellCountForDimension[0] = 4;  // 3 cells x 2 vertices - (2 interior vertices)
  } else if ((spaceDim == 2) && (vertexCount == 4)) {
    // quad
    // first element: LL @ (0,1), UR @ (1,2)
    cellWorkset_host(0, 0, 0) = 0.0;
    cellWorkset_host(0, 0, 1) = 1.0;
    cellWorkset_host(0, 1, 0) = 1.0;
    cellWorkset_host(0, 1, 1) = 1.0;
    cellWorkset_host(0, 2, 0) = 1.0;
    cellWorkset_host(0, 2, 1) = 2.0;
    cellWorkset_host(0, 3, 0) = 0.0;
    cellWorkset_host(0, 3, 1) = 2.0;
    // second element: LL @ (0,0), UR @ (1,1)
    cellWorkset_host(1, 0, 0) = 0.0;
    cellWorkset_host(1, 0, 1) = 0.0;
    cellWorkset_host(1, 1, 0) = 1.0;
    cellWorkset_host(1, 1, 1) = 0.0;
    cellWorkset_host(1, 2, 0) = 1.0;
    cellWorkset_host(1, 2, 1) = 1.0;
    cellWorkset_host(1, 3, 0) = 0.0;
    cellWorkset_host(1, 3, 1) = 1.0;
    // third element: LL @ (1,0), UR @ (2,1)
    cellWorkset_host(2, 0, 0)   = 1.0;
    cellWorkset_host(2, 0, 1)   = 0.0;
    cellWorkset_host(2, 1, 0)   = 2.0;
    cellWorkset_host(2, 1, 1)   = 0.0;
    cellWorkset_host(2, 2, 0)   = 2.0;
    cellWorkset_host(2, 2, 1)   = 1.0;
    cellWorkset_host(2, 3, 0)   = 1.0;
    cellWorkset_host(2, 3, 1)   = 1.0;
    subcellCountForDimension[0] = 12 - 2 - 1 * 2;  // 3 cells x 4 vertices - (2 vertices seen by 2 cells, 1 vertex seen by 3 cells)
    subcellCountForDimension[1] = 12 - 2;          // 3 cells x 4 edges - (2 shared edges)
  } else if ((spaceDim == 3) && (vertexCount == 8)) {
    // hex: same geometry as quad, but extruded in z from 0 to 1
    // first element: LL @ (0,1), UR @ (1,2)
    cellWorkset_host(0, 0, 0) = 0.0;
    cellWorkset_host(0, 0, 1) = 1.0;
    cellWorkset_host(0, 0, 2) = 0.0;
    cellWorkset_host(0, 1, 0) = 1.0;
    cellWorkset_host(0, 1, 1) = 1.0;
    cellWorkset_host(0, 1, 2) = 0.0;
    cellWorkset_host(0, 2, 0) = 1.0;
    cellWorkset_host(0, 2, 1) = 2.0;
    cellWorkset_host(0, 2, 2) = 0.0;
    cellWorkset_host(0, 3, 0) = 0.0;
    cellWorkset_host(0, 3, 1) = 2.0;
    cellWorkset_host(0, 3, 2) = 0.0;
    cellWorkset_host(0, 4, 0) = 0.0;
    cellWorkset_host(0, 4, 1) = 1.0;
    cellWorkset_host(0, 4, 2) = 1.0;
    cellWorkset_host(0, 5, 0) = 1.0;
    cellWorkset_host(0, 5, 1) = 1.0;
    cellWorkset_host(0, 5, 2) = 1.0;
    cellWorkset_host(0, 6, 0) = 1.0;
    cellWorkset_host(0, 6, 1) = 2.0;
    cellWorkset_host(0, 6, 2) = 1.0;
    cellWorkset_host(0, 7, 0) = 0.0;
    cellWorkset_host(0, 7, 1) = 2.0;
    cellWorkset_host(0, 7, 2) = 1.0;
    // second element: LL @ (0,0), UR @ (1,1)
    cellWorkset_host(1, 0, 0) = 0.0;
    cellWorkset_host(1, 0, 1) = 0.0;
    cellWorkset_host(1, 0, 2) = 0.0;
    cellWorkset_host(1, 1, 0) = 1.0;
    cellWorkset_host(1, 1, 1) = 0.0;
    cellWorkset_host(1, 1, 2) = 0.0;
    cellWorkset_host(1, 2, 0) = 1.0;
    cellWorkset_host(1, 2, 1) = 1.0;
    cellWorkset_host(1, 2, 2) = 0.0;
    cellWorkset_host(1, 3, 0) = 0.0;
    cellWorkset_host(1, 3, 1) = 1.0;
    cellWorkset_host(1, 3, 2) = 0.0;
    cellWorkset_host(1, 4, 0) = 0.0;
    cellWorkset_host(1, 4, 1) = 0.0;
    cellWorkset_host(1, 4, 2) = 1.0;
    cellWorkset_host(1, 5, 0) = 1.0;
    cellWorkset_host(1, 5, 1) = 0.0;
    cellWorkset_host(1, 5, 2) = 1.0;
    cellWorkset_host(1, 6, 0) = 1.0;
    cellWorkset_host(1, 6, 1) = 1.0;
    cellWorkset_host(1, 6, 2) = 1.0;
    cellWorkset_host(1, 7, 0) = 0.0;
    cellWorkset_host(1, 7, 1) = 1.0;
    cellWorkset_host(1, 7, 2) = 1.0;
    // third element: LL @ (1,0), UR @ (2,1)
    cellWorkset_host(2, 0, 0)   = 1.0;
    cellWorkset_host(2, 0, 1)   = 0.0;
    cellWorkset_host(2, 0, 2)   = 0.0;
    cellWorkset_host(2, 1, 0)   = 2.0;
    cellWorkset_host(2, 1, 1)   = 0.0;
    cellWorkset_host(2, 1, 2)   = 0.0;
    cellWorkset_host(2, 2, 0)   = 2.0;
    cellWorkset_host(2, 2, 1)   = 1.0;
    cellWorkset_host(2, 2, 2)   = 0.0;
    cellWorkset_host(2, 3, 0)   = 1.0;
    cellWorkset_host(2, 3, 1)   = 1.0;
    cellWorkset_host(2, 3, 2)   = 0.0;
    cellWorkset_host(2, 4, 0)   = 1.0;
    cellWorkset_host(2, 4, 1)   = 0.0;
    cellWorkset_host(2, 4, 2)   = 1.0;
    cellWorkset_host(2, 5, 0)   = 2.0;
    cellWorkset_host(2, 5, 1)   = 0.0;
    cellWorkset_host(2, 5, 2)   = 1.0;
    cellWorkset_host(2, 6, 0)   = 2.0;
    cellWorkset_host(2, 6, 1)   = 1.0;
    cellWorkset_host(2, 6, 2)   = 1.0;
    cellWorkset_host(2, 7, 0)   = 1.0;
    cellWorkset_host(2, 7, 1)   = 1.0;
    cellWorkset_host(2, 7, 2)   = 1.0;
    subcellCountForDimension[0] = 24 - 4 * 1 - 2 * 2;  // 3 cells x 8 vertices - (4 vertices seen by 2 cells, 2 vertices seen by 3 cells)
    subcellCountForDimension[1] = 36 - 6 * 1 - 1 * 2;  // 3 cells x 12 edges - (6 edges seen by 2 cells, 1 edge seen by 3 cells)
    subcellCountForDimension[2] = 18 - 2 * 1;          // 3 cells x 6 faces - (2 faces seen by 2 cells)
  } else {
    // we haven't yet set up test geometry for this CellTopology
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unimplemented test case");
  }

  Kokkos::deep_copy(cellWorkset, cellWorkset_host);
  CellTools::mapToPhysicalFrame(cellDofCoords, refDofCoords, cellWorkset, cellTopo);

  Kokkos::fence();  // mapToPhysicalFrame calls getValues which calls kernels, so fence is required before UVM reads below

  UniqueNumbering globalNumbering(spaceDim, pointTol);
  globalNumbering.getIDs<ArrayScalar, ArrayOrdinal>(cellDofCoords, cellDofIDs);

  int tagOrdSubcellDim = 0;  // where to find the subcell extent within the tag
  int numDofsPerCell   = basis->getCardinality();

  // store ordinals in a set for easy uniquing
  vector<set<int>> ordinalsForSubcellDimensionSet(spaceDim + 1);

  auto cellDofIDs_host    = Kokkos::create_mirror_view(cellDofIDs);
  auto elemToNodeMap_host = Kokkos::create_mirror_view(elemToNodeMap);
  Kokkos::deep_copy(cellDofIDs_host, cellDofIDs);
  for (int cellOrdinal = 0; cellOrdinal < numCells; cellOrdinal++) {
    for (int dofOrdinal = 0; dofOrdinal < numDofsPerCell; dofOrdinal++) {
      int subcellDofDim = basis->getDofTag(dofOrdinal)[tagOrdSubcellDim];
      int globalID      = cellDofIDs_host(cellOrdinal, dofOrdinal);
      ordinalsForSubcellDimensionSet[subcellDofDim].insert(globalID);
      elemToNodeMap_host(cellOrdinal, dofOrdinal) = globalID;
    }
  }
  Kokkos::deep_copy(elemToNodeMap, elemToNodeMap_host);
  ordinalsForSubcellDimension.clear();
  for (int d = 0; d < spaceDim + 1; d++) {
    vector<int> ordinals(ordinalsForSubcellDimensionSet[d].begin(), ordinalsForSubcellDimensionSet[d].end());
    ordinalsForSubcellDimension.push_back(ordinals);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void testUniqueNumbering(Teuchos::FancyOStream &out, bool &success) {
  // simple test with dof coords corresponding to our sample three-element quad mesh with a quadratic basis
  typedef typename Node::device_type::execution_space ES;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef typename Teuchos::ScalarTraits<LocalOrdinal>::magnitudeType OT;  // ordinal type

  using namespace Kokkos::Experimental;
  typedef Kokkos::DynRankView<MT, typename Node::device_type> FC;
  typedef Kokkos::DynRankView<OT, typename Node::device_type> FCO;  // FC of ordinals
  typedef Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES, MT, MT> Basis;

  typename FC::HostMirror physDofCoords;
  typename FCO::HostMirror cellGIDs;
  std::vector<std::vector<LocalOrdinal>> expectedGIDs = {{0, 1, 2, 3, 4, 5, 6, 7, 8},
                                                         {9, 10, 11, 12, 13, 14, 0, 1, 2},
                                                         {11, 15, 16, 14, 17, 18, 2, 19, 20}};
  int numCells = 3, numPointsPerCell = 9, spaceDim = 2;
  resize(physDofCoords, numCells, numPointsPerCell, spaceDim);
  resize(cellGIDs, numCells, numPointsPerCell);

  // x,y coords of lower-left vertices for each cell:
  std::vector<double> x0s = {0.0, 0.0, 1.0};
  std::vector<double> y0s = {1.0, 0.0, 0.0};
  for (int cellOrdinal = 0; cellOrdinal < numCells; cellOrdinal++) {
    double x0 = x0s[cellOrdinal], y0 = y0s[cellOrdinal];
    for (int pointOrdinal = 0; pointOrdinal < numPointsPerCell; pointOrdinal++) {
      double xOffset                              = (pointOrdinal % 3) * 0.5;  // 0, 0.5, 1.0
      double yOffset                              = (pointOrdinal / 3) * 0.5;  // 0, 0.5, 1.0
      physDofCoords(cellOrdinal, pointOrdinal, 0) = x0 + xOffset;
      physDofCoords(cellOrdinal, pointOrdinal, 1) = y0 + yOffset;
    }
  }
  double pointTol = 1e-10;
  UniqueNumbering numbering(spaceDim, pointTol);
  numbering.getIDs(physDofCoords, cellGIDs);

  for (int cellOrdinal = 0; cellOrdinal < numCells; cellOrdinal++) {
    out << "cell " << cellOrdinal << " GIDs: {";
    for (int pointOrdinal = 0; pointOrdinal < numPointsPerCell; pointOrdinal++) {
      out << " " << cellGIDs(cellOrdinal, pointOrdinal);
      if (cellGIDs(cellOrdinal, pointOrdinal) != expectedGIDs[cellOrdinal][pointOrdinal]) {
        success = false;
      }
    }
    out << " }\n";
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, UniqueNumbering, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  ;
  testUniqueNumbering<Scalar, LocalOrdinal, GlobalOrdinal, Node>(out, success);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void testBuildSampleElementToNodeMapThreeElementQuad(Teuchos::FancyOStream &out, bool &success) {
  // simple test with quadratic basis on quads: check that we have the right numbering (hard-coded)
  std::vector<int> vertexGIDs = {0, 2, 6, 8, 9, 11, 16, 20};           // should be 8 of these
  std::vector<int> edgeGIDs   = {1, 3, 5, 7, 10, 12, 14, 15, 18, 19};  // 10 of these
  std::vector<int> cellGIDs   = {4, 13, 17};

  vector<vector<int>> expectedGIDs = {vertexGIDs, edgeGIDs, cellGIDs};

  typedef typename Node::device_type::execution_space ES;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef typename Teuchos::ScalarTraits<LocalOrdinal>::magnitudeType OT;  // ordinal type

  using namespace Kokkos::Experimental;
  typedef Kokkos::DynRankView<MT, typename Node::device_type> FC;
  typedef Kokkos::DynRankView<OT, typename Node::device_type> FCO;  // FC of ordinals
  typedef Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES, MT, MT> Basis;

  int polyOrder = 2, spaceDim = 2;
  RCP<Basis> basis = rcp(new Basis(polyOrder, Intrepid2::EPointType::POINTTYPE_EQUISPACED));
  std::vector<int> subcellCountForDimension;
  std::vector<std::vector<int>> ordinalsForSubcellDimension;
  FCO elementToNodeMap;
  buildSampleElementToNodeMapThreeElements<Basis, ES, FC, FCO>(basis, elementToNodeMap, ordinalsForSubcellDimension,
                                                               subcellCountForDimension);
  // check subcell counts first
  for (int d = 0; d <= spaceDim; d++) {
    if (subcellCountForDimension[d] != int(expectedGIDs[d].size())) {
      success = false;
      out << "subcellCountForDimension[" << d << "] != " << expectedGIDs[d].size() << endl;
    } else {
      int dofCount     = ordinalsForSubcellDimension[d].size();
      bool listsDiffer = ordinalsForSubcellDimension[d].size() != expectedGIDs[d].size();

      if (!listsDiffer) {
        for (int dofOrdinal = 0; dofOrdinal < dofCount; dofOrdinal++) {
          if (ordinalsForSubcellDimension[d][dofOrdinal] != expectedGIDs[d][dofOrdinal]) {
            listsDiffer = true;
            break;
          }
        }
      }
      if (listsDiffer) {
        success = false;
        out << "lists for d=" << d << " differ; expected: {";
        for (LocalOrdinal LID : expectedGIDs[d]) {
          out << " " << LID;
        }
        out << " }; actual: {";
        for (LocalOrdinal LID : ordinalsForSubcellDimension[d]) {
          out << " " << LID;
        }
        out << " }\n";
      }
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, BuildSampleElementToNodeMapThreeElementQuad, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  ;
  testBuildSampleElementToNodeMapThreeElementQuad<Scalar, LocalOrdinal, GlobalOrdinal, Node>(out, success);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node, class Basis, class ExecutionSpace, class ArrayScalar, class ArrayOrdinal>
void testFindSeeds(Xpetra::UnderlyingLib lib, int max_degree, Intrepid2::EPointType ptype, int numRanks, Teuchos::FancyOStream &out, bool &success) {
  // "numRanks" is for emulated parallel execution: we set up several maps, one corresponding to each emulated MPI rank.

  typedef Intrepid2::CellTools<ExecutionSpace> CellTools;
  typedef ArrayScalar FC;
  typedef ArrayOrdinal FCO;

  FCO elementToNodeMap;
  vector<vector<int>> ordinalsForSubcellDimension;
  vector<int> subcellCountForDimension;

  using namespace Teuchos;
  RCP<Comm<int>> serialComm = rcp(new SerialComm<int>());
  typedef Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node> MapFactory;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> serialMapRCP, rowMapRCP, colMapRCP;

  for (int polyOrder = 1; polyOrder < max_degree; polyOrder++) {
    RCP<Basis> basis              = rcp(new Basis(polyOrder, ptype));
    shards::CellTopology cellTopo = basis->getBaseCellTopology();
    int spaceDim                  = cellTopo.getDimension();
    buildSampleElementToNodeMapThreeElements<Basis, ExecutionSpace, FC, FCO>(basis, elementToNodeMap, ordinalsForSubcellDimension,
                                                                             subcellCountForDimension);
    vector<vector<LocalOrdinal>> seeds;

    int numCells               = elementToNodeMap.extent(0);
    int dofsPerCell            = elementToNodeMap.extent(1);
    int maxGID                 = -1;
    auto elementToNodeMap_host = Kokkos::create_mirror_view(elementToNodeMap);
    Kokkos::deep_copy(elementToNodeMap_host, elementToNodeMap);
    for (int cellOrdinal = 0; cellOrdinal < numCells; cellOrdinal++) {
      for (int dofOrdinal = 0; dofOrdinal < dofsPerCell; dofOrdinal++) {
        maxGID = max(maxGID, elementToNodeMap_host(cellOrdinal, dofOrdinal));
      }
    }
    int numElements = maxGID + 1;
    // construct a serial map (claiming all rows and cells)
    serialMapRCP = MapFactory::createLocalMapWithNode(lib, numElements, serialComm);

    int startingGID = 0;
    for (int rank = 0; rank < numRanks; rank++) {
      // simple distribution: take numElements / numRanks on each rank
      int numRankLocalElements = numElements / numRanks;
      int remainder            = numElements - (numElements / numRanks) * numRanks;
      if (remainder > rank) numRankLocalElements++;
      FCO rankLocalElementToNodeMap;
      resize(rankLocalElementToNodeMap, numCells, basis->getCardinality());
      auto rankLocalElementToNodeMap_host = Kokkos::create_mirror_view(rankLocalElementToNodeMap);
      vector<set<LocalOrdinal>> expectedSeedsSets(spaceDim + 1);

      auto isRankLocal = [startingGID, numRankLocalElements](GlobalOrdinal GID) -> bool {
        return (GID >= startingGID) && (GID < startingGID + numRankLocalElements);
      };

      Teuchos::Array<GlobalOrdinal> myRowGIDs(numRankLocalElements);
      for (int i = 0; i < numRankLocalElements; i++) {
        myRowGIDs[i] = startingGID + i;
      }
      // colMap sees all local rows, plus any non-local IDs they talk to:
      Teuchos::Array<GlobalOrdinal> myColGIDs = myRowGIDs;
      set<GlobalOrdinal> offRankGIDs;
      for (int cellOrdinal = 0; cellOrdinal < numCells; cellOrdinal++) {
        // first pass: does this cell contain any GIDs we own?
        bool hasOwnedGIDs = false;
        for (int dofOrdinal = 0; dofOrdinal < dofsPerCell; dofOrdinal++) {
          GlobalOrdinal GID = elementToNodeMap_host(cellOrdinal, dofOrdinal);
          if (isRankLocal(GID)) {
            hasOwnedGIDs = true;
          }
        }
        // second pass: if there are owned GIDs, add any off-rank guys this cell sees
        if (hasOwnedGIDs) {
          for (int dofOrdinal = 0; dofOrdinal < dofsPerCell; dofOrdinal++) {
            GlobalOrdinal GID = elementToNodeMap_host(cellOrdinal, dofOrdinal);
            if (!isRankLocal(GID)) {
              if (GID == -1) {
                cout << "hmmm....\n";
              }
              offRankGIDs.insert(GID);
            }
          }
        }
      }
      // add off-rank GIDs to the end of the colMapGIDs:
      myColGIDs.insert(myColGIDs.end(), offRankGIDs.begin(), offRankGIDs.end());

      GlobalOrdinal indexBase  = 0;
      GlobalOrdinal GO_INVALID = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();
      rowMapRCP                = MapFactory::Build(lib, GO_INVALID, myRowGIDs().getConst(), indexBase, serialComm);
      colMapRCP                = MapFactory::Build(lib, GO_INVALID, myColGIDs().getConst(), indexBase, serialComm);

      // rewrite elementToNodeMap to contain only LIDs, and determine expected seeds
      for (int cellOrdinal = 0; cellOrdinal < numCells; cellOrdinal++) {
        vector<vector<LocalOrdinal>> seedForSubcell(spaceDim + 1);
        for (int d = 0; d <= spaceDim; d++) {
          int subcellCount  = cellTopo.getSubcellCount(d);
          seedForSubcell[d] = vector<LocalOrdinal>(subcellCount, LocalOrdinal(-1));
        }
        for (int dofOrdinal = 0; dofOrdinal < dofsPerCell; dofOrdinal++) {
          GlobalOrdinal GID                                       = elementToNodeMap_host(cellOrdinal, dofOrdinal);
          LocalOrdinal LID                                        = colMapRCP->getLocalElement(GID);
          rankLocalElementToNodeMap_host(cellOrdinal, dofOrdinal) = LID;  // may be -1
          if (LID != -1) {
            int subcdim               = basis->getDofTag(dofOrdinal)[0];
            int subcord               = basis->getDofTag(dofOrdinal)[1];
            LocalOrdinal existingSeed = seedForSubcell[subcdim][subcord];
            if ((existingSeed == (LocalOrdinal)-1) || (LID < existingSeed)) {
              seedForSubcell[subcdim][subcord] = LID;
            }
          }
        }
        for (int d = 0; d <= spaceDim; d++) {
          for (LocalOrdinal seed : seedForSubcell[d]) {
            if (seed != -1) {
              GlobalOrdinal GID = colMapRCP->getGlobalElement(seed);
              if (rowMapRCP->getLocalElement(GID) == seed) {
                expectedSeedsSets[d].insert(seed);
              }
            }
          }
        }
      }

      Kokkos::deep_copy(rankLocalElementToNodeMap, rankLocalElementToNodeMap_host);
      MueLu::MueLuIntrepid::FindGeometricSeedOrdinals<Basis, FCO, LocalOrdinal, GlobalOrdinal, Node>(basis, rankLocalElementToNodeMap,
                                                                                                     seeds, *rowMapRCP, *colMapRCP);

      if (int(seeds.size()) != spaceDim + 1) {
        success = false;
        out << "seeds should have extent = spaceDim + 1; ";
        out << seeds.size() << " != " << spaceDim + 1 << endl;
      } else {
        for (int d = 0; d <= spaceDim; d++) {
          // check that we have exactly one entry for each subcell entity
          int expectedSeedCount = expectedSeedsSets[d].size();
          int seedCount         = seeds[d].size();
          if (expectedSeedCount != seedCount) {
            success = false;
            out << "expected " << expectedSeedCount << " seeds for extent " << d;
            out << "; had " << seedCount << endl;
          }
          // check that each entry belongs to an entity of the correct type
          for (LocalOrdinal localDofOrdinal : seeds[d]) {
            if (expectedSeedsSets[d].find(localDofOrdinal) == expectedSeedsSets[d].end()) {
              success = false;
              out << "Found local dof ordinal " << localDofOrdinal << " as a seed for extent ";
              out << d << ", but did not find it in expectedSeedsSets[" << d << "]\n";
            }
          }
        }
      }
      startingGID += numRankLocalElements;
    }
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node, class Basis, class ExecutionSpace, class ArrayScalar, class ArrayOrdinal>
void testFindSeedsSerial(Xpetra::UnderlyingLib lib, int max_degree, Intrepid2::EPointType ptype, Teuchos::FancyOStream &out, bool &success) {
  typedef Intrepid2::CellTools<ExecutionSpace> CellTools;
  typedef ArrayScalar FC;
  typedef ArrayOrdinal FCO;

  FCO elementToNodeMap;
  vector<vector<int>> ordinalsForSubcellDimension;
  vector<int> subcellCountForDimension;

  using namespace Teuchos;
  RCP<Comm<int>> serialComm = rcp(new SerialComm<int>());
  typedef Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node> MapFactory;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> serialMapRCP;

  for (int polyOrder = 1; polyOrder < max_degree; polyOrder++) {
    RCP<Basis> basis = rcp(new Basis(polyOrder, ptype));
    buildSampleElementToNodeMapThreeElements<Basis, ExecutionSpace, FC, FCO>(basis, elementToNodeMap, ordinalsForSubcellDimension,
                                                                             subcellCountForDimension);
    vector<vector<LocalOrdinal>> seeds;

    // construct a serial map (claiming all rows and cells)
    int numCells               = elementToNodeMap.extent(0);
    int dofsPerCell            = elementToNodeMap.extent(1);
    int maxLID                 = -1;
    auto elementToNodeMap_host = Kokkos::create_mirror_view(elementToNodeMap);
    Kokkos::deep_copy(elementToNodeMap_host, elementToNodeMap);
    for (int cellOrdinal = 0; cellOrdinal < numCells; cellOrdinal++) {
      for (int dofOrdinal = 0; dofOrdinal < dofsPerCell; dofOrdinal++) {
        maxLID = max(maxLID, elementToNodeMap_host(cellOrdinal, dofOrdinal));
      }
    }
    int numElements = maxLID + 1;
    serialMapRCP    = MapFactory::createLocalMapWithNode(lib, numElements, serialComm);

    MueLu::MueLuIntrepid::FindGeometricSeedOrdinals<Basis, FCO, LocalOrdinal, GlobalOrdinal, Node>(basis, elementToNodeMap,
                                                                                                   seeds, *serialMapRCP, *serialMapRCP);

    int spaceDim = basis->getBaseCellTopology().getDimension();
    if (int(seeds.size()) != spaceDim + 1) {
      success = false;
      out << "seeds should have extent = spaceDim + 1; ";
      out << seeds.size() << " != " << spaceDim + 1 << endl;
    } else {
      for (int d = 0; d <= spaceDim; d++) {
        // check that we have exactly one entry for each subcell entity
        int expectedSeedCount = subcellCountForDimension[d];
        if (basis->getDofCount(d, 0) == 0) {
          // no dofs for first subcell of extent d; we assume that all other subcells of extent d also
          // don't have dofs assigned
          expectedSeedCount = 0;
        }
        int seedCount = seeds[d].size();
        if (expectedSeedCount != seedCount) {
          success = false;
          out << "expected " << expectedSeedCount << " seeds for extent " << d;
          out << "; had " << seedCount << endl;
        }
        // check that each entry belongs to an entity of the correct type
        for (LocalOrdinal localDofOrdinal : seeds[d]) {
          if (std::find(ordinalsForSubcellDimension[d].begin(), ordinalsForSubcellDimension[d].end(), localDofOrdinal) == ordinalsForSubcellDimension[d].end()) {
            success = false;
            out << "Found local dof ordinal " << localDofOrdinal << " as a seed for extent ";
            out << d << ", but did not find it in ordinalsForSubcellDimension[" << d << "]\n";
          }
        }
      }
    }
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node, class Basis, class ExecutionSpace, class ArrayScalar, class ArrayOrdinal>
void testFindSeedsParallel(Xpetra::UnderlyingLib lib, int max_degree, Intrepid2::EPointType ptype, Teuchos::FancyOStream &out, bool &success) {
  typedef Intrepid2::CellTools<ExecutionSpace> CellTools;
  typedef ArrayScalar FC;
  typedef ArrayOrdinal FCO;

  for (int rankCount = 1; rankCount <= MAX_RANK_COUNT; rankCount++) {
    out << "running testFindSeedsParallel on " << rankCount << " emulated MPI ranks\n";
    testFindSeeds<LocalOrdinal, GlobalOrdinal, Node, Basis, ExecutionSpace, ArrayScalar, ArrayOrdinal>(lib, max_degree, ptype, rankCount, out, success);
  }
}

/******* End helper methods and classes by Nate ********/

/******* Begin typedefs for FindSeeds tests by Nate ********/
#define FIND_SEEDS_MACROS                                                    \
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;          \
  typedef typename Node::device_type::execution_space ES;                    \
  typedef Kokkos::DynRankView<MT, typename Node::device_type> FC;            \
  typedef Kokkos::DynRankView<LocalOrdinal, typename Node::device_type> FCO; \
  typedef Intrepid2::Basis_HGRAD_LINE_Cn_FEM<ES, MT, MT> LineBasis;          \
  typedef Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES, MT, MT> QuadBasis;          \
  typedef Intrepid2::Basis_HGRAD_HEX_Cn_FEM<ES, MT, MT> HexBasis;

/******* End typedefs for FindSeeds tests by Nate ********/

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, FindSeeds_Equispaced_Line, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  FIND_SEEDS_MACROS;
  testFindSeedsParallel<LocalOrdinal, GlobalOrdinal, Node, LineBasis, ES, FC, FCO>(MueLuTests::TestHelpers::Parameters::getLib(), MAX_LINE_DEGREE, Intrepid2::POINTTYPE_EQUISPACED, out, success);
  testFindSeedsSerial<LocalOrdinal, GlobalOrdinal, Node, LineBasis, ES, FC, FCO>(MueLuTests::TestHelpers::Parameters::getLib(), MAX_LINE_DEGREE, Intrepid2::POINTTYPE_EQUISPACED, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, FindSeeds_Spectral_Line, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  FIND_SEEDS_MACROS;
  const Intrepid2::EPointType POINTTYPE_SPECTRAL = static_cast<Intrepid2::EPointType>(1);  // Not sure why I have to do this...
  testFindSeedsParallel<LocalOrdinal, GlobalOrdinal, Node, LineBasis, ES, FC, FCO>(MueLuTests::TestHelpers::Parameters::getLib(), MAX_LINE_DEGREE, POINTTYPE_SPECTRAL, out, success);
  testFindSeedsSerial<LocalOrdinal, GlobalOrdinal, Node, LineBasis, ES, FC, FCO>(MueLuTests::TestHelpers::Parameters::getLib(), MAX_LINE_DEGREE, POINTTYPE_SPECTRAL, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, FindSeeds_Equispaced_Quad, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  FIND_SEEDS_MACROS;
  testFindSeedsParallel<LocalOrdinal, GlobalOrdinal, Node, QuadBasis, ES, FC, FCO>(MueLuTests::TestHelpers::Parameters::getLib(), MAX_QUAD_DEGREE, Intrepid2::POINTTYPE_EQUISPACED, out, success);
  testFindSeedsSerial<LocalOrdinal, GlobalOrdinal, Node, QuadBasis, ES, FC, FCO>(MueLuTests::TestHelpers::Parameters::getLib(), MAX_QUAD_DEGREE, Intrepid2::POINTTYPE_EQUISPACED, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, FindSeeds_Spectral_Quad, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  FIND_SEEDS_MACROS;
  const Intrepid2::EPointType POINTTYPE_SPECTRAL = static_cast<Intrepid2::EPointType>(1);  // Not sure why I have to do this...
  testFindSeedsParallel<LocalOrdinal, GlobalOrdinal, Node, QuadBasis, ES, FC, FCO>(MueLuTests::TestHelpers::Parameters::getLib(), MAX_QUAD_DEGREE, POINTTYPE_SPECTRAL, out, success);
  testFindSeedsSerial<LocalOrdinal, GlobalOrdinal, Node, QuadBasis, ES, FC, FCO>(MueLuTests::TestHelpers::Parameters::getLib(), MAX_QUAD_DEGREE, POINTTYPE_SPECTRAL, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, FindSeeds_Equispaced_Hex, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  FIND_SEEDS_MACROS;
  testFindSeedsParallel<LocalOrdinal, GlobalOrdinal, Node, HexBasis, ES, FC, FCO>(MueLuTests::TestHelpers::Parameters::getLib(), MAX_QUAD_DEGREE, Intrepid2::POINTTYPE_EQUISPACED, out, success);
  testFindSeedsSerial<LocalOrdinal, GlobalOrdinal, Node, HexBasis, ES, FC, FCO>(MueLuTests::TestHelpers::Parameters::getLib(), MAX_QUAD_DEGREE, Intrepid2::POINTTYPE_EQUISPACED, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, FindSeeds_Spectral_Hex, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  FIND_SEEDS_MACROS;
  const Intrepid2::EPointType POINTTYPE_SPECTRAL = static_cast<Intrepid2::EPointType>(1);  // Not sure why I have to do this...
  testFindSeedsParallel<LocalOrdinal, GlobalOrdinal, Node, HexBasis, ES, FC, FCO>(MueLuTests::TestHelpers::Parameters::getLib(), MAX_HEX_DEGREE, POINTTYPE_SPECTRAL, out, success);
  testFindSeedsSerial<LocalOrdinal, GlobalOrdinal, Node, HexBasis, ES, FC, FCO>(MueLuTests::TestHelpers::Parameters::getLib(), MAX_HEX_DEGREE, POINTTYPE_SPECTRAL, out, success);
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, GetP1NodeInHi, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType MT;
  typedef Kokkos::DynRankView<MT, typename Node::device_type> FC;
  typedef typename Node::device_type::execution_space ES;
  typedef Intrepid2::Basis<ES, MT, MT> Basis;

  out << "version: " << MueLu::Version() << std::endl;

  int max_degree = (Intrepid2::Parameters::MaxOrder < 5) ? Intrepid2::Parameters::MaxOrder : 5;

  {
    // QUAD
    RCP<Basis> lo = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<ES, MT, MT>());

    for (int i = 0; i < max_degree; i++) {
      RCP<Basis> hi = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES, MT, MT>(i, Intrepid2::POINTTYPE_EQUISPACED));
      std::vector<size_t> lo_node_in_hi;
      FC hi_dofCoords;

      MueLu::MueLuIntrepid::IntrepidGetP1NodeInHi<MT, typename Node::device_type>(hi, lo_node_in_hi, hi_dofCoords);
      TEST_EQUALITY((size_t)hi_dofCoords.extent(0), (size_t)hi->getCardinality());
      TEST_EQUALITY((size_t)hi_dofCoords.extent(1), (size_t)hi->getBaseCellTopology().getDimension());
      TEST_EQUALITY(lo_node_in_hi.size(), (size_t)lo->getCardinality());
    }
  }
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, BasisFactory, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType MT;
  typedef typename Node::device_type::execution_space ES;

  out << "version: " << MueLu::Version() << std::endl;

  // QUAD
  int degree;
  {
    bool test = rcp_dynamic_cast<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<ES, MT, MT>>(MueLu::MueLuIntrepid::BasisFactory<MT, ES>("hgrad_quad_c1", degree)) != Teuchos::null;
    TEST_EQUALITY(test, true);
  }
  {
    bool test = rcp_dynamic_cast<Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES, MT, MT>>(MueLu::MueLuIntrepid::BasisFactory<MT, ES>("hgrad_quad_c2", degree)) != Teuchos::null;
    TEST_EQUALITY(test, true);
  }
  {
    bool test = rcp_dynamic_cast<Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES, MT, MT>>(MueLu::MueLuIntrepid::BasisFactory<MT, ES>("hgrad_quad_c3", degree)) != Teuchos::null;
    TEST_EQUALITY(test, true);
  }
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, BuildLoElemToNode, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType MT;
  typedef Kokkos::DynRankView<MT, typename Node::device_type> FC;
  typedef Kokkos::DynRankView<LocalOrdinal, typename Node::device_type> FCi;
  typedef typename Node::device_type::execution_space ES;
  typedef Intrepid2::Basis<ES, MT, MT> Basis;

  out << "version: " << MueLu::Version() << std::endl;
  int max_degree = (Intrepid2::Parameters::MaxOrder < 5) ? Intrepid2::Parameters::MaxOrder : 5;

  {
    // QUAD
    //  A one element test with Kirby-numbered nodes where the top edge is not owned
    RCP<Basis> lo = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<ES, MT, MT>());

    for (int degree = 2; degree < max_degree; degree++) {
      int Nn        = (degree + 1) * (degree + 1);
      RCP<Basis> hi = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES, MT, MT>(degree, Intrepid2::POINTTYPE_EQUISPACED));
      FCi hi_e2n("hi_e2n", 1, Nn), lo_e2n;

      std::vector<bool> hi_owned(Nn, false), lo_owned;
      std::vector<size_t> lo_node_in_hi;
      std::vector<LO> hi_to_lo_map;
      int lo_numOwnedNodes = 0;
      FC hi_dofCoords;
      MueLu::MueLuIntrepid::IntrepidGetP1NodeInHi<MT, typename Node::device_type>(hi, lo_node_in_hi, hi_dofCoords);

      Kokkos::parallel_for(
          "IntrepidPCoarsenFactory,BuildLoElemToNode", Kokkos::RangePolicy<typename Node::device_type::execution_space>(0, Nn), KOKKOS_LAMBDA(int i) {
            hi_e2n(0, i) = i;
          });
      Kokkos::fence();
      for (int i = 0; i < Nn; i++) {
        if (i < Nn - (degree + 1)) hi_owned[i] = true;
      }

      Teuchos::ArrayRCP<const int> is_dirichlet(Nn, false);

      MueLu::MueLuIntrepid::BuildLoElemToNode(hi_e2n, hi_owned, lo_node_in_hi, is_dirichlet, lo_e2n, lo_owned, hi_to_lo_map, lo_numOwnedNodes);

      // Checks
      TEST_EQUALITY(lo_numOwnedNodes, 2);

      size_t num_lo_nodes_located = 0;
      for (size_t i = 0; i < hi_to_lo_map.size(); i++) {
        if (hi_to_lo_map[i] != Teuchos::OrdinalTraits<LO>::invalid())
          num_lo_nodes_located++;
      }
      TEST_EQUALITY(lo_owned.size(), num_lo_nodes_located);
    }
  }  // end QUAD
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, BuildLoElemToNodeWithDirichlet, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType MT;
  typedef Kokkos::DynRankView<MT, typename Node::device_type> FC;
  typedef Kokkos::DynRankView<LocalOrdinal, typename Node::device_type> FCi;
  typedef typename Node::device_type::execution_space ES;
  typedef Intrepid2::Basis<ES, MT, MT> Basis;

  out << "version: " << MueLu::Version() << std::endl;
  int max_degree = (Intrepid2::Parameters::MaxOrder < 5) ? Intrepid2::Parameters::MaxOrder : 5;

  {
    // QUAD
    //  A one element test with Kirby-numbered nodes where the top edge is not owned
    RCP<Basis> lo = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<ES, MT, MT>());

    for (int degree = 2; degree < max_degree; degree++) {
      int Nn        = (degree + 1) * (degree + 1);
      RCP<Basis> hi = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES, MT, MT>(degree, Intrepid2::POINTTYPE_EQUISPACED));
      FCi hi_e2n("hi_e2n", 1, Nn), lo_e2n;

      std::vector<bool> hi_owned(Nn, false), lo_owned;
      std::vector<size_t> lo_node_in_hi;
      std::vector<LO> hi_to_lo_map;
      int lo_numOwnedNodes = 0;
      FC hi_dofCoords;
      MueLu::MueLuIntrepid::IntrepidGetP1NodeInHi<MT, typename Node::device_type>(hi, lo_node_in_hi, hi_dofCoords);

      Kokkos::parallel_for(
          "IntrepidPCoarsenFactory,BuildLoElemToNodeWithDirichlet",
          Kokkos::RangePolicy<typename Node::device_type::execution_space>(0, Nn), KOKKOS_LAMBDA(int i) {
            hi_e2n(0, i) = i;
          });
      Kokkos::fence();
      for (int i = 0; i < Nn; i++) {
        if (i < Nn - (degree + 1)) hi_owned[i] = true;
      }

      // Mark of the very first node as Dirichlet
      Teuchos::ArrayRCP<int> is_dirichlet(Nn, false);
      is_dirichlet[0] = 1;

      MueLu::MueLuIntrepid::BuildLoElemToNode(hi_e2n, hi_owned, lo_node_in_hi, is_dirichlet, lo_e2n, lo_owned, hi_to_lo_map, lo_numOwnedNodes);

      // Checks
      TEST_EQUALITY(lo_numOwnedNodes, 1);

      size_t num_lo_nodes_located = 0;
      for (size_t i = 0; i < hi_to_lo_map.size(); i++) {
        if (hi_to_lo_map[i] != Teuchos::OrdinalTraits<LO>::invalid())
          num_lo_nodes_located++;
      }
      TEST_EQUALITY(lo_owned.size(), num_lo_nodes_located);
    }
  }  // end QUAD
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, GenerateColMapFromImport, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef GlobalOrdinal GO;
  typedef LocalOrdinal LO;

  out << "version: " << MueLu::Version() << std::endl;
  Xpetra::UnderlyingLib lib          = MueLuTests::TestHelpers::Parameters::getLib();
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();
  GO gst_invalid                     = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
  GO lo_invalid                      = Teuchos::OrdinalTraits<LO>::invalid();
  int NumProcs                       = comm->getSize();
  // This function basically takes an existing domain->column importer (the "hi order" guy) and using the hi_to_lo_map,
  // generates "lo order" version.  The domain map is already given to us here, so we just need to make tha column one.

  // We'll test this by starting with some linearly distributed map for the "hi" domain map and a duplicated map as the "hi" column map.
  // We then do a "lo" domain map, which just grabs the first GID on each proc.
  GO numGlobalElements = 100;
  Teuchos::Array<GO> hi_Cols(numGlobalElements);
  for (size_t i = 0; i < (size_t)numGlobalElements; i++)
    hi_Cols[i] = i;

  RCP<Map> hi_domainMap   = MapFactory::Build(lib, numGlobalElements, 0, comm);
  RCP<Map> hi_colMap      = MapFactory::Build(lib, gst_invalid, hi_Cols(), 0, comm);
  RCP<Import> hi_importer = ImportFactory::Build(hi_domainMap, hi_colMap);

  Teuchos::Array<GO> lo_Doms(1);
  lo_Doms[0]            = hi_domainMap->getGlobalElement(0);
  RCP<Map> lo_domainMap = MapFactory::Build(lib, gst_invalid, lo_Doms(), 0, comm);

  // Get the list of the first GID from each proc (aka the lo_colMap)
  std::vector<GO> send_colgids(1);
  std::vector<GO> recv_colgids(NumProcs);
  send_colgids[0] = hi_domainMap->getGlobalElement(0);
  comm->gatherAll(1 * sizeof(GO), (char *)send_colgids.data(), NumProcs * sizeof(GO), (char *)recv_colgids.data());

  // Build the h2l map
  std::vector<LO> hi_to_lo_map(hi_colMap->getLocalNumElements(), lo_invalid);
  for (size_t i = 0; i < (size_t)NumProcs; i++)
    hi_to_lo_map[recv_colgids[i]] = i;

  // Import
  size_t lo_columnMapLength = NumProcs;  // One col per proc
  RCP<const Map> lo_colMap;
  MueLu::MueLuIntrepid::GenerateColMapFromImport(*hi_importer, hi_to_lo_map, *lo_domainMap, lo_columnMapLength, lo_colMap);

  // Final test
  for (size_t i = 0; i < lo_colMap->getLocalNumElements(); i++)
    TEST_EQUALITY(recv_colgids[i], lo_colMap->getGlobalElement(i));
}

/*********************************************************************************************************************/
/* How this guy works:
   num_p1_nodes - number of nodes in the p=1 mesh
   p1_gold_in   - input vector for the lo_basis
   p2_gold_in   - output vector of linear interpolation from lo_basis to hi_basis

*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TestPseudoPoisson(Teuchos::FancyOStream &out, int num_p1_nodes, int degree, std::vector<Scalar> &lo_gold_in, std::vector<Scalar> &hi_gold_out, const std::string &hi_basis, const std::string lo_basis = "hgrad_line_c1") {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef Scalar SC;
  typedef GlobalOrdinal GO;
  typedef LocalOrdinal LO;
  typedef Node NO;
  typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef Kokkos::DynRankView<LocalOrdinal, typename Node::device_type> FCi;

  out << "version: " << MueLu::Version() << std::endl;

  Xpetra::UnderlyingLib lib          = MueLuTests::TestHelpers::Parameters::getLib();
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();
  //    GO gst_invalid = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
  //    GO lo_invalid = Teuchos::OrdinalTraits<LO>::invalid();
  int MyPID = comm->getRank();

  // Setup Levels
  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  // Build a pseudo-poisson test matrix
  FCi elem_to_node;
  RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC, LO, GO, NO>(num_p1_nodes, degree, elem_to_node, lib);
  fineLevel.Set("A", A);
  fineLevel.Set("pcoarsen: element to node map", rcp(&elem_to_node, false));

  // only one NS vector
  LocalOrdinal NSdim         = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->setSeed(846930886);
  nullSpace->randomize();
  fineLevel.Set("Nullspace", nullSpace);

  // ParameterList
  ParameterList Params;
  Params.set("pcoarsen: hi basis", hi_basis);
  Params.set("pcoarsen: lo basis", lo_basis);

  // Build P
  RCP<MueLu::IntrepidPCoarsenFactory<SC, LO, GO, NO>> IPCFact = rcp(new MueLu::IntrepidPCoarsenFactory<SC, LO, GO, NO>());
  IPCFact->SetParameterList(Params);
  coarseLevel.Request("P", IPCFact.get());  // request Ptent
  coarseLevel.Request("Nullspace", IPCFact.get());
  coarseLevel.Request("CoarseMap", IPCFact.get());
  coarseLevel.Request(*IPCFact);
  IPCFact->Build(fineLevel, coarseLevel);

  // Get P
  RCP<Matrix> P;
  coarseLevel.Get("P", P, IPCFact.get());
  RCP<CrsMatrix> Pcrs = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();
  // if(!MyPID) printf("P size = %d x %d\n",(int)P->getRangeMap()->getGlobalNumElements(),(int)P->getDomainMap()->getGlobalNumElements());
  out << "P size = "
      << (int)P->getRangeMap()->getGlobalNumElements()
      << " x "
      << (int)P->getDomainMap()->getGlobalNumElements()
      << std::endl;

  // Sanity
  if ((int)P->getRangeMap()->getGlobalNumElements() != (int)hi_gold_out.size())
    throw std::runtime_error("P range size does not match hi_gold_out");
  if ((int)P->getDomainMap()->getGlobalNumElements() != (int)lo_gold_in.size())
    throw std::runtime_error("P domain size does not match lo_gold_in");

  // Build serial comparison maps
  GO hi_num_global_dofs     = A->getRowMap()->getGlobalNumElements();
  GO hi_num_serial_elements = !MyPID ? hi_num_global_dofs : 0;
  RCP<Map> hi_SerialMap     = MapFactory::Build(lib, hi_num_global_dofs, hi_num_serial_elements, 0, comm);

  GO lo_num_global_dofs     = P->getDomainMap()->getGlobalNumElements();
  GO lo_num_serial_elements = !MyPID ? lo_num_global_dofs : 0;
  RCP<Map> lo_SerialMap     = MapFactory::Build(lib, lo_num_global_dofs, lo_num_serial_elements, 0, comm);

  RCP<Export> lo_importer = ExportFactory::Build(lo_SerialMap, P->getDomainMap());
  RCP<Export> hi_importer = ExportFactory::Build(A->getRowMap(), hi_SerialMap);

  // Allocate some vectors
  RCP<Vector> s_InVec      = VectorFactory::Build(lo_SerialMap);
  RCP<Vector> p_InVec      = VectorFactory::Build(P->getDomainMap());
  RCP<Vector> s_OutVec     = VectorFactory::Build(hi_SerialMap);
  RCP<Vector> s_codeOutput = VectorFactory::Build(hi_SerialMap);
  RCP<Vector> p_codeOutput = VectorFactory::Build(A->getRowMap());

  // Fill serial GOLD vecs on Proc 0
  if (!MyPID) {
    for (size_t i = 0; i < (size_t)lo_gold_in.size(); i++)
      s_InVec->replaceLocalValue(i, lo_gold_in[i]);

    for (size_t i = 0; i < (size_t)hi_gold_out.size(); i++)
      s_OutVec->replaceLocalValue(i, hi_gold_out[i]);
  }

  // Migrate input data
  p_InVec->doExport(*s_InVec, *lo_importer, Xpetra::ADD);

  // Apply P
  P->apply(*p_InVec, *p_codeOutput);

  // Migrate Output data
  s_codeOutput->doExport(*p_codeOutput, *hi_importer, Xpetra::ADD);

  // Compare vs. GOLD
  s_codeOutput->update(-1.0, *s_OutVec, 1.0);
  Teuchos::Array<MT> norm2(1);
  s_codeOutput->norm2(norm2());

  // if(!MyPID) printf("Diff norm = %10.4e\n",norm2[0]);
  out << "Diff norm = " << norm2[0] << std::endl;
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, BuildP_PseudoPoisson_p2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode);
  // GOLD vector collection
  std::vector<Scalar> p2_gold_in  = {/*0,*/ 1, 2, 3, 4, 5, 6, 7, 8, /*9*/};  // Ignore Dirichlet unknowns
  std::vector<Scalar> p2_gold_out = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                     0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5};
  TestPseudoPoisson<Scalar, LocalOrdinal, GlobalOrdinal, Node>(out, 2 + p2_gold_in.size(), 2, p2_gold_in, p2_gold_out, std::string("hgrad_line_c2"));
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, BuildP_PseudoPoisson_p3, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode);
  // GOLD vector collection
  size_t total_num_points = 10;
  int degree              = 3;
  std::vector<Scalar> gold_in(total_num_points - 2);
  std::vector<Scalar> gold_out(total_num_points + (total_num_points - 1) * (degree - 1));
  for (size_t i = 0; i < total_num_points - 2; i++) {
    gold_in[i]      = i + 1;
    gold_out[i + 1] = i + 1;
  }
  gold_out[0]                    = 0.0;
  gold_out[total_num_points - 1] = total_num_points - 1;

  size_t idx = total_num_points;
  for (size_t i = 0; i < total_num_points - 1; i++) {
    for (size_t j = 0; j < (size_t)degree - 1; j++) {
      gold_out[idx] = i + ((double)j + 1) / degree;
      idx++;
    }
  }

  TestPseudoPoisson<Scalar, LocalOrdinal, GlobalOrdinal, Node>(out, 2 + gold_in.size(), 3, gold_in, gold_out, std::string("hgrad_line_c3"));
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, BuildP_PseudoPoisson_p4, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode);

  // GOLD vector collection
  size_t total_num_points = 10;
  int degree              = 4;

  std::vector<Scalar> gold_in(total_num_points - 2);
  std::vector<Scalar> gold_out(total_num_points + (total_num_points - 1) * (degree - 1));
  for (size_t i = 0; i < total_num_points - 2; i++) {
    gold_in[i]      = i + 1;
    gold_out[i + 1] = i + 1;
  }
  gold_out[0]                    = 0.0;
  gold_out[total_num_points - 1] = total_num_points - 1;

  size_t idx = total_num_points;
  for (size_t i = 0; i < total_num_points - 1; i++) {
    for (size_t j = 0; j < (size_t)degree - 1; j++) {
      gold_out[idx] = i + ((double)j + 1) / degree;
      idx++;
    }
  }

  TestPseudoPoisson<Scalar, LocalOrdinal, GlobalOrdinal, Node>(out, 2 + gold_in.size(), 4, gold_in, gold_out, std::string("hgrad_line_c4"));
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, CreatePreconditioner_p2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode);
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
#if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#endif
#if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#endif

  typedef Scalar SC;
  typedef GlobalOrdinal GO;
  typedef LocalOrdinal LO;
  typedef Node NO;
  typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef Kokkos::DynRankView<LocalOrdinal, typename Node::device_type> FCi;

  out << "version: " << MueLu::Version() << std::endl;
  using Teuchos::RCP;
  int degree = 2;
  std::string hi_basis("hgrad_line_c2");

  Xpetra::UnderlyingLib lib          = TestHelpers::Parameters::getLib();
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  GO num_nodes = 972;
  // Build a pseudo-poisson test matrix
  FCi elem_to_node;
  RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC, LO, GO, NO>(num_nodes, degree, elem_to_node, lib);

  // Normalized RHS
  RCP<MultiVector> RHS1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  RHS1->setSeed(846930886);
  RHS1->randomize();
  Teuchos::Array<MT> norms(1);
  RHS1->norm2(norms);
  RHS1->scale(1 / norms[0]);

  // Zero initial guess
  RCP<MultiVector> X1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  X1->putScalar(Teuchos::ScalarTraits<SC>::zero());

  // ParameterList
  ParameterList Params, level0;
  Params.set("multigrid algorithm", "pcoarsen");
  //    Params.set("rap: fix zero diagonals",true);
  Params.set("pcoarsen: hi basis", hi_basis);
  Params.set("pcoarsen: lo basis", "hgrad_line_c1");
  Params.set("verbosity", "high");
  Params.set("max levels", 2);
  //     if(lib==Xpetra::UseEpetra) Params.set("coarse: type","RELAXATION");// FIXME remove when we sort out the OAZ issue
  Params.set("coarse: max size", 100);
  level0.set("pcoarsen: element to node map", rcp(&elem_to_node, false));
  Params.set("level 0", level0);
  Params.set("verbosity", "none");

  // Build hierarchy
  RCP<Hierarchy> tH = MueLu::CreateXpetraPreconditioner<SC, LO, GO, NO>(A, Params);
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, CreatePreconditioner_p3, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode);
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
#if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#endif
#if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#endif

  typedef Scalar SC;
  typedef GlobalOrdinal GO;
  typedef LocalOrdinal LO;
  typedef Node NO;
  typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef Kokkos::DynRankView<LocalOrdinal, typename Node::device_type> FCi;

  out << "version: " << MueLu::Version() << std::endl;
  using Teuchos::RCP;
  int degree = 3;
  std::string hi_basis("hgrad_line_c3");

  Xpetra::UnderlyingLib lib          = TestHelpers::Parameters::getLib();
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  GO num_nodes = 972;
  // Build a pseudo-poisson test matrix
  FCi elem_to_node;
  RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC, LO, GO, NO>(num_nodes, degree, elem_to_node, lib);

  // Normalized RHS
  RCP<MultiVector> RHS1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  RHS1->setSeed(846930886);
  RHS1->randomize();
  Teuchos::Array<MT> norms(1);
  RHS1->norm2(norms);
  RHS1->scale(1 / norms[0]);

  // Zero initial guess
  RCP<MultiVector> X1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  X1->putScalar(Teuchos::ScalarTraits<SC>::zero());

  // ParameterList
  ParameterList Params, level0;
  Params.set("multigrid algorithm", "pcoarsen");
  //    Params.set("rap: fix zero diagonals",true);
  Params.set("pcoarsen: hi basis", hi_basis);
  Params.set("pcoarsen: lo basis", "hgrad_line_c1");
  Params.set("verbosity", "high");
  Params.set("max levels", 2);
  Params.set("coarse: max size", 100);
  level0.set("pcoarsen: element to node map", rcp(&elem_to_node, false));
  Params.set("level 0", level0);
  Params.set("verbosity", "none");

  // Build hierarchy
  RCP<Hierarchy> tH = MueLu::CreateXpetraPreconditioner<SC, LO, GO, NO>(A, Params);
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, CreatePreconditioner_p4, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode);
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
#if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#endif
#if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#endif

  typedef Scalar SC;
  typedef GlobalOrdinal GO;
  typedef LocalOrdinal LO;
  typedef Node NO;
  typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef Kokkos::DynRankView<LocalOrdinal, typename Node::device_type> FCi;

  out << "version: " << MueLu::Version() << std::endl;
  using Teuchos::RCP;
  int degree = 4;
  std::string hi_basis("hgrad_line_c4");

  Xpetra::UnderlyingLib lib          = TestHelpers::Parameters::getLib();
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  GO num_nodes = 972;
  // Build a pseudo-poisson test matrix
  FCi elem_to_node;
  RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC, LO, GO, NO>(num_nodes, degree, elem_to_node, lib);

  // Normalized RHS
  RCP<MultiVector> RHS1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  RHS1->setSeed(846930886);
  RHS1->randomize();
  Teuchos::Array<MT> norms(1);
  RHS1->norm2(norms);
  RHS1->scale(1 / norms[0]);

  // Zero initial guess
  RCP<MultiVector> X1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  X1->putScalar(Teuchos::ScalarTraits<SC>::zero());

  // ParameterList
  ParameterList Params, level0;
  Params.set("multigrid algorithm", "pcoarsen");
  //    Params.set("rap: fix zero diagonals",true);
  Params.set("pcoarsen: hi basis", hi_basis);
  Params.set("pcoarsen: lo basis", "hgrad_line_c1");
  Params.set("verbosity", "high");
  Params.set("max levels", 2);
  Params.set("coarse: max size", 100);
  level0.set("pcoarsen: element to node map", rcp(&elem_to_node, false));
  Params.set("level 0", level0);
  Params.set("verbosity", "none");

  // Build hierarchy
  RCP<Hierarchy> tH = MueLu::CreateXpetraPreconditioner<SC, LO, GO, NO>(A, Params);
}

/*********************************************************************************************************************/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class Basis>
bool test_representative_basis(Teuchos::FancyOStream &out, const std::string &name, Intrepid2::EPointType ptype, int max_degree) {
#include "MueLu_UseShortNames.hpp"
  typedef Scalar SC;
  typedef GlobalOrdinal GO;
  typedef LocalOrdinal LO;
  typedef Node NO;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef Intrepid2::CellTools<typename Node::device_type::execution_space> CellTools;
  typedef Kokkos::DynRankView<MT, typename Node::device_type> FC;

  out << "version: " << MueLu::Version() << std::endl;

  bool success               = true;
  int combinationTestedCount = 0;

  // Contruct a container that has the reference coordinates for the domain cell topology
  RCP<Basis> linearBasis        = rcp(new Basis(1, ptype));
  shards::CellTopology cellTopo = linearBasis->getBaseCellTopology();
  int spaceDim                  = cellTopo.getDimension();
  int vertexCount               = cellTopo.getVertexCount();
  FC refCellVertices, refCellVertex, physCellVertices, physCellVerticesPermuted;
  /*
   How many physical cells do we need to test adequately?  I believe the answer is just 2, even in 3D;
   by running through all the symmetries we should have all the edge-to-edge pairings covered.
   */
  int numCells                = 2;
  double xTranslationForCell1 = 2.0;  // shift to the right by 2
  double pointTol             = 1e-12;

  resize(physCellVertices, numCells, vertexCount, spaceDim);
  resize(physCellVerticesPermuted, numCells, vertexCount, spaceDim);
  resize(refCellVertices, vertexCount, spaceDim);
  resize(refCellVertex, 3);

  // regardless of spatial extent, CellTools::getReferenceVertex() populates 3 slots
  for (int vertexOrdinal = 0; vertexOrdinal < vertexCount; vertexOrdinal++) {
    CellTools::getReferenceVertex(refCellVertex, cellTopo, vertexOrdinal);
    // UVM used here, accessing vertex coordinates on host that were populated on device.

    Kokkos::parallel_for(
        Kokkos::RangePolicy<typename Node::device_type::execution_space>(0, spaceDim), KOKKOS_LAMBDA(int d) {
          refCellVertices(vertexOrdinal, d) = refCellVertex(d);
          //      cout << "refCellVertices(" << vertexOrdinal << "," << d << ") = " << refCellVertex(d) << endl;
          // cell 0 is just the reference cell:
          physCellVertices(0, vertexOrdinal, d) = refCellVertex(d);
          // cell 1 is the reference cell, except that the x coords get translated
          // NOTE: this will need to change to support non-hypercube topologies
          physCellVertices(1, vertexOrdinal, d) = refCellVertex(d) + ((d == 0) ? xTranslationForCell1 : 0);
        });
    Kokkos::fence();
  }
  Symmetries cellSymmetries = Symmetries::shardsSymmetries(cellTopo);
  int symmetryCount         = cellSymmetries.getPermutationCount();

  for (int highPolyDegree = 1; highPolyDegree < max_degree; highPolyDegree++) {
    FC hi_DofCoords, lo_DofCoords;
    RCP<Basis> hi = rcp(new Basis(highPolyDegree, ptype));
    Kokkos::resize(hi_DofCoords, hi->getCardinality(), hi->getBaseCellTopology().getDimension());
    hi->getDofCoords(hi_DofCoords);
    auto hi_DofCoords_host = Kokkos::create_mirror_view(hi_DofCoords);
    Kokkos::deep_copy(hi_DofCoords_host, hi_DofCoords);

    // we'll want to create a global numbering for both high and low order bases
    // --> we make a lambda function that accepts FC with dof coords as argument
    auto getTwoCellNumbering = [pointTol, numCells, spaceDim, xTranslationForCell1](const typename FC::HostMirror &dofCoords) -> UniqueNumbering {
      int dofsPerCell = dofCoords.extent(0);

      vector<double> coords(spaceDim);
      UniqueNumbering numbering(spaceDim, pointTol);

      // number the guys on cell 0 first, then the ones on cell 1:
      for (int cellOrdinal = 0; cellOrdinal < numCells; cellOrdinal++) {
        for (int dofOrdinal = 0; dofOrdinal < dofsPerCell; dofOrdinal++) {
          for (int d = 0; d < spaceDim; d++) {
            if ((d == 0) && (cellOrdinal == 1)) {
              coords[d] = dofCoords(dofOrdinal, d) + xTranslationForCell1;
            } else {
              coords[d] = dofCoords(dofOrdinal, d);
            }
          }
          numbering.getGlobalID(coords);
        }
      }
      return numbering;
    };

    UniqueNumbering hiNumbering = getTwoCellNumbering(hi_DofCoords_host);
    out << "Total dof count two cells of degree " << highPolyDegree << ": ";
    out << hiNumbering.totalCount() << endl;

    for (int lowPolyDegree = 1; lowPolyDegree < highPolyDegree; lowPolyDegree++) {
      out << "Testing with high order " << highPolyDegree << ", low order " << lowPolyDegree << endl;
      RCP<Basis> lo = rcp(new Basis(lowPolyDegree, ptype));
      Kokkos::resize(lo_DofCoords, lo->getCardinality(), lo->getBaseCellTopology().getDimension());
      lo->getDofCoords(lo_DofCoords);
      auto lo_DofCoords_host = Kokkos::create_mirror_view(lo_DofCoords);
      Kokkos::deep_copy(lo_DofCoords_host, lo_DofCoords);
      UniqueNumbering loNumbering = getTwoCellNumbering(lo_DofCoords_host);

      // print out the high/low global numbering along the x=1 interface:
      out << "Low-order global IDs along intercell interface:\n";
      for (int lowOrdinal = 0; lowOrdinal < lo->getCardinality(); lowOrdinal++) {
        vector<double> coords(lo_DofCoords.extent(1));
        for (int d = 0; d < int(lo_DofCoords.extent(1)); d++) {
          coords[d] = lo_DofCoords_host(lowOrdinal, d);
        }
        if (coords[0] == 1.0) {
          int globalOrdinal = loNumbering.getGlobalID(coords);
          out << globalOrdinal << ": (";
          for (int d = 0; d < int(coords.size()) - 1; d++) {
            out << coords[d] << ",";
          }
          out << coords[coords.size() - 1] << ")\n";
        }
      }
      // print out the high/low global numbering along the x=1 interface:
      out << "High-order global IDs along intercell interface:\n";
      for (int highOrdinal = 0; highOrdinal < hi->getCardinality(); highOrdinal++) {
        vector<double> coords(hi_DofCoords.extent(1));
        for (int d = 0; d < int(hi_DofCoords.extent(1)); d++) {
          coords[d] = hi_DofCoords_host(highOrdinal, d);
        }
        if (coords[0] == 1.0) {
          int globalOrdinal = hiNumbering.getGlobalID(coords);
          out << globalOrdinal << ": (";
          for (int d = 0; d < int(coords.size()) - 1; d++) {
            out << coords[d] << ",";
          }
          out << coords[coords.size() - 1] << ")\n";
        }
      }

      // Get the candidates
      double threshold = 1e-10;
      std::vector<std::vector<size_t>> candidates;
      MueLu::MueLuIntrepid::GenerateRepresentativeBasisNodes<Basis, FC>(*lo, hi_DofCoords, threshold, candidates);

      // Correctness Test 1: Make sure that there are no duplicates in the representative lists / no low DOF has no candidates
      std::vector<bool> is_candidate(hi_DofCoords.extent(0), false);
      bool no_doubles = true;
      for (int lowOrderDof = 0; no_doubles && lowOrderDof < (int)candidates.size(); lowOrderDof++) {
        if (candidates[lowOrderDof].size() == 0) no_doubles = false;  // this low DOF has no candidates!
        for (int l = 0; l < (int)candidates[lowOrderDof].size(); l++)
          if (is_candidate[candidates[lowOrderDof][l]] == false)
            is_candidate[candidates[lowOrderDof][l]] = true;
          else {
            no_doubles = false;
            break;
          }  // this high-order dof was already claimed by an earlier low DOF!
      }
      if (!no_doubles) {
        out << "ERROR: " << name << " The 'no duplicates' test fails w/ lo/hi = " << lowPolyDegree << "/" << highPolyDegree << std::endl;
        return false;
      }

      // Correctness Test 2: Try 2 elements, in all possible relative orientations, and confirm that the
      //                     "lowest global ordinal" tie-breaker always returns the same thing for both neighbors
      auto physCellVertices_host         = Kokkos::create_mirror_view(physCellVertices);
      auto physCellVerticesPermuted_host = Kokkos::create_mirror_view(physCellVerticesPermuted);
      Kokkos::deep_copy(physCellVertices_host, physCellVertices);
      for (int permOrdinal0 = 0; permOrdinal0 < symmetryCount; permOrdinal0++) {
        vector<int> perm0 = cellSymmetries.getPermutation(permOrdinal0);
        for (int vertexOrdinal = 0; vertexOrdinal < vertexCount; vertexOrdinal++) {
          int mappedVertexOrdinal = perm0[vertexOrdinal];
          for (int d = 0; d < spaceDim; d++) {
            physCellVerticesPermuted_host(0, vertexOrdinal, d) = physCellVertices_host(0, mappedVertexOrdinal, d);
          }
        }

        // brute force search for the side of cell shared with neighbor
        // this is the one that has points with x coordinates equal to 1.0
        // we'll want to do this once for cell 0, and once for cell 1, so we make it a lambda
        // (NOTE: this will need to change for non-hypercube topology support)
        auto searchForX1Side = [cellTopo, physCellVerticesPermuted_host, spaceDim](int cellOrdinal) -> int {
          // Line<2> gives wrong answers for getSideCount() and getNodeCount(), so we handle 1D case separately:
          if (spaceDim == 1) {
            int sideCount = 2;
            for (int sideVertexOrdinal = 0; sideVertexOrdinal < sideCount; sideVertexOrdinal++) {
              int cellVertexOrdinal = sideVertexOrdinal;
              if (physCellVerticesPermuted_host(cellOrdinal, cellVertexOrdinal, 0) == 1.0) {
                return sideVertexOrdinal;
              }
            }
            return -1;
          }
          int sideCount = (int)cellTopo.getSideCount();
          for (int sideOrdinal = 0; sideOrdinal < sideCount; sideOrdinal++) {
            int sideVertexCount = cellTopo.getNodeCount(spaceDim - 1, sideOrdinal);
            bool matchFound     = true;
            for (int sideVertexOrdinal = 0; sideVertexOrdinal < sideVertexCount; sideVertexOrdinal++) {
              int cellVertexOrdinal = cellTopo.getNodeMap(spaceDim - 1, sideOrdinal, sideVertexOrdinal);
              if (physCellVerticesPermuted_host(cellOrdinal, cellVertexOrdinal, 0) != 1.0) {
                matchFound = false;
                break;
              }
            }
            if (matchFound) {
              return sideOrdinal;
            }
          }
          return -1;
        };

        int cell0Side = searchForX1Side(0);
        //        out << "cell 0 side is " << cell0Side << endl;

        for (int permOrdinal1 = 0; permOrdinal1 < symmetryCount; permOrdinal1++) {
          vector<int> perm1 = cellSymmetries.getPermutation(permOrdinal1);
          for (int vertexOrdinal = 0; vertexOrdinal < vertexCount; vertexOrdinal++) {
            int mappedVertexOrdinal = perm1[vertexOrdinal];
            for (int d = 0; d < spaceDim; d++) {
              physCellVerticesPermuted_host(1, vertexOrdinal, d) = physCellVertices_host(1, mappedVertexOrdinal, d);
            }
          }
          // get the mapped dof coords for lo and high bases:
          FC lo_physDofCoords, hi_physDofCoords;
          Kokkos::resize(lo_physDofCoords, numCells, lo->getCardinality(), cellTopo.getDimension());
          Kokkos::resize(hi_physDofCoords, numCells, hi->getCardinality(), cellTopo.getDimension());

          Kokkos::deep_copy(physCellVerticesPermuted, physCellVerticesPermuted_host);
          CellTools::mapToPhysicalFrame(lo_physDofCoords, lo_DofCoords, physCellVerticesPermuted, cellTopo);
          CellTools::mapToPhysicalFrame(hi_physDofCoords, hi_DofCoords, physCellVerticesPermuted, cellTopo);
          auto lo_physDofCoords_host = Kokkos::create_mirror_view(lo_physDofCoords);
          auto hi_physDofCoords_host = Kokkos::create_mirror_view(hi_physDofCoords);
          Kokkos::deep_copy(lo_physDofCoords_host, lo_physDofCoords);
          Kokkos::deep_copy(hi_physDofCoords_host, hi_physDofCoords);
          Kokkos::fence();  // mapToPhysicalFrame calls getValues which calls kernels, so fence is required before UVM reads below

          int cell1Side = searchForX1Side(1);

          /*
           We want to verify that the neighboring cells (as permuted) agree on their
           representation of the low-order basis in their selection of high-order basis
           members.  It suffices to check that the *set* of candidate high-order basis functions
           is the same; then for any global numbering the "lowest global number" tie-breaker
           will select the same function.

           The logic is this.  We have:
           - a global numbering for the physical dof coords of the coarse discretization
           - a global numbering for the physical dof coords of the fine discretization
           - a local map from lo to hi on the reference cell
           - a map from local dofs on each cell to the physical dof coords for the {lo|hi} basis

           On each cell, then, we can map from the local low dof ordinal to the local high dof ordinal.
           At the same time, we can map the local ordinals for high and low to global ordinals, giving
           us a mapping from global low ordinals to global high ordinals.  We then can compare the mapping
           to verify that the two cells agree.
           */

          auto constructMap = [lo, lo_physDofCoords_host, &loNumbering, candidates, hi_physDofCoords_host, &hiNumbering, spaceDim](int cellOrdinal, int sideOrdinal) -> map<int, set<int>> {
            map<int, set<int>> globalLowToHighMap;

            vector<int> cell_loDofOrdinals = localDofOrdinalsForSide(lo, sideOrdinal);
            for (int lowLocalOrdinal : cell_loDofOrdinals) {
              vector<double> loCoords(spaceDim);
              for (int d = 0; d < spaceDim; d++) {
                loCoords[d] = lo_physDofCoords_host(cellOrdinal, lowLocalOrdinal, d);
              }
              int lowGlobalNumber              = loNumbering.getGlobalID(loCoords);
              vector<size_t> highLocalOrdinals = candidates[lowLocalOrdinal];
              for (int highLocalOrdinal : highLocalOrdinals) {
                vector<double> hiCoords(spaceDim);
                for (int d = 0; d < spaceDim; d++) {
                  hiCoords[d] = hi_physDofCoords_host(cellOrdinal, highLocalOrdinal, d);
                }
                int highGlobalNumber = hiNumbering.getGlobalID(hiCoords);
                globalLowToHighMap[lowGlobalNumber].insert(highGlobalNumber);
              }
            }
            return globalLowToHighMap;
          };

          map<int, set<int>> cell0_mapping = constructMap(0, cell0Side);
          map<int, set<int>> cell1_mapping = constructMap(1, cell1Side);

          // verify the two maps are identical:
          if (cell0_mapping.size() != cell1_mapping.size()) {
            success = false;
            out << "cell 0 and cell 1 mapping do not have the same number of low global ordinal entries!\n";
          } else {
            auto cell0_it = cell0_mapping.begin();
            auto cell1_it = cell1_mapping.begin();
            while (cell0_it != cell0_mapping.end()) {
              if (cell0_it->first != cell1_it->first) {
                success = false;
                out << "cell 0 and cell 1 have different low global ordinal entries; ";
                out << cell0_it->first << " != " << cell1_it->first << endl;
              } else if (cell0_it->second.size() != cell1_it->second.size()) {
                success = false;
                out << "cell 0 and cell 1 have a different number of high global ordinal entries for low global ordinal ";
                out << cell0_it->first << cell0_it->second.size() << " != " << cell1_it->second.size() << endl;
              } else {
                // check that the set contents are identical
                auto cell0_set_it = cell0_it->second.begin();
                auto cell1_set_it = cell1_it->second.begin();
                bool setsMatch    = true;
                while (cell0_set_it != cell0_it->second.end()) {
                  if (*cell0_set_it != *cell1_set_it) {
                    setsMatch = false;
                    break;
                  }
                  cell0_set_it++;
                  cell1_set_it++;
                }
                if (!setsMatch) {
                  int lowGlobalOrdinal = cell0_it->first;
                  success              = false;
                  out << "with permutation selections " << permOrdinal0 << " and " << permOrdinal1 << ", ";
                  out << "cell 0 and cell 1 maps differ in the high global entries for low global ordinal ";
                  out << lowGlobalOrdinal << ": { ";
                  for (int highGlobalOrdinal : cell0_it->second) {
                    out << highGlobalOrdinal << " ";
                  }
                  out << "} != { ";
                  for (int highGlobalOrdinal : cell1_it->second) {
                    out << highGlobalOrdinal << " ";
                  }
                  out << "}\n";
                }
              }

              cell0_it++;
              cell1_it++;
            }
          }

          combinationTestedCount++;
        }
      }
    }
  }
  cout << "Tested " << combinationTestedCount << " lo/hi, two-cell permutation combinations.\n";
  return success;
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_LINE_Equispaced, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef typename Node::device_type::execution_space ES;
  typedef Kokkos::DynRankView<MT, typename Node::device_type> FC;
  typedef Intrepid2::Basis_HGRAD_LINE_Cn_FEM<ES, MT, MT> Basis;

  bool rv = test_representative_basis<Scalar, LocalOrdinal, GlobalOrdinal, Node, Basis>(out, " GenerateRepresentativeBasisNodes_LINE_EQUISPACED", Intrepid2::POINTTYPE_EQUISPACED, MAX_LINE_DEGREE);
  TEST_EQUALITY(rv, true);
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_QUAD_Equispaced, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef typename Node::device_type::execution_space ES;
  typedef Kokkos::DynRankView<MT, typename Node::device_type> FC;
  typedef Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES, MT, MT> Basis;

  bool rv = test_representative_basis<Scalar, LocalOrdinal, GlobalOrdinal, Node, Basis>(out, " GenerateRepresentativeBasisNodes_QUAD_EQUISPACED", Intrepid2::POINTTYPE_EQUISPACED, MAX_QUAD_DEGREE);
  TEST_EQUALITY(rv, true);
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_QUAD_Spectral, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef typename Node::device_type::execution_space ES;
  typedef Kokkos::DynRankView<MT, typename Node::device_type> FC;
  typedef Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES, MT, MT> Basis;

  const Intrepid2::EPointType POINTTYPE_SPECTRAL = static_cast<Intrepid2::EPointType>(1);  // Not sure why I have to do this...
  bool rv                                        = test_representative_basis<Scalar, LocalOrdinal, GlobalOrdinal, Node, Basis>(out, " GenerateRepresentativeBasisNodes_QUAD_SPECTRAL", POINTTYPE_SPECTRAL, MAX_QUAD_DEGREE);
  TEST_EQUALITY(rv, true);
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_HEX_Equispaced, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef typename Node::device_type::execution_space ES;
  typedef Kokkos::DynRankView<MT, typename Node::device_type> FC;
  typedef Intrepid2::Basis_HGRAD_HEX_Cn_FEM<ES, MT, MT> Basis;

  bool rv = test_representative_basis<Scalar, LocalOrdinal, GlobalOrdinal, Node, Basis>(out, " GenerateRepresentativeBasisNodes_HEX_EQUISPACED", Intrepid2::POINTTYPE_EQUISPACED, MAX_HEX_DEGREE);
  TEST_EQUALITY(rv, true);
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_HEX_Spectral, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef typename Node::device_type::execution_space ES;
  typedef Kokkos::DynRankView<MT, typename Node::device_type> FC;
  typedef Intrepid2::Basis_HGRAD_HEX_Cn_FEM<ES, MT, MT> Basis;

  const Intrepid2::EPointType POINTTYPE_SPECTRAL = static_cast<Intrepid2::EPointType>(1);  // Not sure why I have to do this...
  bool rv                                        = test_representative_basis<Scalar, LocalOrdinal, GlobalOrdinal, Node, Basis>(out, " GenerateRepresentativeBasisNodes_HEX_SPECTRAL", POINTTYPE_SPECTRAL, MAX_HEX_DEGREE);
  TEST_EQUALITY(rv, true);
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, GenerateLoNodeInHighViaGIDs_QUAD_pn_to_p1, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // NOTE: We need more tests for this that do pn to pm
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef GlobalOrdinal GO;
  typedef LocalOrdinal LO;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;

  typedef typename Node::device_type::execution_space ES;
  typedef Kokkos::DynRankView<MT, typename Node::device_type> FC;
  typedef Kokkos::DynRankView<LO, typename Node::device_type> FCi;
  typedef Intrepid2::Basis_HGRAD_QUAD_C1_FEM<ES, MT, MT> LoBasis;
  typedef Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES, MT, MT> HiBasis;

  out << "version: " << MueLu::Version() << std::endl;
  Xpetra::UnderlyingLib lib          = MueLuTests::TestHelpers::Parameters::getLib();
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  int max_degree   = MAX_QUAD_DEGREE;
  double threshold = 1e-10;
  GO gst_invalid   = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();

  // Lo
  RCP<LoBasis> lo = rcp(new LoBasis());
  size_t numLo    = lo->getCardinality();

  for (int i = 0; i < max_degree; i++) {
    RCP<HiBasis> hi = rcp(new HiBasis(i, Intrepid2::POINTTYPE_EQUISPACED));
    size_t numHi    = hi->getCardinality();

    // The quad-only stuff
    std::vector<size_t> lo_node_in_hi;
    FC hi_dofCoords;
    MueLu::MueLuIntrepid::IntrepidGetP1NodeInHi<MT, typename Node::device_type>(hi, lo_node_in_hi, hi_dofCoords);
    FCi hi_e2n("hi_e2n", 1, numHi);
    FCi lo_e2n("lo_e2n", 1, numLo);

    // Dummy elem2node map
    Kokkos::parallel_for(
        "IntrepidPCoarsenFactory,GenerateLoNodeInHighViaGIDs_QUAD_pn_to_p1",
        Kokkos::RangePolicy<typename Node::device_type::execution_space>(0, numHi), KOKKOS_LAMBDA(int j) {
          hi_e2n(0, j) = j;
        });
    Kokkos::fence();
    Teuchos::Array<GO> hi_colids(numHi);
    for (size_t j = 0; j < numHi; j++) {
      hi_colids[j] = j;
    }

    // Dummy column map
    RCP<const Map> hi_colMap = MapFactory::Build(lib, gst_invalid, hi_colids(), 0, comm);

    // The dynamic stuff
    std::vector<std::vector<size_t>> candidates;
    MueLu::MueLuIntrepid::GenerateRepresentativeBasisNodes<LoBasis, FC>(*lo, hi_dofCoords, threshold, candidates);
    MueLu::MueLuIntrepid::GenerateLoNodeInHiViaGIDs<LO, GO, Node, FCi>(candidates, hi_e2n, hi_colMap, lo_e2n);

    // Compare and make sure we're cool
    bool node_diff   = false;
    auto lo_e2n_host = Kokkos::create_mirror_view(lo_e2n);
    Kokkos::deep_copy(lo_e2n_host, lo_e2n);
    for (size_t j = 0; j < numLo; j++)
      if (lo_node_in_hi[j] != (size_t)lo_e2n_host(0, j)) node_diff = true;
#if 0
      printf("[%d] Comparison = ",i);
      for(size_t j=0; j<numLo; j++)
        printf("%d|%d ",(int)lo_node_in_hi[j],(int)lo_e2n(0,j));
      printf("\n");
#endif
    TEST_EQUALITY(node_diff, false);
  }
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, BuildLoElemToNodeViaRepresentatives_QUAD_pn_to_p1, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // NOTE: We need more tests for this that do pn to pm
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType MT;
  typedef Kokkos::DynRankView<MT, typename Node::device_type> FC;
  typedef Kokkos::DynRankView<LocalOrdinal, typename Node::device_type> FCi;
  typedef typename Node::device_type::execution_space ES;
  typedef Intrepid2::Basis<ES, MT, MT> Basis;
  GO gst_invalid = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();

  out << "version: " << MueLu::Version() << std::endl;
  Xpetra::UnderlyingLib lib          = MueLuTests::TestHelpers::Parameters::getLib();
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  int max_degree = Intrepid2::Parameters::MaxOrder;

  {
    // QUAD
    //  A one element test with Kirby-numbered nodes where the top edge is not owned
    RCP<Basis> lo = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<ES, MT, MT>());

    for (int degree = 2; degree < max_degree; degree++) {
      int Nn        = (degree + 1) * (degree + 1);
      RCP<Basis> hi = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES, MT, MT>(degree, Intrepid2::POINTTYPE_EQUISPACED));
      FCi hi_e2n("hi_e2n", 1, Nn), lo_e2n, lo_e2n_mk2;

      std::vector<bool> hi_owned(Nn, false), lo_owned, lo_owned_mk2;
      std::vector<size_t> lo_node_in_hi;
      std::vector<LO> hi_to_lo_map, hi_to_lo_map_mk2;
      int lo_numOwnedNodes = 0, lo_numOwnedNodes_mk2 = 0;
      FC hi_dofCoords;

      // El2node / ownership / colmap
      Kokkos::parallel_for(
          "IntrepidPCoarsenFactory,BuildLoElemToNodeViaRepresentatives_QUAD_pn_to_p1",
          Kokkos::RangePolicy<typename Node::device_type::execution_space>(0, Nn), KOKKOS_LAMBDA(int i) {
            hi_e2n(0, i) = i;
          });
      Kokkos::fence();
      Teuchos::Array<GO> hi_colids(Nn);
      for (int i = 0; i < Nn; i++) {
        hi_colids[i] = i;
        if (i < Nn - (degree + 1)) hi_owned[i] = true;
      }

      Teuchos::ArrayRCP<const int> is_dirichlet(Nn, 0);

      /*** Do stuff the injection way ***/
      MueLu::MueLuIntrepid::IntrepidGetP1NodeInHi<MT, typename Node::device_type>(hi, lo_node_in_hi, hi_dofCoords);
      MueLu::MueLuIntrepid::BuildLoElemToNode(hi_e2n, hi_owned, lo_node_in_hi, is_dirichlet, lo_e2n, lo_owned, hi_to_lo_map, lo_numOwnedNodes);

      /*** Do stuff the representative way ***/
      RCP<const Map> hi_colMap = MapFactory::Build(lib, gst_invalid, hi_colids(), 0, comm);
      FCi lo_elemToHiRepresentativeNode;
      double threshold = 1e-10;
      std::vector<std::vector<size_t>> candidates;
      MueLu::MueLuIntrepid::GenerateRepresentativeBasisNodes<Basis, FC>(*lo, hi_dofCoords, threshold, candidates);
      MueLu::MueLuIntrepid::GenerateLoNodeInHiViaGIDs(candidates, hi_e2n, hi_colMap, lo_elemToHiRepresentativeNode);
      MueLu::MueLuIntrepid::BuildLoElemToNodeViaRepresentatives(hi_e2n, hi_owned, lo_elemToHiRepresentativeNode, lo_e2n_mk2, lo_owned_mk2, hi_to_lo_map_mk2, lo_numOwnedNodes_mk2);

      // Compare stuff
      TEST_EQUALITY(lo_numOwnedNodes, 2);
      TEST_EQUALITY(lo_numOwnedNodes_mk2, 2);

      size_t num_lo_nodes_located = 0;
      for (size_t i = 0; i < hi_to_lo_map.size(); i++) {
        if (hi_to_lo_map[i] != Teuchos::OrdinalTraits<LO>::invalid())
          num_lo_nodes_located++;
      }
      TEST_EQUALITY(lo_owned.size(), num_lo_nodes_located);
      TEST_EQUALITY(lo_owned_mk2.size(), num_lo_nodes_located);

      auto lo_e2n_host     = Kokkos::create_mirror_view(lo_e2n);
      auto lo_e2n_mk2_host = Kokkos::create_mirror_view(lo_e2n_mk2);
      Kokkos::deep_copy(lo_e2n_host, lo_e2n);
      Kokkos::deep_copy(lo_e2n_mk2_host, lo_e2n_mk2);
      for (size_t i = 0; i < lo_e2n.extent(0); i++)
        for (size_t j = 0; j < lo_e2n.extent(1); j++)
          TEST_EQUALITY(lo_e2n_host(i, j), lo_e2n_mk2_host(i, j));

      for (size_t i = 0; i < (size_t)lo_owned.size(); i++)
        TEST_EQUALITY(lo_owned[i], lo_owned_mk2[i]);

      TEST_EQUALITY(hi_to_lo_map.size(), hi_to_lo_map_mk2.size());
      for (size_t i = 0; i < (size_t)hi_to_lo_map.size(); i++)
        TEST_EQUALITY(hi_to_lo_map[i], hi_to_lo_map_mk2[i]);
    }
  }  // end QUAD
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, BuildP_PseudoPoisson_LINE_p3_to_p2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode);
  int hi_degree = 3;
  int lo_degree = 2;
  // GOLD vector collection
  // Note: Vectors are exodus-ordered, not Kirby-ordered
  size_t num_p1_points = 10;
  size_t num_hi_points = num_p1_points + (num_p1_points - 1) * (hi_degree - 1);
  size_t num_lo_points = num_p1_points + (num_p1_points - 1) * (lo_degree - 1);

  std::vector<Scalar> lo_gold_in(num_lo_points);
  std::vector<Scalar> hi_gold_out(num_hi_points);
  for (size_t i = 0; i < num_p1_points; i++) {
    lo_gold_in[i]  = i;
    hi_gold_out[i] = i;
  }

  size_t idx = num_p1_points;
  for (size_t i = 0; i < num_p1_points - 1; i++) {
    for (size_t j = 0; j < (size_t)hi_degree - 1; j++) {
      hi_gold_out[idx] = i + ((double)j + 1) / hi_degree;
      idx++;
    }
  }

  idx = num_p1_points;
  for (size_t i = 0; i < num_p1_points - 1; i++) {
    for (size_t j = 0; j < (size_t)lo_degree - 1; j++) {
      lo_gold_in[idx] = i + ((double)j + 1) / lo_degree;
      idx++;
    }
  }

  TestPseudoPoisson<Scalar, LocalOrdinal, GlobalOrdinal, Node>(out, num_p1_points, hi_degree, lo_gold_in, hi_gold_out, std::string("hgrad_line_c3"), std::string("hgrad_line_c2"));
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, BuildP_PseudoPoisson_LINE_p4_to_p3, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode);
  int hi_degree = 4;
  int lo_degree = 3;
  // GOLD vector collection
  // Note: Vectors are exodus-ordered, not Kirby-ordered
  size_t num_p1_points = 10;
  size_t num_hi_points = num_p1_points + (num_p1_points - 1) * (hi_degree - 1);
  size_t num_lo_points = num_p1_points + (num_p1_points - 1) * (lo_degree - 1);

  std::vector<Scalar> lo_gold_in(num_lo_points);
  std::vector<Scalar> hi_gold_out(num_hi_points);
  for (size_t i = 0; i < num_p1_points; i++) {
    lo_gold_in[i]  = i;
    hi_gold_out[i] = i;
  }

  size_t idx = num_p1_points;
  for (size_t i = 0; i < num_p1_points - 1; i++) {
    for (size_t j = 0; j < (size_t)hi_degree - 1; j++) {
      hi_gold_out[idx] = i + ((double)j + 1) / hi_degree;
      idx++;
    }
  }

  idx = num_p1_points;
  for (size_t i = 0; i < num_p1_points - 1; i++) {
    for (size_t j = 0; j < (size_t)lo_degree - 1; j++) {
      lo_gold_in[idx] = i + ((double)j + 1) / lo_degree;
      idx++;
    }
  }

  TestPseudoPoisson<Scalar, LocalOrdinal, GlobalOrdinal, Node>(out, num_p1_points, hi_degree, lo_gold_in, hi_gold_out, std::string("hgrad_line_c4"), std::string("hgrad_line_c3"));
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, BuildP_PseudoPoisson_LINE_p4_to_p2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode);
  int hi_degree = 4;
  int lo_degree = 2;
  // GOLD vector collection
  // Note: Vectors are exodus-ordered, not Kirby-ordered
  size_t num_p1_points = 10;
  size_t num_hi_points = num_p1_points + (num_p1_points - 1) * (hi_degree - 1);
  size_t num_lo_points = num_p1_points + (num_p1_points - 1) * (lo_degree - 1);

  std::vector<Scalar> lo_gold_in(num_lo_points);
  std::vector<Scalar> hi_gold_out(num_hi_points);
  for (size_t i = 0; i < num_p1_points; i++) {
    lo_gold_in[i]  = i;
    hi_gold_out[i] = i;
  }

  size_t idx = num_p1_points;
  for (size_t i = 0; i < num_p1_points - 1; i++) {
    for (size_t j = 0; j < (size_t)hi_degree - 1; j++) {
      hi_gold_out[idx] = i + ((double)j + 1) / hi_degree;
      idx++;
    }
  }

  idx = num_p1_points;
  for (size_t i = 0; i < num_p1_points - 1; i++) {
    for (size_t j = 0; j < (size_t)lo_degree - 1; j++) {
      lo_gold_in[idx] = i + ((double)j + 1) / lo_degree;
      idx++;
    }
  }

  TestPseudoPoisson<Scalar, LocalOrdinal, GlobalOrdinal, Node>(out, num_p1_points, hi_degree, lo_gold_in, hi_gold_out, std::string("hgrad_line_c4"), std::string("hgrad_line_c2"));
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, CreatePreconditioner_p3_to_p2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode);
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
#if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#endif
#if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#endif

  typedef Scalar SC;
  typedef GlobalOrdinal GO;
  typedef LocalOrdinal LO;
  typedef Node NO;
  typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef Kokkos::DynRankView<LocalOrdinal, typename Node::device_type> FCi;

  out << "version: " << MueLu::Version() << std::endl;
  using Teuchos::RCP;
  int degree = 3;
  std::string hi_basis("hgrad_line_c3");

  Xpetra::UnderlyingLib lib          = TestHelpers::Parameters::getLib();
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  GO num_nodes = 972;
  // Build a pseudo-poisson test matrix
  FCi elem_to_node;
  RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC, LO, GO, NO>(num_nodes, degree, elem_to_node, lib);

  // Normalized RHS
  RCP<MultiVector> RHS1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  RHS1->setSeed(846930886);
  RHS1->randomize();
  Teuchos::Array<MT> norms(1);
  RHS1->norm2(norms);
  RHS1->scale(1 / norms[0]);

  // Zero initial guess
  RCP<MultiVector> X1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  X1->putScalar(Teuchos::ScalarTraits<SC>::zero());

  // ParameterList
  ParameterList Params, level0;
  Params.set("multigrid algorithm", "pcoarsen");
  //    Params.set("rap: fix zero diagonals",true);
  Params.set("pcoarsen: hi basis", hi_basis);
  Params.set("pcoarsen: lo basis", "hgrad_line_c2");
  Params.set("verbosity", "high");
  Params.set("max levels", 2);
  Params.set("coarse: max size", 100);
  level0.set("pcoarsen: element to node map", rcp(&elem_to_node, false));
  Params.set("level 0", level0);
  Params.set("verbosity", "none");

  // Build hierarchy
  RCP<Hierarchy> tH = MueLu::CreateXpetraPreconditioner<SC, LO, GO, NO>(A, Params);
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, CreatePreconditioner_p4_to_p3, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode);
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
#if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#endif
#if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#endif

  typedef Scalar SC;
  typedef GlobalOrdinal GO;
  typedef LocalOrdinal LO;
  typedef Node NO;
  typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef Kokkos::DynRankView<LocalOrdinal, typename Node::device_type> FCi;

  out << "version: " << MueLu::Version() << std::endl;
  using Teuchos::RCP;
  int degree = 4;
  std::string hi_basis("hgrad_line_c4");

  Xpetra::UnderlyingLib lib          = TestHelpers::Parameters::getLib();
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  GO num_nodes = 972;
  // Build a pseudo-poisson test matrix
  FCi elem_to_node;
  RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC, LO, GO, NO>(num_nodes, degree, elem_to_node, lib);

  // Normalized RHS
  RCP<MultiVector> RHS1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  RHS1->setSeed(846930886);
  RHS1->randomize();
  Teuchos::Array<MT> norms(1);
  RHS1->norm2(norms);
  RHS1->scale(1 / norms[0]);

  // Zero initial guess
  RCP<MultiVector> X1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  X1->putScalar(Teuchos::ScalarTraits<SC>::zero());

  // ParameterList
  ParameterList Params, level0;
  Params.set("multigrid algorithm", "pcoarsen");
  //    Params.set("rap: fix zero diagonals",true);
  Params.set("pcoarsen: hi basis", hi_basis);
  Params.set("pcoarsen: lo basis", "hgrad_line_c3");
  Params.set("verbosity", "high");
  Params.set("max levels", 2);
  Params.set("coarse: max size", 100);
  level0.set("pcoarsen: element to node map", rcp(&elem_to_node, false));
  Params.set("level 0", level0);
  Params.set("verbosity", "none");

  // Build hierarchy
  RCP<Hierarchy> tH = MueLu::CreateXpetraPreconditioner<SC, LO, GO, NO>(A, Params);
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, CreatePreconditioner_p4_to_p2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode);
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
#if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#endif
#if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#endif

  typedef Scalar SC;
  typedef GlobalOrdinal GO;
  typedef LocalOrdinal LO;
  typedef Node NO;
  typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef Kokkos::DynRankView<LocalOrdinal, typename Node::device_type> FCi;

  out << "version: " << MueLu::Version() << std::endl;
  using Teuchos::RCP;
  int degree = 4;
  std::string hi_basis("hgrad_line_c4");

  Xpetra::UnderlyingLib lib          = TestHelpers::Parameters::getLib();
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  GO num_nodes = 972;
  // Build a pseudo-poisson test matrix
  FCi elem_to_node;
  RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC, LO, GO, NO>(num_nodes, degree, elem_to_node, lib);

  // Normalized RHS
  RCP<MultiVector> RHS1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  RHS1->setSeed(846930886);
  RHS1->randomize();
  Teuchos::Array<MT> norms(1);
  RHS1->norm2(norms);
  RHS1->scale(1 / norms[0]);

  // Zero initial guess
  RCP<MultiVector> X1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  X1->putScalar(Teuchos::ScalarTraits<SC>::zero());

  // ParameterList
  ParameterList Params, level0;
  Params.set("multigrid algorithm", "pcoarsen");
  //    Params.set("rap: fix zero diagonals",true);
  Params.set("pcoarsen: hi basis", hi_basis);
  Params.set("pcoarsen: lo basis", "hgrad_line_c2");
  Params.set("verbosity", "high");
  Params.set("max levels", 2);
  Params.set("coarse: max size", 100);
  level0.set("pcoarsen: element to node map", rcp(&elem_to_node, false));
  Params.set("level 0", level0);
  Params.set("verbosity", "none");

  // Build hierarchy
  RCP<Hierarchy> tH = MueLu::CreateXpetraPreconditioner<SC, LO, GO, NO>(A, Params);
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, CreatePreconditioner_p4_to_p3_to_p2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode);
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
#if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#endif
#if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#endif

  typedef Scalar SC;
  typedef GlobalOrdinal GO;
  typedef LocalOrdinal LO;
  typedef Node NO;
  typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef Kokkos::DynRankView<LocalOrdinal, typename Node::device_type> FCi;

  out << "version: " << MueLu::Version() << std::endl;
  using Teuchos::RCP;
  int degree = 4;
  std::string hi_basis("hgrad_line_c4");

  Xpetra::UnderlyingLib lib          = TestHelpers::Parameters::getLib();
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  GO num_nodes = 972;
  // Build a pseudo-poisson test matrix
  FCi elem_to_node;
  RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC, LO, GO, NO>(num_nodes, degree, elem_to_node, lib);

  // Normalized RHS
  RCP<MultiVector> RHS1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  RHS1->setSeed(846930886);
  RHS1->randomize();
  Teuchos::Array<MT> norms(1);
  RHS1->norm2(norms);
  RHS1->scale(1 / norms[0]);

  // Zero initial guess
  RCP<MultiVector> X1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  X1->putScalar(Teuchos::ScalarTraits<SC>::zero());

  // ParameterList
  ParameterList Params, level0, level1, level2, mmm;
  Params.set("multigrid algorithm", "pcoarsen");
  Params.set("verbosity", "high");
  Params.set("max levels", 3);
  Params.set("coarse: max size", 5);

  level0.set("pcoarsen: element to node map", rcp(&elem_to_node, false));
  Params.set("level 0", level0);

  level1.set("pcoarsen: hi basis", hi_basis);
  level1.set("pcoarsen: lo basis", "hgrad_line_c3");
  Params.set("level 1", level1);

  level2.set("pcoarsen: hi basis", "hgrad_line_c3");
  level2.set("pcoarsen: lo basis", "hgrad_line_c2");
  Params.set("level 2", level2);
  Params.set("verbosity", "none");

#if 0
    // DEBUG
    ParameterList dump;
    dump.set("A","{0,1,2}");
    dump.set("P","{0,1}");
    dump.set("R","{0,1}");
    dump.set("pcoarsen: element to node map","{0,1,2}");
    Params.set("export data",dump);
#endif

  // Build hierarchy
  RCP<Hierarchy> tH = MueLu::CreateXpetraPreconditioner<SC, LO, GO, NO>(A, Params);
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, CreatePreconditioner_p2_to_p1_sa, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode);
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
#if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#endif
#if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#endif

  typedef Scalar SC;
  typedef GlobalOrdinal GO;
  typedef LocalOrdinal LO;
  typedef Node NO;
  typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef Kokkos::DynRankView<LocalOrdinal, typename Node::device_type> FCi;

  out << "version: " << MueLu::Version() << std::endl;
  using Teuchos::RCP;
  int degree = 2;
  std::string hi_basis("hgrad_line_c2");

  Xpetra::UnderlyingLib lib          = TestHelpers::Parameters::getLib();
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  GO num_nodes = 972;
  // Build a pseudo-poisson test matrix
  FCi elem_to_node;
  RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC, LO, GO, NO>(num_nodes, degree, elem_to_node, lib);

  // Normalized RHS
  RCP<MultiVector> RHS1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  RHS1->setSeed(846930886);
  RHS1->randomize();
  Teuchos::Array<MT> norms(1);
  RHS1->norm2(norms);
  RHS1->scale(1 / norms[0]);

  // Zero initial guess
  RCP<MultiVector> X1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  X1->putScalar(Teuchos::ScalarTraits<SC>::zero());

  // ParameterList
  ParameterList Params, level0, level1, level2;
  Params.set("multigrid algorithm", "pcoarsen");
  Params.set("rap: fix zero diagonals", true);
  Params.set("verbosity", "high");
  Params.set("max levels", 3);
  Params.set("coarse: max size", 100);

  level0.set("pcoarsen: element to node map", rcp(&elem_to_node, false));
  Params.set("level 0", level0);

  level1.set("multigrid algorithm", "pcoarsen");
  level1.set("pcoarsen: hi basis", hi_basis);
  level1.set("pcoarsen: lo basis", "hgrad_line_c1");
  Params.set("level 1", level1);

  level2.set("multigrid algorithm", "sa");
  Params.set("level 2", level2);
  Params.set("verbosity", "none");

  // Build hierarchy
  RCP<Hierarchy> tH = MueLu::CreateXpetraPreconditioner<SC, LO, GO, NO>(A, Params);
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, CreatePreconditioner_p3_to_p2_to_p1_sa_manual, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  ;
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode);
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
#if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#endif
#if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#endif

  typedef Scalar SC;
  typedef GlobalOrdinal GO;
  typedef LocalOrdinal LO;
  typedef Node NO;
  typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef Kokkos::DynRankView<LocalOrdinal, typename Node::device_type> FCi;

  out << "version: " << MueLu::Version() << std::endl;
  using Teuchos::RCP;
  int degree = 3;
  std::string hi_basis("hgrad_line_c3");

  Xpetra::UnderlyingLib lib          = TestHelpers::Parameters::getLib();
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  GO num_nodes = 972;
  // Build a pseudo-poisson test matrix
  FCi elem_to_node;
  RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC, LO, GO, NO>(num_nodes, degree, elem_to_node, lib);

  // Normalized RHS
  RCP<MultiVector> RHS1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  RHS1->setSeed(846930886);
  RHS1->randomize();
  Teuchos::Array<MT> norms(1);
  RHS1->norm2(norms);
  RHS1->scale(1 / norms[0]);

  // Zero initial guess
  RCP<MultiVector> X1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  X1->putScalar(Teuchos::ScalarTraits<SC>::zero());

  // ParameterList
  ParameterList Params, level0, level1, level2, level3;
  Params.set("multigrid algorithm", "pcoarsen");
  Params.set("rap: fix zero diagonals", true);
  Params.set("verbosity", "high");
  Params.set("max levels", 4);
  Params.set("coarse: max size", 100);

  level0.set("pcoarsen: element to node map", rcp(&elem_to_node, false));
  Params.set("level 0", level0);

  level1.set("pcoarsen: hi basis", "hgrad_line_c3");
  level1.set("pcoarsen: lo basis", "hgrad_line_c2");
  Params.set("level 1", level1);

  level2.set("pcoarsen: hi basis", "hgrad_line_c2");
  level2.set("pcoarsen: lo basis", "hgrad_line_c1");
  Params.set("level 2", level2);

  level3.set("multigrid algorithm", "sa");
  Params.set("level 3", level3);
  Params.set("verbosity", "none");

  // Build hierarchy
  RCP<Hierarchy> tH = MueLu::CreateXpetraPreconditioner<SC, LO, GO, NO>(A, Params);
}

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, CreatePreconditioner_p3_to_p2_to_p1_sa_schedule, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode);
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
#if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#endif
#if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#endif

  typedef Scalar SC;
  typedef GlobalOrdinal GO;
  typedef LocalOrdinal LO;
  typedef Node NO;
  typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef Kokkos::DynRankView<LocalOrdinal, typename Node::device_type> FCi;

  out << "version: " << MueLu::Version() << std::endl;
  using Teuchos::RCP;
  int degree = 3;
  std::string hi_basis("hgrad_line_c3");

  Xpetra::UnderlyingLib lib          = TestHelpers::Parameters::getLib();
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  GO num_nodes = 972;
  // Build a pseudo-poisson test matrix
  FCi elem_to_node;
  RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC, LO, GO, NO>(num_nodes, degree, elem_to_node, lib);

  // Normalized RHS
  RCP<MultiVector> RHS1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  RHS1->setSeed(846930886);
  RHS1->randomize();
  Teuchos::Array<MT> norms(1);
  RHS1->norm2(norms);
  RHS1->scale(1 / norms[0]);

  // Zero initial guess
  RCP<MultiVector> X1 = MultiVectorFactory::Build(A->getRowMap(), 1);
  X1->putScalar(Teuchos::ScalarTraits<SC>::zero());

  // ParameterList
  ParameterList Params, level0;
  Params.set("multigrid algorithm", "pcoarsen");
  Params.set("rap: fix zero diagonals", true);
  Params.set("verbosity", "high");
  Params.set("max levels", 5);
  Params.set("coarse: max size", 10);

  level0.set("pcoarsen: element to node map", rcp(&elem_to_node, false));
  Params.set("level 0", level0);

  Params.set("pcoarsen: element", "hgrad_line_c");
  Params.set("pcoarsen: schedule", "{3,2,1}");
  Params.set("verbosity", "none");

  // Build hierarchy
  RCP<Hierarchy> tH = MueLu::CreateXpetraPreconditioner<SC, LO, GO, NO>(A, Params);
}

/*********************************************************************************************************************/

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                                                            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, GetP1NodeInHi, Scalar, LO, GO, Node)                                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, BasisFactory, Scalar, LO, GO, Node)                                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, BuildLoElemToNode, Scalar, LO, GO, Node)                                 \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, BuildLoElemToNodeWithDirichlet, Scalar, LO, GO, Node)                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, GenerateColMapFromImport, Scalar, LO, GO, Node)                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, BuildP_PseudoPoisson_p2, Scalar, LO, GO, Node)                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, BuildP_PseudoPoisson_p3, Scalar, LO, GO, Node)                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, BuildP_PseudoPoisson_p4, Scalar, LO, GO, Node)                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p2, Scalar, LO, GO, Node)                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p3, Scalar, LO, GO, Node)                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p4, Scalar, LO, GO, Node)                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_LINE_Equispaced, Scalar, LO, GO, Node)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_QUAD_Equispaced, Scalar, LO, GO, Node)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_QUAD_Spectral, Scalar, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_HEX_Equispaced, Scalar, LO, GO, Node)   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_HEX_Spectral, Scalar, LO, GO, Node)     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, GenerateLoNodeInHighViaGIDs_QUAD_pn_to_p1, Scalar, LO, GO, Node)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, BuildLoElemToNodeViaRepresentatives_QUAD_pn_to_p1, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, BuildP_PseudoPoisson_LINE_p3_to_p2, Scalar, LO, GO, Node)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, BuildP_PseudoPoisson_LINE_p4_to_p3, Scalar, LO, GO, Node)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, BuildP_PseudoPoisson_LINE_p4_to_p2, Scalar, LO, GO, Node)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p3_to_p2, Scalar, LO, GO, Node)                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p4_to_p3, Scalar, LO, GO, Node)                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p4_to_p2, Scalar, LO, GO, Node)                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p4_to_p3_to_p2, Scalar, LO, GO, Node)               \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p2_to_p1_sa, Scalar, LO, GO, Node)                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p3_to_p2_to_p1_sa_manual, Scalar, LO, GO, Node)     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p3_to_p2_to_p1_sa_schedule, Scalar, LO, GO, Node)   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, BuildSampleElementToNodeMapThreeElementQuad, Scalar, LO, GO, Node)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, UniqueNumbering, Scalar, LO, GO, Node)                                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, FindSeeds_Equispaced_Line, Scalar, LO, GO, Node)                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, FindSeeds_Spectral_Line, Scalar, LO, GO, Node)                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, FindSeeds_Equispaced_Quad, Scalar, LO, GO, Node)                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, FindSeeds_Spectral_Quad, Scalar, LO, GO, Node)                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, FindSeeds_Equispaced_Hex, Scalar, LO, GO, Node)                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, FindSeeds_Spectral_Hex, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
#endif
