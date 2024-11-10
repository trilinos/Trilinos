// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_GLOBALLEXICOGRAPHICINDEXMANAGER_DEF_HPP_
#define MUELU_GLOBALLEXICOGRAPHICINDEXMANAGER_DEF_HPP_

#include <MueLu_GlobalLexicographicIndexManager_decl.hpp>
#include <Xpetra_MapFactory.hpp>

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    GlobalLexicographicIndexManager(const RCP<const Teuchos::Comm<int> > comm, const bool coupled,
                                    const int NumDimensions, const int interpolationOrder,
                                    const Array<GO> GFineNodesPerDir,
                                    const Array<LO> LFineNodesPerDir, const Array<LO> CoarseRate,
                                    const GO MinGlobalIndex)
  : IndexManager(comm, coupled, false, NumDimensions, interpolationOrder, GFineNodesPerDir, LFineNodesPerDir) {
  // Load coarse rate, being careful about formating.
  for (int dim = 0; dim < 3; ++dim) {
    if (dim < this->numDimensions) {
      if (CoarseRate.size() == 1) {
        this->coarseRate[dim] = CoarseRate[0];
      } else if (CoarseRate.size() == this->numDimensions) {
        this->coarseRate[dim] = CoarseRate[dim];
      }
    } else {
      this->coarseRate[dim] = 1;
    }
  }

  {
    GO tmp                = 0;
    this->startIndices[2] = MinGlobalIndex / (this->gFineNodesPerDir[1] * this->gFineNodesPerDir[0]);
    tmp                   = MinGlobalIndex % (this->gFineNodesPerDir[1] * this->gFineNodesPerDir[0]);
    this->startIndices[1] = tmp / this->gFineNodesPerDir[0];
    this->startIndices[0] = tmp % this->gFineNodesPerDir[0];

    for (int dim = 0; dim < 3; ++dim) {
      this->startIndices[dim + 3] = this->startIndices[dim] + this->lFineNodesPerDir[dim] - 1;
    }
  }

  this->computeMeshParameters();
  computeGlobalCoarseParameters();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    computeGlobalCoarseParameters() {
  this->gNumCoarseNodes10 = this->gCoarseNodesPerDir[0] * this->gCoarseNodesPerDir[1];
  this->gNumCoarseNodes   = this->gNumCoarseNodes10 * this->gCoarseNodesPerDir[2];
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getGhostedNodesData(const RCP<const Map> fineMap,
                        Array<LO>& ghostedNodeCoarseLIDs, Array<int>& ghostedNodeCoarsePIDs, Array<GO>& ghostedNodeCoarseGIDs) const {
  ghostedNodeCoarseLIDs.resize(this->getNumLocalGhostedNodes());
  ghostedNodeCoarsePIDs.resize(this->getNumLocalGhostedNodes());
  ghostedNodeCoarseGIDs.resize(this->numGhostedNodes);

  // Find the GIDs, LIDs and PIDs of the coarse points on the fine mesh and coarse
  // mesh as this data will be used to fill vertex2AggId and procWinner vectors.
  Array<GO> lCoarseNodeCoarseGIDs(this->lNumCoarseNodes),
      lCoarseNodeFineGIDs(this->lNumCoarseNodes);
  Array<GO> ghostedCoarseNodeFineGIDs(this->numGhostedNodes);
  Array<LO> ghostedCoarseNodeCoarseIndices(3), ghostedCoarseNodeFineIndices(3), ijk(3);
  LO currentIndex = -1, currentCoarseIndex = -1;
  for (ijk[2] = 0; ijk[2] < this->ghostedNodesPerDir[2]; ++ijk[2]) {
    for (ijk[1] = 0; ijk[1] < this->ghostedNodesPerDir[1]; ++ijk[1]) {
      for (ijk[0] = 0; ijk[0] < this->ghostedNodesPerDir[0]; ++ijk[0]) {
        currentIndex                        = ijk[2] * this->numGhostedNodes10 + ijk[1] * this->ghostedNodesPerDir[0] + ijk[0];
        ghostedCoarseNodeCoarseIndices[0]   = this->startGhostedCoarseNode[0] + ijk[0];
        ghostedCoarseNodeCoarseIndices[1]   = this->startGhostedCoarseNode[1] + ijk[1];
        ghostedCoarseNodeCoarseIndices[2]   = this->startGhostedCoarseNode[2] + ijk[2];
        GO myCoarseGID                      = ghostedCoarseNodeCoarseIndices[0] + ghostedCoarseNodeCoarseIndices[1] * this->gCoarseNodesPerDir[0] + ghostedCoarseNodeCoarseIndices[2] * this->gNumCoarseNodes10;
        ghostedNodeCoarseGIDs[currentIndex] = myCoarseGID;
        GO myGID = 0, factor[3] = {};
        factor[2] = this->gNumFineNodes10;
        factor[1] = this->gFineNodesPerDir[0];
        factor[0] = 1;
        for (int dim = 0; dim < 3; ++dim) {
          if (dim < this->numDimensions) {
            if (this->startIndices[dim] - this->offsets[dim] + ijk[dim] * this->coarseRate[dim] < this->gFineNodesPerDir[dim] - 1) {
              myGID += (this->startIndices[dim] - this->offsets[dim] + ijk[dim] * this->coarseRate[dim]) * factor[dim];
            } else {
              myGID += (this->startIndices[dim] - this->offsets[dim] + (ijk[dim] - 1) * this->coarseRate[dim] + this->endRate[dim]) * factor[dim];
            }
          }
        }
        // lbv 02-08-2018:
        // This check is simplistic and should be replaced by a condition that checks
        // if the local tuple of the current index is wihin the range of local nodes
        // or not in the range of ghosted nodes.
        if ((!this->ghostInterface[0] || ijk[0] != 0) &&
            (!this->ghostInterface[2] || ijk[1] != 0) &&
            (!this->ghostInterface[4] || ijk[2] != 0) &&
            (!this->ghostInterface[1] || ijk[0] != this->ghostedNodesPerDir[0] - 1) &&
            (!this->ghostInterface[3] || ijk[1] != this->ghostedNodesPerDir[1] - 1) &&
            (!this->ghostInterface[5] || ijk[2] != this->ghostedNodesPerDir[2] - 1)) {
          // this->getGhostedNodeFineLID(ijk[0], ijk[1], ijk[2], coarseNodeFineLID);
          if (this->interpolationOrder_ == 0) {
            currentCoarseIndex = 0;
            if (this->ghostInterface[4]) {
              currentCoarseIndex += (ijk[2] - 1) * this->lNumCoarseNodes10;
            } else {
              currentCoarseIndex += ijk[2] * this->lNumCoarseNodes10;
            }
            if (this->ghostInterface[2]) {
              currentCoarseIndex += (ijk[1] - 1) * this->getLocalCoarseNodesInDir(0);
            } else {
              currentCoarseIndex += ijk[1] * this->getLocalCoarseNodesInDir(0);
            }
            if (this->ghostInterface[0]) {
              currentCoarseIndex += ijk[0] - 1;
            } else {
              currentCoarseIndex += ijk[0];
            }
          } else {
            this->getGhostedNodeCoarseLID(ijk[0], ijk[1], ijk[2], currentCoarseIndex);
          }

          lCoarseNodeCoarseGIDs[currentCoarseIndex] = myCoarseGID;
          lCoarseNodeFineGIDs[currentCoarseIndex]   = myGID;
        }
        ghostedCoarseNodeFineGIDs[currentIndex] = myGID;
      }
    }
  }

  RCP<const Map> coarseMap = Xpetra::MapFactory<LO, GO, NO>::Build(fineMap->lib(),
                                                                   this->gNumCoarseNodes,
                                                                   lCoarseNodeCoarseGIDs(),
                                                                   fineMap->getIndexBase(),
                                                                   fineMap->getComm());

  coarseMap->getRemoteIndexList(ghostedNodeCoarseGIDs(),
                                ghostedNodeCoarsePIDs(),
                                ghostedNodeCoarseLIDs());

}  // End getGhostedMeshData

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodesData(const RCP<const Map> fineCoordinatesMap,
                       Array<GO>& coarseNodeCoarseGIDs,
                       Array<GO>& coarseNodeFineGIDs) const {
  // Allocate sufficient storage space for outputs
  coarseNodeCoarseGIDs.resize(this->getNumLocalCoarseNodes());
  coarseNodeFineGIDs.resize(this->getNumLocalCoarseNodes());

  // Load all the GIDs on the fine mesh
  ArrayView<const GO> fineNodeGIDs = fineCoordinatesMap->getLocalElementList();

  Array<GO> coarseStartIndices(3);
  GO tmp;
  for (int dim = 0; dim < 3; ++dim) {
    coarseStartIndices[dim] = this->startIndices[dim] / this->coarseRate[dim];
    tmp                     = this->startIndices[dim] % this->coarseRate[dim];
    if (tmp > 0) {
      ++coarseStartIndices[dim];
    }
  }

  // Extract the fine LIDs of the coarse nodes and store the corresponding GIDs
  LO fineLID;
  Array<LO> lCoarseIndices(3);
  Array<GO> gCoarseIndices(3);
  for (LO coarseLID = 0; coarseLID < this->getNumLocalCoarseNodes(); ++coarseLID) {
    this->getCoarseNodeLocalTuple(coarseLID,
                                  lCoarseIndices[0],
                                  lCoarseIndices[1],
                                  lCoarseIndices[2]);
    getCoarseNodeFineLID(lCoarseIndices[0], lCoarseIndices[1], lCoarseIndices[2], fineLID);
    coarseNodeFineGIDs[coarseLID] = fineNodeGIDs[fineLID];

    // Get Coarse Global IJK
    for (int dim = 0; dim < 3; dim++) {
      gCoarseIndices[dim] = coarseStartIndices[dim] + lCoarseIndices[dim];
    }
    getCoarseNodeGID(gCoarseIndices[0],
                     gCoarseIndices[1],
                     gCoarseIndices[2],
                     coarseNodeCoarseGIDs[coarseLID]);
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
std::vector<std::vector<GlobalOrdinal> > GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseMeshData() const {
  std::vector<std::vector<GO> > coarseMeshData;
  return coarseMeshData;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getFineNodeGlobalTuple(const GO myGID, GO& i, GO& j, GO& k) const {
  GO tmp;
  k   = myGID / this->gNumFineNodes10;
  tmp = myGID % this->gNumFineNodes10;
  j   = tmp / this->gFineNodesPerDir[0];
  i   = tmp % this->gFineNodesPerDir[0];
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getFineNodeLocalTuple(const LO myLID, LO& i, LO& j, LO& k) const {
  LO tmp;
  k   = myLID / this->lNumFineNodes10;
  tmp = myLID % this->lNumFineNodes10;
  j   = tmp / this->lFineNodesPerDir[0];
  i   = tmp % this->lFineNodesPerDir[0];
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getFineNodeGhostedTuple(const LO myLID, LO& i, LO& j, LO& k) const {
  LO tmp;
  k   = myLID / this->lNumFineNodes10;
  tmp = myLID % this->lNumFineNodes10;
  j   = tmp / this->lFineNodesPerDir[0];
  i   = tmp % this->lFineNodesPerDir[0];

  k += this->offsets[2];
  j += this->offsets[1];
  i += this->offsets[0];
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getFineNodeGID(const GO i, const GO j, const GO k, GO& myGID) const {
  myGID = k * this->gNumFineNodes10 + j * this->gFineNodesPerDir[0] + i;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getFineNodeLID(const LO i, const LO j, const LO k, LO& myLID) const {
  myLID = k * this->lNumFineNodes10 + j * this->lFineNodesPerDir[0] + i;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodeGlobalTuple(const GO myGID, GO& i, GO& j, GO& k) const {
  GO tmp;
  k   = myGID / this->gNumCoarseNodes10;
  tmp = myGID % this->gNumCoarseNodes10;
  j   = tmp / this->gCoarseNodesPerDir[0];
  i   = tmp % this->gCoarseNodesPerDir[0];
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodeLocalTuple(const LO myLID, LO& i, LO& j, LO& k) const {
  LO tmp;
  k   = myLID / this->lNumCoarseNodes10;
  tmp = myLID % this->lNumCoarseNodes10;
  j   = tmp / this->lCoarseNodesPerDir[0];
  i   = tmp % this->lCoarseNodesPerDir[0];
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodeGID(const GO i, const GO j, const GO k, GO& myGID) const {
  myGID = k * this->gNumCoarseNodes10 + j * this->gCoarseNodesPerDir[0] + i;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodeLID(const LO i, const LO j, const LO k, LO& myLID) const {
  myLID = k * this->lNumCoarseNodes10 + j * this->lCoarseNodesPerDir[0] + i;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodeGhostedLID(const LO i, const LO j, const LO k, LO& myLID) const {
  myLID = k * this->numGhostedNodes10 + j * this->ghostedNodesPerDir[0] + i;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodeFineLID(const LO i, const LO j, const LO k, LO& myLID) const {
  // Assumptions: (i,j,k) is a tuple on the coarse mesh
  //              myLID is the corresponding local ID on the fine mesh
  const LO multiplier[3] = {1, this->lFineNodesPerDir[0], this->lNumFineNodes10};
  const LO indices[3]    = {i, j, k};

  myLID = 0;
  for (int dim = 0; dim < 3; ++dim) {
    if ((indices[dim] == this->getLocalCoarseNodesInDir(dim) - 1) && this->meshEdge[2 * dim + 1]) {
      // We are dealing with the last node on the mesh in direction dim
      // so we can simply use the number of nodes on the fine mesh in that direction
      myLID += (this->getLocalFineNodesInDir(dim) - 1) * multiplier[dim];
    } else {
      myLID += (indices[dim] * this->getCoarseningRate(dim) + this->getCoarseNodeOffset(dim)) * multiplier[dim];
    }
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getGhostedNodeFineLID(const LO i, const LO j, const LO k, LO& myLID) const {
  LO itmp = i - (this->offsets[0] > 0 ? 1 : 0);
  LO jtmp = j - (this->offsets[1] > 0 ? 1 : 0);
  LO ktmp = k - (this->offsets[2] > 0 ? 1 : 0);
  myLID   = 0;
  if (ktmp * this->coarseRate[2] < this->lFineNodesPerDir[2]) {
    myLID += ktmp * this->coarseRate[2] * this->lNumCoarseNodes10;
  } else {
    myLID += (this->lFineNodesPerDir[2] - 1) * this->lNumCoarseNodes10;
  }

  if (jtmp * this->coarseRate[1] < this->lFineNodesPerDir[1]) {
    myLID += jtmp * this->coarseRate[1] * this->lFineNodesPerDir[0];
  } else {
    myLID += (this->lFineNodesPerDir[1] - 1) * this->lFineNodesPerDir[1];
  }

  if (itmp * this->coarseRate[0] < this->lFineNodesPerDir[0]) {
    myLID += itmp * this->coarseRate[0];
  } else {
    myLID += this->lFineNodesPerDir[0] - 1;
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getGhostedNodeCoarseLID(const LO i, const LO j, const LO k, LO& myLID) const {
  LO itmp = i - (this->offsets[0] > 0 ? 1 : 0);
  LO jtmp = j - (this->offsets[1] > 0 ? 1 : 0);
  LO ktmp = k - (this->offsets[2] > 0 ? 1 : 0);
  myLID   = ktmp * this->lNumCoarseNodes10 + jtmp * this->lCoarseNodesPerDir[0] + itmp;
}

}  // namespace MueLu

#endif /* MUELU_GLOBALLEXICOGRAPHICINDEXMANAGER_DEF_HPP_ */
