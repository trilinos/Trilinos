// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_LOCALLEXICOGRAPHICINDEXMANAGER_DEF_HPP_
#define MUELU_LOCALLEXICOGRAPHICINDEXMANAGER_DEF_HPP_

#include <MueLu_LocalLexicographicIndexManager_decl.hpp>
#include <Xpetra_MapFactory.hpp>

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    LocalLexicographicIndexManager(const RCP<const Teuchos::Comm<int> > comm, const bool coupled,
                                   const int NumDimensions, const int interpolationOrder,
                                   const int MyRank, const int NumRanks,
                                   const Array<GO> GFineNodesPerDir, const Array<LO> LFineNodesPerDir,
                                   const Array<LO> CoarseRate, const Array<GO> MeshData)
  : IndexManager(comm, coupled, false, NumDimensions, interpolationOrder, GFineNodesPerDir, LFineNodesPerDir)
  , myRank(MyRank)
  , numRanks(NumRanks) {
  // Allocate data based on user input
  meshData.resize(numRanks);
  rankIndices.resize(numRanks);
  coarseMeshData.resize(numRanks);

  // Load coarse rate, being careful about formating
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

  // Load meshData for local lexicographic case
  for (int rank = 0; rank < numRanks; ++rank) {
    meshData[rank].resize(10);
    for (int entry = 0; entry < 10; ++entry) {
      meshData[rank][entry] = MeshData[10 * rank + entry];
    }
  }

  if (this->coupled_) {
    myBlock = meshData[myRank][2];
    sortLocalLexicographicData();
  }

  // Start simple parameter calculation
  myRankIndex = rankIndices[myRank];
  for (int dim = 0; dim < 3; ++dim) {
    this->startIndices[dim]     = meshData[myRankIndex][2 * dim + 3];
    this->startIndices[dim + 3] = meshData[myRankIndex][2 * dim + 4];
  }

  this->computeMeshParameters();
  computeGlobalCoarseParameters();
  computeCoarseLocalLexicographicData();
}  // Constructor

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    computeGlobalCoarseParameters() {
  this->gNumCoarseNodes10 = this->gCoarseNodesPerDir[0] * this->gCoarseNodesPerDir[1];
  this->gNumCoarseNodes   = this->gNumCoarseNodes10 * this->gCoarseNodesPerDir[2];
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getGhostedNodesData(const RCP<const Map> /* fineMap */,
                        Array<LO>& ghostedNodeCoarseLIDs,
                        Array<int>& ghostedNodeCoarsePIDs,
                        Array<GO>& ghostedNodeCoarseGIDs) const {
  // First we allocated memory for the outputs
  ghostedNodeCoarseLIDs.resize(this->getNumLocalGhostedNodes());
  ghostedNodeCoarsePIDs.resize(this->getNumLocalGhostedNodes());
  ghostedNodeCoarseGIDs.resize(this->numGhostedNodes);

  // Now the tricky part starts, the coarse nodes / ghosted coarse nodes need to be imported.
  // This requires finding what their GID on the fine mesh is. They need to be ordered
  // lexicographically to allow for fast sweeps through the mesh.

  // We loop over all ghosted coarse nodes by increasing global lexicographic order
  Array<LO> ghostedCoarseNodeCoarseIndices(3), ghostedCoarseNodeFineIndices(3);
  Array<LO> lCoarseNodeCoarseIndices(3);
  Array<GO> lCoarseNodeCoarseGIDs(this->lNumCoarseNodes);
  LO currentIndex = -1, countCoarseNodes = 0;
  for (int k = 0; k < this->ghostedNodesPerDir[2]; ++k) {
    for (int j = 0; j < this->ghostedNodesPerDir[1]; ++j) {
      for (int i = 0; i < this->ghostedNodesPerDir[0]; ++i) {
        currentIndex                      = k * this->numGhostedNodes10 + j * this->ghostedNodesPerDir[0] + i;
        ghostedCoarseNodeCoarseIndices[0] = this->startGhostedCoarseNode[0] + i;
        ghostedCoarseNodeFineIndices[0]   = ghostedCoarseNodeCoarseIndices[0] * this->coarseRate[0];
        if (ghostedCoarseNodeFineIndices[0] > this->gFineNodesPerDir[0] - 1) {
          ghostedCoarseNodeFineIndices[0] = this->gFineNodesPerDir[0] - 1;
        }
        ghostedCoarseNodeCoarseIndices[1] = this->startGhostedCoarseNode[1] + j;
        ghostedCoarseNodeFineIndices[1]   = ghostedCoarseNodeCoarseIndices[1] * this->coarseRate[1];
        if (ghostedCoarseNodeFineIndices[1] > this->gFineNodesPerDir[1] - 1) {
          ghostedCoarseNodeFineIndices[1] = this->gFineNodesPerDir[1] - 1;
        }
        ghostedCoarseNodeCoarseIndices[2] = this->startGhostedCoarseNode[2] + k;
        ghostedCoarseNodeFineIndices[2]   = ghostedCoarseNodeCoarseIndices[2] * this->coarseRate[2];
        if (ghostedCoarseNodeFineIndices[2] > this->gFineNodesPerDir[2] - 1) {
          ghostedCoarseNodeFineIndices[2] = this->gFineNodesPerDir[2] - 1;
        }

        GO myGID = -1, myCoarseGID = -1;
        LO myLID = -1, myPID = -1, myCoarseLID = -1;
        getGIDLocalLexicographic(i, j, k, ghostedCoarseNodeFineIndices, myGID, myPID, myLID);

        int rankIndex = rankIndices[myPID];
        for (int dim = 0; dim < 3; ++dim) {
          if (dim < this->numDimensions) {
            lCoarseNodeCoarseIndices[dim] = ghostedCoarseNodeCoarseIndices[dim] - coarseMeshData[rankIndex][3 + 2 * dim];
          }
        }
        LO myRankIndexCoarseNodesInDir0 = coarseMeshData[rankIndex][4] - coarseMeshData[rankIndex][3] + 1;
        LO myRankIndexCoarseNodes10     = (coarseMeshData[rankIndex][6] - coarseMeshData[rankIndex][5] + 1) * myRankIndexCoarseNodesInDir0;
        myCoarseLID                     = lCoarseNodeCoarseIndices[2] * myRankIndexCoarseNodes10 + lCoarseNodeCoarseIndices[1] * myRankIndexCoarseNodesInDir0 + lCoarseNodeCoarseIndices[0];
        myCoarseGID                     = myCoarseLID + coarseMeshData[rankIndex][9];

        ghostedNodeCoarseLIDs[currentIndex] = myCoarseLID;
        ghostedNodeCoarsePIDs[currentIndex] = myPID;
        ghostedNodeCoarseGIDs[currentIndex] = myCoarseGID;

        if (myPID == myRank) {
          lCoarseNodeCoarseGIDs[countCoarseNodes] = myCoarseGID;
          ++countCoarseNodes;
        }
      }
    }
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodesData(const RCP<const Map> fineCoordinatesMap,
                       Array<GO>& coarseNodeCoarseGIDs,
                       Array<GO>& coarseNodeFineGIDs) const {
  // Allocate sufficient storage space for outputs
  coarseNodeCoarseGIDs.resize(this->getNumLocalCoarseNodes());
  coarseNodeFineGIDs.resize(this->getNumLocalCoarseNodes());

  // Load all the GIDs on the fine mesh
  ArrayView<const GO> fineNodeGIDs = fineCoordinatesMap->getLocalElementList();

  Array<GO> coarseStartIndices(3);
  for (int dim = 0; dim < 3; ++dim) {
    coarseStartIndices[dim] = this->coarseMeshData[myRankIndex][2 * dim + 3];
  }

  // Extract the fine LIDs of the coarse nodes and store the corresponding GIDs
  LO fineLID;
  for (LO coarseLID = 0; coarseLID < this->getNumLocalCoarseNodes(); ++coarseLID) {
    Array<LO> coarseIndices(3), fineIndices(3), gCoarseIndices(3);
    this->getCoarseNodeLocalTuple(coarseLID,
                                  coarseIndices[0],
                                  coarseIndices[1],
                                  coarseIndices[2]);
    getCoarseNodeFineLID(coarseIndices[0], coarseIndices[1], coarseIndices[2], fineLID);
    coarseNodeFineGIDs[coarseLID] = fineNodeGIDs[fineLID];

    LO myRankIndexCoarseNodesInDir0 = coarseMeshData[myRankIndex][4] - coarseMeshData[myRankIndex][3] + 1;
    LO myRankIndexCoarseNodes10     = (coarseMeshData[myRankIndex][6] - coarseMeshData[myRankIndex][5] + 1) * myRankIndexCoarseNodesInDir0;
    LO myCoarseLID                  = coarseIndices[2] * myRankIndexCoarseNodes10 + coarseIndices[1] * myRankIndexCoarseNodesInDir0 + coarseIndices[0];
    GO myCoarseGID                  = myCoarseLID + coarseMeshData[myRankIndex][9];
    coarseNodeCoarseGIDs[coarseLID] = myCoarseGID;
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getGIDLocalLexicographic(const LO iGhosted, const LO jGhosted, const LO kGhosted,
                             const Array<LO> coarseNodeFineIndices,
                             GO& myGID, LO& myPID, LO& myLID) const {
  LO ni = -1, nj = -1, li = -1, lj = -1, lk = -1;
  LO myRankGuess = myRankIndex;
  // We try to make a logical guess as to which PID owns the current coarse node
  if (iGhosted == 0 && this->ghostInterface[0]) {
    --myRankGuess;
  } else if ((iGhosted == this->ghostedNodesPerDir[0] - 1) && this->ghostInterface[1]) {
    ++myRankGuess;
  }
  if (jGhosted == 0 && this->ghostInterface[2]) {
    myRankGuess -= pi;
  } else if ((jGhosted == this->ghostedNodesPerDir[1] - 1) && this->ghostInterface[3]) {
    myRankGuess += pi;
  }
  if (kGhosted == 0 && this->ghostInterface[4]) {
    myRankGuess -= pj * pi;
  } else if ((kGhosted == this->ghostedNodesPerDir[2] - 1) && this->ghostInterface[5]) {
    myRankGuess += pj * pi;
  }
  if (coarseNodeFineIndices[0] >= meshData[myRankGuess][3] && coarseNodeFineIndices[0] <= meshData[myRankGuess][4] && coarseNodeFineIndices[1] >= meshData[myRankGuess][5] && coarseNodeFineIndices[1] <= meshData[myRankGuess][6] && coarseNodeFineIndices[2] >= meshData[myRankGuess][7] && coarseNodeFineIndices[2] <= meshData[myRankGuess][8] && myRankGuess < numRanks - 1) {
    myPID = meshData[myRankGuess][0];
    ni    = meshData[myRankGuess][4] - meshData[myRankGuess][3] + 1;
    nj    = meshData[myRankGuess][6] - meshData[myRankGuess][5] + 1;
    li    = coarseNodeFineIndices[0] - meshData[myRankGuess][3];
    lj    = coarseNodeFineIndices[1] - meshData[myRankGuess][5];
    lk    = coarseNodeFineIndices[2] - meshData[myRankGuess][7];
    myLID = lk * nj * ni + lj * ni + li;
    myGID = meshData[myRankGuess][9] + myLID;
  } else {  // The guess failed, let us use the heavy artilery: std::find_if()
    // It could be interesting to monitor how many times this branch of the code gets
    // used as it is far more expensive than the above one...
    auto nodeRank = std::find_if(myBlockStart, myBlockEnd,
                                 [coarseNodeFineIndices](const std::vector<GO>& vec) {
                                   if (coarseNodeFineIndices[0] >= vec[3] && coarseNodeFineIndices[0] <= vec[4] && coarseNodeFineIndices[1] >= vec[5] && coarseNodeFineIndices[1] <= vec[6] && coarseNodeFineIndices[2] >= vec[7] && coarseNodeFineIndices[2] <= vec[8]) {
                                     return true;
                                   } else {
                                     return false;
                                   }
                                 });
    myPID         = (*nodeRank)[0];
    ni            = (*nodeRank)[4] - (*nodeRank)[3] + 1;
    nj            = (*nodeRank)[6] - (*nodeRank)[5] + 1;
    li            = coarseNodeFineIndices[0] - (*nodeRank)[3];
    lj            = coarseNodeFineIndices[1] - (*nodeRank)[5];
    lk            = coarseNodeFineIndices[2] - (*nodeRank)[7];
    myLID         = lk * nj * ni + lj * ni + li;
    myGID         = (*nodeRank)[9] + myLID;
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    sortLocalLexicographicData() {
  std::sort(meshData.begin(), meshData.end(),
            [](const std::vector<GO>& a, const std::vector<GO>& b) -> bool {
              // The below function sorts ranks by blockID, kmin, jmin and imin
              if (a[2] < b[2]) {
                return true;
              } else if (a[2] == b[2]) {
                if (a[7] < b[7]) {
                  return true;
                } else if (a[7] == b[7]) {
                  if (a[5] < b[5]) {
                    return true;
                  } else if (a[5] == b[5]) {
                    if (a[3] < b[3]) {
                      return true;
                    }
                  }
                }
              }
              return false;
            });

  numBlocks = meshData[numRanks - 1][2] + 1;
  // Find the range of the current block
  myBlockStart = std::lower_bound(meshData.begin(), meshData.end(), myBlock - 1,
                                  [](const std::vector<GO>& vec, const GO val) -> bool {
                                    return (vec[2] < val) ? true : false;
                                  });
  myBlockEnd   = std::upper_bound(meshData.begin(), meshData.end(), myBlock,
                                  [](const GO val, const std::vector<GO>& vec) -> bool {
                                  return (val < vec[2]) ? true : false;
                                });
  // Assuming that i,j,k and ranges are split in pi, pj and pk processors
  // we search for these numbers as they will allow us to find quickly the PID of processors
  // owning ghost nodes.
  auto myKEnd = std::upper_bound(myBlockStart, myBlockEnd, (*myBlockStart)[3],
                                 [](const GO val, const std::vector<GO>& vec) -> bool {
                                   return (val < vec[7]) ? true : false;
                                 });
  auto myJEnd = std::upper_bound(myBlockStart, myKEnd, (*myBlockStart)[3],
                                 [](const GO val, const std::vector<GO>& vec) -> bool {
                                   return (val < vec[5]) ? true : false;
                                 });
  pi          = std::distance(myBlockStart, myJEnd);
  pj          = std::distance(myBlockStart, myKEnd) / pi;
  pk          = std::distance(myBlockStart, myBlockEnd) / (pj * pi);

  // We also look for the index of the local rank in the current block.
  const int MyRank = myRank;
  myRankIndex      = std::distance(meshData.begin(),
                                   std::find_if(myBlockStart, myBlockEnd,
                                                [MyRank](const std::vector<GO>& vec) -> bool {
                                             return (vec[0] == MyRank) ? true : false;
                                           }));
  // We also construct a mapping of rank to rankIndex in the meshData vector,
  // this will allow us to access data quickly later on.
  for (int rankIndex = 0; rankIndex < numRanks; ++rankIndex) {
    rankIndices[meshData[rankIndex][0]] = rankIndex;
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    computeCoarseLocalLexicographicData() {
  Array<LO> rankOffset(3);
  for (int rank = 0; rank < numRanks; ++rank) {
    coarseMeshData[rank].resize(10);
    coarseMeshData[rank][0] = meshData[rank][0];
    coarseMeshData[rank][1] = meshData[rank][1];
    coarseMeshData[rank][2] = meshData[rank][2];
    for (int dim = 0; dim < 3; ++dim) {
      coarseMeshData[rank][3 + 2 * dim] = meshData[rank][3 + 2 * dim] / this->coarseRate[dim];
      if (meshData[rank][3 + 2 * dim] % this->coarseRate[dim] > 0) {
        ++coarseMeshData[rank][3 + 2 * dim];
      }
      coarseMeshData[rank][3 + 2 * dim + 1] = meshData[rank][3 + 2 * dim + 1] / this->coarseRate[dim];
      if (meshData[rank][3 + 2 * dim + 1] == this->gFineNodesPerDir[dim] - 1 &&
          meshData[rank][3 + 2 * dim + 1] % this->coarseRate[dim] > 0) {
        // this->endRate[dim] < this->coarseRate[dim]) {
        ++coarseMeshData[rank][3 + 2 * dim + 1];
      }
    }
    if (rank > 0) {
      coarseMeshData[rank][9] = coarseMeshData[rank - 1][9] + (coarseMeshData[rank - 1][8] - coarseMeshData[rank - 1][7] + 1) * (coarseMeshData[rank - 1][6] - coarseMeshData[rank - 1][5] + 1) * (coarseMeshData[rank - 1][4] - coarseMeshData[rank - 1][3] + 1);
    }
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
std::vector<std::vector<GlobalOrdinal> > LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseMeshData() const { return coarseMeshData; }

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getFineNodeGlobalTuple(const GO /* myGID */, GO& /* i */, GO& /* j */, GO& /* k */) const {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getFineNodeLocalTuple(const LO myLID, LO& i, LO& j, LO& k) const {
  LO tmp;
  k   = myLID / this->lNumFineNodes10;
  tmp = myLID % this->lNumFineNodes10;
  j   = tmp / this->lFineNodesPerDir[0];
  i   = tmp % this->lFineNodesPerDir[0];
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
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
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getFineNodeGID(const GO /* i */, const GO /* j */, const GO /* k */, GO& /* myGID */) const {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getFineNodeLID(const LO /* i */, const LO /* j */, const LO /* k */, LO& /* myLID */) const {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodeGlobalTuple(const GO /* myGID */, GO& /* i */, GO& /* j */, GO& /* k */) const {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodeLocalTuple(const LO myLID, LO& i, LO& j, LO& k) const {
  LO tmp;
  k   = myLID / this->lNumCoarseNodes10;
  tmp = myLID % this->lNumCoarseNodes10;
  j   = tmp / this->lCoarseNodesPerDir[0];
  i   = tmp % this->lCoarseNodesPerDir[0];
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodeGID(const GO /* i */, const GO /* j */, const GO /* k */, GO& /* myGID */) const {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodeLID(const LO /* i */, const LO /* j */, const LO /* k */, LO& /* myLID */) const {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodeGhostedLID(const LO i, const LO j, const LO k, LO& myLID) const {
  myLID = k * this->numGhostedNodes10 + j * this->ghostedNodesPerDir[0] + i;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
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
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getGhostedNodeFineLID(const LO /* i */, const LO /* j */, const LO /* k */, LO& /* myLID */) const {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getGhostedNodeCoarseLID(const LO /* i */, const LO /* j */, const LO /* k */, LO& /* myLID */) const {
}

}  // namespace MueLu

#endif /* MUELU_LOCALLEXICOGRAPHICINDEXMANAGER_DEF_HPP_ */
