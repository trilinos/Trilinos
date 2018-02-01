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
//                    Jonathan Hu        (jhu@sandia.gov)
//                    Ray Tuminaro       (rstumin@sandia.gov)
//                    Luc Berger-Vergoat (lberge@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_LOCALLEXICOGRAPHICINDEXMANAGER_DEF_HPP_
#define MUELU_LOCALLEXICOGRAPHICINDEXMANAGER_DEF_HPP_

#include <MueLu_LocalLexicographicIndexManager_decl.hpp>

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  LocalLexicographicIndexManager(const int NumDimensions, const int MyRank, const int NumRanks,
                                 const Array<GO> GFineNodesPerDir, const Array<LO> LFineNodesPerDir,
                                 const Array<LO> CoarseRate, const Array<GO> MeshData) :
    IndexManager(NumDimensions, GFineNodesPerDir, LFineNodesPerDir), myRank(MyRank),
    numRanks(NumRanks) {

    // Load coarse rate, being careful about formating
    for(int dim = 0; dim < 3; ++dim) {
      if(dim < this->numDimensions) {
        if(CoarseRate.size() == 1) {
          this->coarseRate[dim] = CoarseRate[0];
        } else if(CoarseRate.size() == this->numDimensions) {
          this->coarseRate[dim] = CoarseRate[dim];
        }
      } else {
        this->coarseRate[dim] = 1;
      }
    }

    // Load meshData for local lexicographic case
    meshData.resize(numRanks);
    for(int rank = 0; rank < numRanks; ++rank) {
      meshData[rank].resize(10);
      for(int entry = 0; entry < 10; ++entry) {
        meshData[rank][entry] = MeshData[10*rank + entry];
      }
    }

    myBlock = meshData[myRank][2];
    sortLocalLexicographicData();

    // Start simple parameter calculation
    myRankIndex = rankIndices[myRank];
    for(int dim = 0; dim < 3; ++dim) {
      this->startIndices[dim]     = meshData[myRankIndex][2*dim + 3];
      this->startIndices[dim + 3] = meshData[myRankIndex][2*dim + 4];
    }

    this->computeMeshParameters();

  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  printMeshInfo() {

  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  sortLocalLexicographicData() {

    std::sort(meshData.begin(), meshData.end(),
              [](const std::vector<GO>& a, const std::vector<GO>& b)->bool {
                // The below function sorts ranks by blockID, kmin, jmin and imin
                if(a[2] < b[2]) {
                  return true;
                } else if(a[2] == b[2]) {
                  if(a[7] < b[7]) {
                    return true;
                  } else if(a[7] == b[7]) {
                    if(a[5] < b[5]) {
                      return true;
                    } else if(a[5] == b[5]) {
                      if(a[3] < b[3]) {return true;}
                    }
                  }
                }
                return false;
              });

    numBlocks = meshData[numRanks - 1][2] + 1;
    // Find the range of the current block
    myBlockStart = std::lower_bound(meshData.begin(), meshData.end(), myBlock - 1,
                                    [] (const std::vector<GO>& vec, const GO val)->bool {
                                      return (vec[2] < val) ? true : false;
                                    });
    myBlockEnd = std::upper_bound(meshData.begin(), meshData.end(), myBlock,
                                  [] (const GO val, const std::vector<GO>& vec)->bool {
                                    return (val < vec[2]) ? true : false;
                                  });
    // Assuming that i,j,k and ranges are split in pi, pj and pk processors
    // we search for these numbers as they will allow us to find quickly the PID of processors
    // owning ghost nodes.
    auto myKEnd = std::upper_bound(myBlockStart, myBlockEnd, (*myBlockStart)[3],
                                   [] (const GO val, const std::vector<GO>& vec)->bool {
                                     return (val < vec[7]) ? true : false;
                                   });
    auto myJEnd = std::upper_bound(myBlockStart, myKEnd, (*myBlockStart)[3],
                                   [] (const GO val, const std::vector<GO>& vec)->bool {
                                     return (val < vec[5]) ? true : false;
                                   });
    pi = std::distance(myBlockStart, myJEnd);
    pj = std::distance(myBlockStart, myKEnd) / pi;
    pk = std::distance(myBlockStart, myBlockEnd) / (pj*pi);

    // We also look for the index of the local rank in the current block.
    const int MyRank = myRank;
    myRankIndex = std::distance(meshData.begin(),
                                std::find_if(myBlockStart, myBlockEnd,
                                             [MyRank] (const std::vector<GO>& vec)->bool {
                                               return (vec[0] == MyRank) ? true : false;
                                             })
                                );
    // We also construct a mapping of rank to rankIndex in the meshData vector,
    // this will allow us to access data quickly later on.
    rankIndices.resize(numRanks);
    for(int rankIndex = 0; rankIndex < numRanks; ++rankIndex) {
      rankIndices[meshData[rankIndex][0]] = rankIndex;
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  computeCoarseLocalLexicographicData() {
    coarseMeshData.resize(numRanks);
    Array<LO> rankOffset(3);
    for(int rank = 0; rank < numRanks; ++rank) {
      coarseMeshData[rank].resize(10);
      coarseMeshData[rank][0] = meshData[rank][0];
      coarseMeshData[rank][1] = meshData[rank][1];
      coarseMeshData[rank][2] = meshData[rank][2];
      for(int dim = 0; dim < 3; ++dim) {
        coarseMeshData[rank][3 + 2*dim] = meshData[rank][3 + 2*dim] / this->coarseRate[dim];
        if(meshData[rank][3 + 2*dim] % this->coarseRate[dim] > 0) {
          ++coarseMeshData[rank][3 + 2*dim];
        }
        coarseMeshData[rank][3 + 2*dim + 1] = meshData[rank][3 + 2*dim + 1] / this->coarseRate[dim];
        if(meshData[rank][3 + 2*dim + 1] == this->gFineNodesPerDir[dim] - 1 &&
           this->endRate[dim] < this->coarseRate[dim]) {
          ++coarseMeshData[rank][3 + 2*dim + 1];
        }
      }
      if(rank > 0) {
        coarseMeshData[rank][9] = coarseMeshData[rank - 1][9]
          + (coarseMeshData[rank - 1][8] - coarseMeshData[rank - 1][7] + 1)
          * (coarseMeshData[rank - 1][6] - coarseMeshData[rank - 1][5] + 1)
          * (coarseMeshData[rank - 1][4] - coarseMeshData[rank - 1][3] + 1);
      }
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getFineNodeGlobalTuple(const GO myGID, GO& i, GO& j, GO& k) const {
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getFineNodeLocalTuple(const LO myLID, LO& i, LO& j, LO& k) const {
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getFineNodeGhostedTuple(const LO myLID, LO& i, LO& j, LO& k) const {
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getFineNodeGID(const GO i, const GO j, const GO k, GO& myGID) const {
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getFineNodeLID(const LO i, const LO j, const LO k, LO& myLID) const {
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getCoarseNodeGlobalTuple(const GO myGID, GO& i, GO& j, GO& k) const {
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getCoarseNodeLocalTuple(const LO myLID, LO& i, LO& j, LO& k) const {
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getCoarseNodeGID(const GO i, const GO j, const GO k, GO& myGID) const {
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getCoarseNodeLID(const LO i, const LO j, const LO k, LO& myLID) const {
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getCoarseNodeGhostedLID(const LO i, const LO j, const LO k, LO& myLID) const {
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getCoarseNodeFineLID(const LO i, const LO j, const LO k, LO& myLID) const {
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getGhostedNodeFineLID(const LO i, const LO j, const LO k, LO& myLID) const {
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getGhostedNodeCoarseLID(const LO i, const LO j, const LO k, LO& myLID) const {
  }

} //namespace MueLu

#endif /* MUELU_LOCALLEXICOGRAPHICINDEXMANAGER_DEF_HPP_ */
