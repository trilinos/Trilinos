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
                                  const GO MinGlobalIndex) :
    IndexManager(comm, coupled, NumDimensions, interpolationOrder, GFineNodesPerDir, LFineNodesPerDir) {

    // Load coarse rate, being careful about formating.
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

    {
      GO tmp = 0;
      this->startIndices[2]= MinGlobalIndex / (this->gFineNodesPerDir[1]*this->gFineNodesPerDir[0]);
      tmp                  = MinGlobalIndex % (this->gFineNodesPerDir[1]*this->gFineNodesPerDir[0]);
      this->startIndices[1]= tmp / this->gFineNodesPerDir[0];
      this->startIndices[0]= tmp % this->gFineNodesPerDir[0];

      for(int dim = 0; dim < 3; ++dim) {
        this->startIndices[dim + 3] = this->startIndices[dim] + this->lFineNodesPerDir[dim] - 1;
      }
    }

    this->computeMeshParameters();

  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getGhostedNodesData(const RCP<const Map> fineMap, RCP<const Map> coarseMap,
                     Array<LO>& ghostedNodeCoarseLIDs, Array<int>& ghostedNodeCoarsePIDs) const {

    ghostedNodeCoarseLIDs.resize(this->getNumLocalGhostedNodes());
    ghostedNodeCoarsePIDs.resize(this->getNumLocalGhostedNodes());

    // Find the GIDs, LIDs and PIDs of the coarse points on the fine mesh and coarse
    // mesh as this data will be used to fill vertex2AggId and procWinner vectors.
    Array<GO> lCoarseNodeCoarseGIDs(this->lNumCoarseNodes),
      lCoarseNodeFineGIDs(this->lNumCoarseNodes);
    Array<GO> ghostedNodeCoarseGIDs(this->numGhostedNodes),
      ghostedCoarseNodeFineGIDs(this->numGhostedNodes);
    Array<LO> ghostedCoarseNodeCoarseIndices(3), ghostedCoarseNodeFineIndices(3), ijk(3);
    LO currentIndex = -1, coarseNodeFineLID = -1, currentCoarseIndex = -1;
    for(ijk[2] = 0; ijk[2] < this->ghostedNodesPerDir[2]; ++ijk[2]) {
      for(ijk[1] = 0; ijk[1] < this->ghostedNodesPerDir[1]; ++ijk[1]) {
        for(ijk[0] = 0; ijk[0] < this->ghostedNodesPerDir[0]; ++ijk[0]) {
          currentIndex = ijk[2]*this->numGhostedNodes10 + ijk[1]*this->ghostedNodesPerDir[0] + ijk[0];
          ghostedCoarseNodeCoarseIndices[0] = this->startGhostedCoarseNode[0] + ijk[0];
          ghostedCoarseNodeCoarseIndices[1] = this->startGhostedCoarseNode[1] + ijk[1];
          ghostedCoarseNodeCoarseIndices[2] = this->startGhostedCoarseNode[2] + ijk[2];
          GO myCoarseGID = ghostedCoarseNodeCoarseIndices[0]
            + ghostedCoarseNodeCoarseIndices[1]*this->gCoarseNodesPerDir[0]
            + ghostedCoarseNodeCoarseIndices[2]*this->gNumCoarseNodes10;
          ghostedNodeCoarseGIDs[currentIndex] = myCoarseGID;
          GO myGID = 0, factor[3] = {};
          factor[2] = this->gNumFineNodes10;
          factor[1] = this->gFineNodesPerDir[0];
          factor[0] = 1;
          for(int dim = 0; dim < 3; ++dim) {
            if(dim < this->numDimensions) {
              if(this->startIndices[dim] - this->offsets[dim] + ijk[dim]*this->coarseRate[dim]
                 < this->gFineNodesPerDir[dim] - 1) {
                myGID += (this->startIndices[dim] - this->offsets[dim]
                          + ijk[dim]*this->coarseRate[dim])*factor[dim];
              } else {
                myGID += (this->startIndices[dim] - this->offsets[dim] + (ijk[dim] - 1)
                          *this->coarseRate[dim] + this->endRate[dim])*factor[dim];
              }
            }
          }
          // lbv 02-08-2018:
          // This check is simplistic and should be replaced by a condition that checks
          // if the local tuple of the current index is wihin the range of local nodes
          // or not in the range of ghosted nodes.
          if((!this->ghostInterface[0] || ijk[0] != 0) &&
             (!this->ghostInterface[2] || ijk[1] != 0) &&
             (!this->ghostInterface[4] || ijk[2] != 0) &&
             (!this->ghostInterface[1] || ijk[0] != this->ghostedNodesPerDir[0] - 1) &&
             (!this->ghostInterface[3] || ijk[1] != this->ghostedNodesPerDir[1] - 1) &&
             (!this->ghostInterface[5] || ijk[2] != this->ghostedNodesPerDir[2] - 1)) {

            // this->getGhostedNodeFineLID(ijk[0], ijk[1], ijk[2], coarseNodeFineLID);
            if(this->interpolationOrder_ == 0) {
              currentCoarseIndex = 0;
              if(this->ghostInterface[4]) {
                currentCoarseIndex += (ijk[2] - 1)*this->lNumCoarseNodes10;
              } else {
                currentCoarseIndex += ijk[2]*this->lNumCoarseNodes10;
              }
              if(this->ghostInterface[2]) {
                currentCoarseIndex += (ijk[1] - 1)*this->getLocalCoarseNodesInDir(0);
              } else {
                currentCoarseIndex += ijk[1]*this->getLocalCoarseNodesInDir(0);
              }
              if(this->ghostInterface[0]) {
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

    coarseMap = Xpetra::MapFactory<LO,GO,NO>::Build (fineMap->lib(),
                                                     this->gNumCoarseNodes,
                                                     lCoarseNodeCoarseGIDs(),
                                                     fineMap->getIndexBase(),
                                                     fineMap->getComm());

    // Array<int> ghostedCoarseNodeCoarsePIDs(this->numGhostedNodes);
    // Array<LO>  ghostedCoarseNodeCoarseLIDs(this->numGhostedNodes);
    coarseMap->getRemoteIndexList(ghostedNodeCoarseGIDs(),
                                  ghostedNodeCoarsePIDs(),
                                  ghostedNodeCoarseLIDs());

  } // End getGhostedMeshData

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
    myGID = k*this->gNumFineNodes10 + j*this->gFineNodesPerDir[0] + i;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getFineNodeLID(const LO i, const LO j, const LO k, LO& myLID) const {
    myLID = k*this->lNumFineNodes10 + j*this->lFineNodesPerDir[0] + i;
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
    myGID = k*this->gNumCoarseNodes10 + j*this->gCoarseNodesPerDir[0] + i;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getCoarseNodeLID(const LO i, const LO j, const LO k, LO& myLID) const {
    myLID = k*this->lNumCoarseNodes10 + j*this->lCoarseNodesPerDir[0] + i;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getCoarseNodeGhostedLID(const LO i, const LO j, const LO k, LO& myLID) const {
    myLID = k*this->numGhostedNodes10 + j*this->ghostedNodesPerDir[0] + i;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getCoarseNodeFineLID(const LO i, const LO j, const LO k, LO& myLID) const {
    myLID = 0;
    if(k*this->coarseRate[2] < this->lFineNodesPerDir[2]) {
      myLID += k*this->coarseRate[2]*this->lNumCoarseNodes10;
    } else {
      myLID += (this->lFineNodesPerDir[2] - 1)*this->lNumCoarseNodes10;
    }

    if(j*this->coarseRate[1] < this->lFineNodesPerDir[1]) {
      myLID += j*this->coarseRate[1]*this->lFineNodesPerDir[0];
    } else {
      myLID += (this->lFineNodesPerDir[1] - 1)*this->lFineNodesPerDir[1];
    }

    if(i*this->coarseRate[0] < this->lFineNodesPerDir[0]) {
      myLID += i*this->coarseRate[0];
    } else {
      myLID += this->lFineNodesPerDir[0] - 1;
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  getGhostedNodeFineLID(const LO i, const LO j, const LO k, LO& myLID) const {
    LO itmp = i - (this->offsets[0] > 0 ? 1 : 0);
    LO jtmp = j - (this->offsets[1] > 0 ? 1 : 0);
    LO ktmp = k - (this->offsets[2] > 0 ? 1 : 0);
    myLID = 0;
    if(ktmp*this->coarseRate[2] < this->lFineNodesPerDir[2]) {
      myLID += ktmp*this->coarseRate[2]*this->lNumCoarseNodes10;
    } else {
      myLID += (this->lFineNodesPerDir[2] - 1)*this->lNumCoarseNodes10;
    }

    if(jtmp*this->coarseRate[1] < this->lFineNodesPerDir[1]) {
      myLID += jtmp*this->coarseRate[1]*this->lFineNodesPerDir[0];
    } else {
      myLID += (this->lFineNodesPerDir[1] - 1)*this->lFineNodesPerDir[1];
    }

    if(itmp*this->coarseRate[0] < this->lFineNodesPerDir[0]) {
      myLID += itmp*this->coarseRate[0];
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
    myLID = ktmp*this->lNumCoarseNodes10 + jtmp*this->lCoarseNodesPerDir[0] + itmp;
  }

} //namespace MueLu

#endif /* MUELU_GLOBALLEXICOGRAPHICINDEXMANAGER_DEF_HPP_ */
