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

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
  GlobalLexicographicIndexManager(const int NumDimensions, const Array<GO> GFineNodesPerDir,
                                  const Array<LO> LFineNodesPerDir, const Array<LO> CoarseRate) :
    IndexManager(NumDimensions, GFineNodesPerDir, LFineNodesPerDir) {

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
    myLID = k*this->numGhostedCoarseNodes10 + j*this->ghostedCoarseNodesPerDir[0] + i;
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
