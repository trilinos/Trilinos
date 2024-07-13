// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_UNCOUPLEDINDEXMANAGER_DEF_HPP_
#define MUELU_UNCOUPLEDINDEXMANAGER_DEF_HPP_

#include <Xpetra_MapFactory.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <MueLu_UncoupledIndexManager_decl.hpp>

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    UncoupledIndexManager(const RCP<const Teuchos::Comm<int> > comm, const bool coupled,
                          const int NumDimensions, const int interpolationOrder,
                          const int MyRank, const int NumRanks,
                          const Array<GO> GFineNodesPerDir, const Array<LO> LFineNodesPerDir,
                          const Array<LO> CoarseRate, const bool singleCoarsePoint)
  : IndexManager(comm, coupled, singleCoarsePoint, NumDimensions, interpolationOrder,
                 Array<GO>(3, -1), LFineNodesPerDir)
  , myRank(MyRank)
  , numRanks(NumRanks) {
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

  this->computeMeshParameters();
  this->gNumCoarseNodes10 = Teuchos::OrdinalTraits<GO>::invalid();
  this->gNumCoarseNodes   = Teuchos::OrdinalTraits<GO>::invalid();
}  // Constructor

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    computeGlobalCoarseParameters() {
  GO input[1] = {as<GO>(this->lNumCoarseNodes)}, output[1] = {0};
  Teuchos::reduceAll(*(this->comm_), Teuchos::REDUCE_SUM, 1, input, output);
  this->gNumCoarseNodes = output[0];
}  // computeGlobalCoarseParameters

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getGhostedNodesData(const RCP<const Map> /* fineMap */,
                        Array<LO>& ghostedNodeCoarseLIDs,
                        Array<int>& ghostedNodeCoarsePIDs,
                        Array<GO>& /* ghostedNodeCoarseGIDs */) const {
  // First we allocate memory for the outputs
  ghostedNodeCoarseLIDs.resize(this->getNumLocalGhostedNodes());
  ghostedNodeCoarsePIDs.resize(this->getNumLocalGhostedNodes());
  // In the uncoupled case the data required is trivial to provide!
  for (LO idx = 0; idx < this->getNumLocalGhostedNodes(); ++idx) {
    ghostedNodeCoarseLIDs[idx] = idx;
    ghostedNodeCoarsePIDs[idx] = myRank;
  }
}  // getGhostedNodesData

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodesData(const RCP<const Map> fineCoordinatesMap,
                       Array<GO>& coarseNodeCoarseGIDs,
                       Array<GO>& coarseNodeFineGIDs) const {
  // Allocate sufficient amount of storage in output arrays
  coarseNodeCoarseGIDs.resize(this->getNumLocalCoarseNodes());
  coarseNodeFineGIDs.resize(this->getNumLocalCoarseNodes());

  // Load all the GIDs on the fine mesh
  ArrayView<const GO> fineNodeGIDs = fineCoordinatesMap->getLocalElementList();

  // Extract the fine LIDs of the coarse nodes and store the corresponding GIDs
  LO fineLID;
  for (LO coarseLID = 0; coarseLID < this->getNumLocalCoarseNodes(); ++coarseLID) {
    Array<LO> coarseIndices(3), fineIndices(3);
    this->getCoarseNodeLocalTuple(coarseLID,
                                  coarseIndices[0],
                                  coarseIndices[1],
                                  coarseIndices[2]);
    for (int dim = 0; dim < 3; ++dim) {
      if (coarseIndices[dim] == this->lCoarseNodesPerDir[dim] - 1) {
        if (this->lCoarseNodesPerDir[dim] == 1) {
          fineIndices[dim] = 0;
        } else {
          fineIndices[dim] = this->lFineNodesPerDir[dim] - 1;
        }
      } else {
        fineIndices[dim] = coarseIndices[dim] * this->coarseRate[dim];
      }
    }

    fineLID                       = fineIndices[2] * this->lNumFineNodes10 + fineIndices[1] * this->lFineNodesPerDir[0] + fineIndices[0];
    coarseNodeFineGIDs[coarseLID] = fineNodeGIDs[fineLID];
  }
}  // getCoarseNodesData

template <class LocalOrdinal, class GlobalOrdinal, class Node>
std::vector<std::vector<GlobalOrdinal> > UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseMeshData() const {
  std::vector<std::vector<GO> > coarseMeshData;
  return coarseMeshData;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getFineNodeGlobalTuple(const GO /* myGID */, GO& /* i */, GO& /* j */, GO& /* k */) const {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getFineNodeLocalTuple(const LO myLID, LO& i, LO& j, LO& k) const {
  LO tmp;
  k   = myLID / this->lNumFineNodes10;
  tmp = myLID % this->lNumFineNodes10;
  j   = tmp / this->lFineNodesPerDir[0];
  i   = tmp % this->lFineNodesPerDir[0];
}  // getFineNodeLocalTuple

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getFineNodeGhostedTuple(const LO myLID, LO& i, LO& j, LO& k) const {
  LO tmp;
  k   = myLID / this->lNumFineNodes10;
  tmp = myLID % this->lNumFineNodes10;
  j   = tmp / this->lFineNodesPerDir[0];
  i   = tmp % this->lFineNodesPerDir[0];

  k += this->offsets[2];
  j += this->offsets[1];
  i += this->offsets[0];
}  // getFineNodeGhostedTuple

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getFineNodeGID(const GO /* i */, const GO /* j */, const GO /* k */, GO& /* myGID */) const {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getFineNodeLID(const LO /* i */, const LO /* j */, const LO /* k */, LO& /* myLID */) const {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodeGlobalTuple(const GO /* myGID */, GO& /* i */, GO& /* j */, GO& /* k */) const {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodeLocalTuple(const LO myLID, LO& i, LO& j, LO& k) const {
  LO tmp;
  k   = myLID / this->lNumCoarseNodes10;
  tmp = myLID % this->lNumCoarseNodes10;
  j   = tmp / this->lCoarseNodesPerDir[0];
  i   = tmp % this->lCoarseNodesPerDir[0];
}  // getCoarseNodeLocalTuple

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodeGID(const GO /* i */, const GO /* j */, const GO /* k */, GO& /* myGID */) const {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodeLID(const LO /* i */, const LO /* j */, const LO /* k */, LO& /* myLID */) const {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodeGhostedLID(const LO i, const LO j, const LO k, LO& myLID) const {
  myLID = k * this->numGhostedNodes10 + j * this->ghostedNodesPerDir[0] + i;
}  // getCoarseNodeGhostedLID

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getCoarseNodeFineLID(const LO /* i */, const LO /* j */, const LO /* k */, LO& /* myLID */) const {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getGhostedNodeFineLID(const LO /* i */, const LO /* j */, const LO /* k */, LO& /* myLID */) const {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    getGhostedNodeCoarseLID(const LO /* i */, const LO /* j */, const LO /* k */, LO& /* myLID */) const {
}

}  // namespace MueLu

#endif /* MUELU_UNCOUPLEDINDEXMANAGER_DEF_HPP_ */
