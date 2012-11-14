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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_AMALGAMATIONFACTORY_DEF_HPP
#define MUELU_AMALGAMATIONFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>

#include "MueLu_AmalgamationFactory.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AmalgamationFactory(RCP<const FactoryBase> AFact)
: AFact_(AFact)
  {
  }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {

  currentLevel.DeclareInput("A", AFact_.get(),    this); // sub-block from blocked A

}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const
{
  FactoryMonitor m(*this, "Build", currentLevel);

  RCP<Matrix> A = currentLevel.Get< RCP<Matrix> >("A", AFact_.get());

  LocalOrdinal  fullblocksize = 1;         // block dim for fixed size blocks
  GlobalOrdinal offset = 0;          // global offset of dof gids
  LocalOrdinal blockid = -1;         // block id in strided map
  LocalOrdinal nStridedOffset = 0;   // DOF offset for strided block id "blockid" (default = 0)
  LocalOrdinal stridedblocksize = fullblocksize; // size of strided block id "blockid" (default = fullblocksize, only if blockid!=-1 stridedblocksize <= fullblocksize)

  // 1) check for blocking/striding information
  if(A->IsView("stridedMaps") &&
      Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps")) != Teuchos::null) {
    Xpetra::viewLabel_t oldView = A->SwitchToView("stridedMaps"); // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap()) == Teuchos::null,Exceptions::BadCast,"MueLu::CoalesceFactory::Build: cast to strided row map failed.");
    fullblocksize = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getFixedBlockSize(); // TODO shorten code
    offset   = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getOffset();
    blockid  = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getStridedBlockId();
    if (blockid > -1) {
      std::vector<size_t> stridingInfo = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getStridingData();
      for(size_t j=0; j<Teuchos::as<size_t>(blockid); j++) {
        nStridedOffset += stridingInfo[j];
      }
      stridedblocksize = Teuchos::as<LocalOrdinal>(stridingInfo[blockid]);
    } else {
      stridedblocksize = fullblocksize;
    }
    oldView = A->SwitchToView(oldView);
    GetOStream(Statistics0, 0) << "AmalagamationFactory::Build():" << " found fullblocksize=" << fullblocksize << " and stridedblocksize=" << stridedblocksize << " from strided maps. offset=" << offset << std::endl;
    /*std::cout << "fullblocksize: " << fullblocksize << std::endl;
      std::cout << "offset: " << offset << std::endl;
      std::cout << "blockid: " << blockid << std::endl;
      std::cout << "nStridedOffset: " << nStridedOffset << std::endl;
      std::cout << "stridedblocksize: " << stridedblocksize << std::endl;*/
  } else GetOStream(Warnings0, 0) << "AmalagamationFactory::Build(): no striding information available. Use blockdim=1 with offset=0" << std::endl;
  // TODO: maybe no striding information on coarser levels -> misuse nullspace vector?

  // 2) prepare maps for amalgamated graph of A and
  //    setup unamalgamation information

  RCP<std::vector<GlobalOrdinal> > gNodeIds; // contains global node ids on current proc
  gNodeIds = Teuchos::rcp(new std::vector<GlobalOrdinal>);
  gNodeIds->empty();

  // in nodegid2dofgids_ for each node on the current proc a vector of
  // the corresponding DOFs gids is stored.
  // The map contains all nodes the current proc has connections to (including
  // nodes that are stored on other procs when there are off-diagonal entries in A)
  nodegid2dofgids_ = Teuchos::rcp(new std::map<GlobalOrdinal,std::vector<GlobalOrdinal> >);

  // extract information from overlapping column map of A
  Teuchos::RCP<const Map> colMap = A->getColMap();
  GlobalOrdinal cnt_amalRows = 0; // counts number of nodes (rows in amalgamated matrix) on current proc
  LocalOrdinal nColEle = Teuchos::as<LocalOrdinal>(A->getColMap()->getNodeNumElements());
  for(LocalOrdinal i=0; i<nColEle;i++) {
    // get global DOF id
    GlobalOrdinal gDofId = colMap->getGlobalElement(i);

    // translate DOFGid to node id
    GlobalOrdinal gNodeId = DOFGid2NodeId(gDofId, A, fullblocksize, offset);

    // gblockid -> gDofId/lDofId
    if(nodegid2dofgids_->count(gNodeId) == 0) {

      // current column DOF gDofId belongs to a node that has not been added
      // to nodeid2dofgids_ yet. Do it now and add ALL DOFs of node gNodeId to
      // unamalgamation information.
      // Note: we use offset and fullblocksize, ie. information from strided maps indirectly
      std::vector<GlobalOrdinal> DOFs;

      DOFs.reserve(stridedblocksize);
      for(LocalOrdinal k=0; k<stridedblocksize; k++) {
        DOFs.push_back(offset + gNodeId*fullblocksize + nStridedOffset + k);
      }

      (*nodegid2dofgids_)[gNodeId] = DOFs;

      if(A->getRowMap()->isNodeGlobalElement(gDofId)) {
        gNodeIds->push_back(gNodeId);
        cnt_amalRows++; // new local block row in amalgamated matrix graph
      }
    }
  }

  // store (un)amalgamation information on current level
  RCP<AmalgamationInfo> amalgamationData = rcp(new AmalgamationInfo());
  amalgamationData->SetAmalgamationParams(nodegid2dofgids_);
  amalgamationData->SetNodeGIDVector(gNodeIds);
  amalgamationData->SetNumberOfNodes(cnt_amalRows);
  currentLevel.Set("UnAmalgamationInfo", amalgamationData, this);
}

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const GlobalOrdinal AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DOFGid2NodeId(GlobalOrdinal gid, const RCP<Matrix>& A, LocalOrdinal blockSize, const GlobalOrdinal offset) {
  GlobalOrdinal globalblockid = ((GlobalOrdinal) gid - offset) / blockSize;
  return globalblockid;
}

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ComputeUnamalgamatedAggregateSizes(const Aggregates & aggregates, const AmalgamationInfo & amalgInfo, Teuchos::ArrayRCP<LocalOrdinal> & aggSizes) {
  // we expect the aggSizes array to be initialized as follows
  // aggSizes = Teuchos::ArrayRCP<LO>(nAggregates_,0);
  // furthermore we suppose the (un)amalgamation info to be set (even for 1 dof per node examples)

  int myPid = aggregates.GetMap()->getComm()->getRank();
  Teuchos::ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);
  Teuchos::ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  LO size = procWinner.size();

  //for (LO i = 0; i< aggregates.GetNumAggregates(); ++i) aggSizes[i] = 0;
  for (LO lnode = 0; lnode < size; ++lnode) {
    LO myAgg = vertex2AggId[lnode];
    if (procWinner[lnode] == myPid) {
      GO gnodeid = aggregates.GetMap()->getGlobalElement(lnode);

      std::vector<GO> gDofIds = (*(amalgInfo.GetGlobalAmalgamationParams()))[gnodeid];
      aggSizes[myAgg] += Teuchos::as<LO>(gDofIds.size());
    }
  }
}

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::UnamalgamateAggregates(const Aggregates& aggregates, const AmalgamationInfo& amalgInfo, const Teuchos::ArrayRCP<LocalOrdinal> & aggSizes, Teuchos::ArrayRCP<Teuchos::ArrayRCP<GlobalOrdinal> > & aggToRowMap) {
  int myPid = aggregates.GetMap()->getComm()->getRank();
  Teuchos::ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);
  Teuchos::ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  LO size = procWinner.size();

  // initialize array aggToRowMap with empty arrays for each aggregate (with correct aggSize)
  LO t = 0;
  for (typename ArrayRCP<ArrayRCP<GO> >::iterator a2r = aggToRowMap.begin(); a2r!=aggToRowMap.end(); ++a2r) {
    *a2r = ArrayRCP<GO>(aggSizes[t++]);
  }

  // count, how many dofs have been recorded for each aggregate
  ArrayRCP<LO> numDofs(aggregates.GetNumAggregates(),0); // empty array with number of Dofs for each aggregate

  for (LO lnode = 0; lnode < size; ++lnode) {
    LO myAgg = vertex2AggId[lnode];
    if (procWinner[lnode] == myPid) {
      GO gnodeid = aggregates.GetMap()->getGlobalElement(lnode);
      std::vector<GO> gDofIds = (*(amalgInfo.GetGlobalAmalgamationParams()))[gnodeid];
      LO gDofIds_size = Teuchos::as<LO>(gDofIds.size());
      for (LO gDofId=0; gDofId < gDofIds_size; gDofId++) {
        aggToRowMap[ myAgg ][ numDofs[myAgg] ] = gDofIds[gDofId]; // fill aggToRowMap structure
        ++(numDofs[myAgg]);
      }
    }
  }
  // todo plausibility check: entry numDofs[k] == aggToRowMap[k].size()

}

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ComputeUnamalgamatedImportDofMap(const Aggregates& aggregates, const AmalgamationInfo& amalgInfo, const Teuchos::ArrayRCP<LocalOrdinal> & aggSizes) {
  Teuchos::RCP<const Map> nodeMap = aggregates.GetMap(); //aggregates.GetVertex2AggId();

  Teuchos::RCP<std::vector<GO> > myDofGids = Teuchos::rcp(new std::vector<GO>);
  LO nodeElements = Teuchos::as<LO>(nodeMap->getNodeNumElements());
  for(LO n = 0; n<nodeElements; n++) {
    GO gnodeid = (GO) nodeMap->getGlobalElement(n);
    std::vector<GO> gDofIds = (*(amalgInfo.GetGlobalAmalgamationParams()))[gnodeid];
    for(typename std::vector<GO>::iterator gDofIdsIt = gDofIds.begin(); gDofIdsIt != gDofIds.end(); gDofIdsIt++) {
      myDofGids->push_back(*gDofIdsIt);
    }
  }

  Teuchos::ArrayRCP<GO> arr_myDofGids = Teuchos::arcp( myDofGids );
  Teuchos::RCP<Map> importDofMap = MapFactory::Build(aggregates.GetMap()->lib(),
      Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), arr_myDofGids(),
      aggregates.GetMap()->getIndexBase(), aggregates.GetMap()->getComm());
  return importDofMap;
}

} //namespace MueLu

#endif /* MUELU_SUBBLOCKUNAMALGAMATIONFACTORY_DEF_HPP */

