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
/*
 * MueLu_AmalgamationInfo_def.hpp
 *
 *  Created on: Mar 28, 2012
 *      Author: wiesner
 */

#ifndef MUELU_AMALGAMATIONINFO_KOKKOS_DEF_HPP_
#define MUELU_AMALGAMATIONINFO_KOKKOS_DEF_HPP_

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Vector.hpp>

#include "MueLu_Exceptions.hpp"

#include "MueLu_AmalgamationInfo_kokkos_decl.hpp"
#include "MueLu_Aggregates_kokkos.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void AmalgamationInfo_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  UnamalgamateAggregates(const Aggregates_kokkos& aggregates,
                         Teuchos::ArrayRCP<LocalOrdinal>& aggStart,
                         Teuchos::ArrayRCP<GlobalOrdinal>& aggToRowMap) const {

    int myPid = aggregates.GetMap()->getComm()->getRank();
    Teuchos::ArrayView<const GO> nodeGlobalElts = aggregates.GetMap()->getLocalElementList();
    Teuchos::ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);
    Teuchos::ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    LO size = procWinner.size();
    GO numAggregates = aggregates.GetNumAggregates();

    std::vector<LO> sizes(numAggregates);
    if (stridedblocksize_ == 1) {
      for (LO lnode = 0; lnode < size; ++lnode) {
        LO myAgg = vertex2AggId[lnode];
        if (procWinner[lnode] == myPid)
          sizes[myAgg] += 1;
      }
    } else {
      for (LO lnode = 0; lnode < size; ++lnode) {
        LO myAgg = vertex2AggId[lnode];
        if (procWinner[lnode] == myPid) {
          GO gnodeid = nodeGlobalElts[lnode];
          for (LocalOrdinal k = 0; k < stridedblocksize_; k++) {
            GlobalOrdinal gDofIndex = ComputeGlobalDOF(gnodeid,k);
            if (columnMap_->isNodeGlobalElement(gDofIndex))
              sizes[myAgg] += 1;
          }
        }
      }
    }
    aggStart = ArrayRCP<LO>(numAggregates+1,0);
    aggStart[0]=0;
    for (GO i=0; i<numAggregates; ++i) {
      aggStart[i+1] = aggStart[i] + sizes[i];
    }
    aggToRowMap = ArrayRCP<GO>(aggStart[numAggregates],0);

    // count, how many dofs have been recorded for each aggregate so far
    Array<LO> numDofs(numAggregates, 0); // empty array with number of Dofs for each aggregate

    if (stridedblocksize_ == 1) {
      for (LO lnode = 0; lnode < size; ++lnode) {
        LO myAgg = vertex2AggId[lnode];
        if (procWinner[lnode] == myPid) {
          aggToRowMap[ aggStart[myAgg] + numDofs[myAgg] ] = ComputeGlobalDOF(nodeGlobalElts[lnode]);
          ++(numDofs[myAgg]);
        }
      }
    } else {
      for (LO lnode = 0; lnode < size; ++lnode) {
        LO myAgg = vertex2AggId[lnode];

        if (procWinner[lnode] == myPid) {
          GO gnodeid = nodeGlobalElts[lnode];
          for (LocalOrdinal k = 0; k < stridedblocksize_; k++) {
            GlobalOrdinal gDofIndex = ComputeGlobalDOF(gnodeid,k);
            if (columnMap_->isNodeGlobalElement(gDofIndex)) {
              aggToRowMap[ aggStart[myAgg] + numDofs[myAgg] ] = gDofIndex;
              ++(numDofs[myAgg]);
            }
          }
        }
      }
    }
    // todo plausibility check: entry numDofs[k] == aggToRowMap[k].size()

  } //UnamalgamateAggregates

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void AmalgamationInfo_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  UnamalgamateAggregatesLO(const Aggregates_kokkos& aggregates,
                           Teuchos::ArrayRCP<LO>& aggStart,
                           Teuchos::ArrayRCP<LO>& aggToRowMap) const {

    int myPid = aggregates.GetMap()->getComm()->getRank();
    Teuchos::ArrayView<const GO> nodeGlobalElts = aggregates.GetMap()->getLocalElementList();

    Teuchos::ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()  ->getDataNonConst(0);
    Teuchos::ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    const GO numAggregates             = aggregates.GetNumAggregates();


    // FIXME: Do we need to compute size here? Or can we use existing?
    LO size = procWinner.size();

    std::vector<LO> sizes(numAggregates);
    if (stridedblocksize_ == 1) {
      for (LO lnode = 0; lnode < size; lnode++)
        if (procWinner[lnode] == myPid)
          sizes[vertex2AggId[lnode]]++;
    } else {
      for (LO lnode = 0; lnode < size; lnode++)
        if (procWinner[lnode] == myPid) {
          GO nodeGID = nodeGlobalElts[lnode];

          for (LO k = 0; k < stridedblocksize_; k++) {
            GO GID = ComputeGlobalDOF(nodeGID, k);
            if (columnMap_->isNodeGlobalElement(GID))
              sizes[vertex2AggId[lnode]]++;
          }
        }
    }

    aggStart = ArrayRCP<LO>(numAggregates+1); // FIXME: useless initialization with zeros
    aggStart[0] = 0;
    for (GO i = 0; i < numAggregates; i++)
      aggStart[i+1] = aggStart[i] + sizes[i];

    aggToRowMap = ArrayRCP<LO>(aggStart[numAggregates], 0);

    // count, how many dofs have been recorded for each aggregate so far
    Array<LO> numDofs(numAggregates, 0); // empty array with number of DOFs for each aggregate
    if (stridedblocksize_ == 1) {
      for (LO lnode = 0; lnode < size; ++lnode)
        if (procWinner[lnode] == myPid) {
          LO myAgg = vertex2AggId[lnode];
          aggToRowMap[aggStart[myAgg] + numDofs[myAgg]] = lnode;
          numDofs[myAgg]++;
        }
    } else {
      for (LO lnode = 0; lnode < size; ++lnode)
        if (procWinner[lnode] == myPid) {
          LO myAgg = vertex2AggId[lnode];
          GO nodeGID = nodeGlobalElts[lnode];

          for (LO k = 0; k < stridedblocksize_; k++) {
            GO GID = ComputeGlobalDOF(nodeGID, k);
            if (columnMap_->isNodeGlobalElement(GID)) {
              aggToRowMap[aggStart[myAgg] + numDofs[myAgg]] = lnode*stridedblocksize_ + k;
              numDofs[myAgg]++;
            }
          }
        }
    }
    // todo plausibility check: entry numDofs[k] == aggToRowMap[k].size()

  } //UnamalgamateAggregates

  /////////////////////////////////////////////////////////////////////////////

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > AmalgamationInfo_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  ComputeUnamalgamatedImportDofMap(const Aggregates_kokkos& aggregates) const {

    Teuchos::RCP<const Map> nodeMap = aggregates.GetMap();

    Teuchos::RCP<std::vector<GO> > myDofGids = Teuchos::rcp(new std::vector<GO>);
    Teuchos::ArrayView<const GO> gEltList = nodeMap->getLocalElementList();
    LO nodeElements = Teuchos::as<LO>(nodeMap->getLocalNumElements());
    if (stridedblocksize_ == 1) {
      for (LO n = 0; n<nodeElements; n++) {
        GlobalOrdinal gDofIndex = ComputeGlobalDOF(gEltList[n]);
        myDofGids->push_back(gDofIndex);
      }
    } else {
      for (LO n = 0; n<nodeElements; n++) {
        for (LocalOrdinal k = 0; k < stridedblocksize_; k++) {
          GlobalOrdinal gDofIndex = ComputeGlobalDOF(gEltList[n],k);
          if (columnMap_->isNodeGlobalElement(gDofIndex))
            myDofGids->push_back(gDofIndex);
        }
      }
    }

    Teuchos::ArrayRCP<GO> arr_myDofGids = Teuchos::arcp( myDofGids );
    Teuchos::RCP<Map> importDofMap = MapFactory::Build(aggregates.GetMap()->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), arr_myDofGids(), aggregates.GetMap()->getIndexBase(), aggregates.GetMap()->getComm());
    return importDofMap;
  }

  /////////////////////////////////////////////////////////////////////////////

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  GlobalOrdinal AmalgamationInfo_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  ComputeGlobalDOF(GlobalOrdinal const &gNodeID, LocalOrdinal const &k) const {
    // here, the assumption is, that the node map has the same indexBase as the dof map
    //                            this is the node map index base                    this is the dof map index base
    GlobalOrdinal gDofIndex = offset_ + (gNodeID-indexBase_)*fullblocksize_ + nStridedOffset_ + k + indexBase_;
    return gDofIndex;
  }

} //namespace


#endif /* MUELU_AMALGAMATIONINFO_KOKKOS_DEF_HPP_ */
