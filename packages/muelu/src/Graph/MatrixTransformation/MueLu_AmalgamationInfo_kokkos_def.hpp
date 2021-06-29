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

#ifdef HAVE_MUELU_KOKKOS_REFACTOR
#include "MueLu_AmalgamationInfo_kokkos_decl.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void AmalgamationInfo_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  UnamalgamateAggregates(const Aggregates_kokkos& aggregates,
                         Kokkos::View<LO*, memory_space>& aggStart,
                         Kokkos::View<GO*, memory_space>& aggToRowMap) const {

    int myPid = aggregates.GetMap()->getComm()->getRank();
    auto nodeLocalMap = aggregates.GetMap()->getLocalMap();
    auto procWinner   = aggregates.GetProcWinner()  ->template getLocalView<memory_space>();
    auto vertex2AggId = aggregates.GetVertex2AggId()->template getLocalView<memory_space>();
    LO size = static_cast<LO>(procWinner.extent(0));
    GO numAggregates = aggregates.GetNumAggregates();

    Kokkos::View<LO*, memory_space> sizes("sizes", numAggregates);
    if (stridedblocksize_ == 1) {
      Kokkos::parallel_for("compute aggregates sizes",
                           Kokkos::RangePolicy<execution_space>(0, size),
                           KOKKOS_LAMBDA(const LO& lnode) {
                             LO myAgg = vertex2AggId(lnode, 0);
                             if (procWinner(lnode, 0) == myPid) {
                               Kokkos::atomic_increment(&sizes(myAgg));
                             }
                           });
    } else {
      Kokkos::parallel_for("compute aggregates sizes",
                           Kokkos::RangePolicy<execution_space>(0, size),
                           KOKKOS_LAMBDA(const LO& lnode) {
                             LO myAgg = vertex2AggId(lnode, 0);
                             if (procWinner(lnode, 0) == myPid) {
                               // GO gnodeid = nodeLocalMap.getGlobalElement(lnode);
                               for (LO k = 0; k < stridedblocksize_; ++k) {
                                 // // Note: LBV on 10/11/19, removing the test
                                 // // below as there is no way to refactor it
                                 // // correctly at the moment since the method
                                 // // is not inlined for GPU kernels...
                                 // GO gDofIndex = ComputeGlobalDOF(gnodeid, k);
                                 // if (columnMap_->isNodeGlobalElement(gDofIndex)) {
                                 Kokkos::atomic_increment(&sizes(myAgg));
                                 // }
                               }
                             }
                           });
    }

    aggStart = Kokkos::View<LO*, memory_space>("aggStart", numAggregates+1);
    Kokkos::parallel_scan("compute aggregate starts",
                          Kokkos::RangePolicy<execution_space>(0, numAggregates+1),
                          KOKKOS_LAMBDA(const LO aggIdx, LO& update, const bool final) {
                            if(final) {aggStart(aggIdx) = update;}
                            update += sizes(aggIdx);
                          });
    aggToRowMap = Kokkos::View<GO*, memory_space>("aggToRowMap", aggStart(numAggregates));

    // count, how many dofs have been recorded for each aggregate so far
    Kokkos::View<LO*, memory_space> numDofs("number of dofs per aggregate", numAggregates); // empty array with number of Dofs for each aggregate

    if (stridedblocksize_ == 1) {
      Kokkos::parallel_for("compute addToRowMap",
                           Kokkos::RangePolicy<execution_space>(0, size),
                           KOKKOS_LAMBDA(const LO lnode) {
                             LO myAgg = vertex2AggId(lnode, 0);
                             if (procWinner(lnode, 0) == myPid) {
                               LO myNumDofs = Kokkos::atomic_fetch_add(&numDofs(myAgg), 1);
                               aggToRowMap(aggStart(myAgg) + myNumDofs)
                                 = ComputeGlobalDOF(nodeLocalMap.getGlobalElement(lnode));
                             }
                           });
    } else {
      Kokkos::parallel_for("compute addToRowMap",
                           Kokkos::RangePolicy<execution_space>(0, size),
                           KOKKOS_LAMBDA(const LO lnode) {
                             LO myAgg = vertex2AggId(lnode, 0);
                             if (procWinner(lnode, 0) == myPid) {
                               GO gnodeid = nodeLocalMap.getGlobalElement(lnode);
                               for (LocalOrdinal k = 0; k < stridedblocksize_; k++) {
                                 // // Note: LBV on 10/11/19, removing the test
                                 // // below as there is no way to refactor it
                                 // // correctly at the moment since the method
                                 // // is not inlined for GPU kernels...
                                 GlobalOrdinal gDofIndex = ComputeGlobalDOF(gnodeid, k);
                                 // if (columnMap_->isNodeGlobalElement(gDofIndex)) {
                                 LO myNumDofs = Kokkos::atomic_fetch_add(&numDofs(myAgg), 1);
                                 aggToRowMap( aggStart(myAgg) + myNumDofs ) = gDofIndex;
                                 // }
                               }
                             }
                           });
    }
    // todo plausibility check: entry numDofs[k] == aggToRowMap[k].size()

  } //UnamalgamateAggregates

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void AmalgamationInfo_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  UnamalgamateAggregatesLO(const Aggregates_kokkos& aggregates,
                           Kokkos::View<LO*, memory_space>& aggStart,
                           Kokkos::View<LO*, memory_space>& aggToRowMap) const {

    const int myPid        = aggregates.GetMap()->getComm()->getRank();
    auto procWinner        = aggregates.GetProcWinner()  ->template getLocalView<memory_space>();
    auto vertex2AggId      = aggregates.GetVertex2AggId()->template getLocalView<memory_space>();
    auto nodeLocalMap      = aggregates.GetMap()->getLocalMap();
    const GO numAggregates = aggregates.GetNumAggregates();


    // FIXME: Do we need to compute size here? Or can we use existing?
    LO size = static_cast<LO>(procWinner.extent(0));

    Kokkos::View<LO*, memory_space> sizes("aggregate sizes", numAggregates);
    if (stridedblocksize_ == 1) {
      Kokkos::parallel_for("compute aggregates sizes",
                           Kokkos::RangePolicy<execution_space>(0, size),
                           KOKKOS_LAMBDA(const LO lnode) {
                             if (procWinner(lnode, 0) == myPid) {
                               Kokkos::atomic_increment(&sizes(vertex2AggId(lnode, 0)));
                             }
                           });
    } else {
      Kokkos::parallel_for("compute aggregates sizes",
                           Kokkos::RangePolicy<execution_space>(0, size),
                           KOKKOS_LAMBDA(const LO lnode) {
                             if (procWinner(lnode, 0) == myPid) {
                               // GO nodeGID = nodeLocalMap.getGlobalElement(lnode);

                               for (LO k = 0; k < stridedblocksize_; ++k) {
                                 // GO GID = ComputeGlobalDOF(nodeGID, k);
                                 // if (columnMap_->isNodeGlobalElement(GID)) {
                                 Kokkos::atomic_increment(&sizes(vertex2AggId(lnode, 0)));
                                 // }
                               }
                             }
                           });
    }

    aggStart = Kokkos::View<LO*, memory_space>("aggStart", numAggregates+1);
    Kokkos::parallel_scan("compute aggregate starts",
                          Kokkos::RangePolicy<execution_space>(0, numAggregates+1),
                          KOKKOS_LAMBDA(const LO aggIdx, LO& update, const bool final) {
                            if(final) {aggStart(aggIdx) = update;}
                            update += sizes(aggIdx);
                          });

    aggToRowMap = Kokkos::View<LO*, memory_space>("aggToRowMap", aggStart(numAggregates));

    // count, how many dofs have been recorded for each aggregate so far
    Kokkos::View<LO*, memory_space> numDofs("number of dofs per aggregate", numAggregates); // empty array with number of DOFs for each aggregate
    if (stridedblocksize_ == 1) {
      Kokkos::parallel_for("compute addToRowMap",
                           Kokkos::RangePolicy<execution_space>(0, size),
                           KOKKOS_LAMBDA(const LO lnode) {
                             if (procWinner(lnode, 0) == myPid) {
                               LO myAgg = vertex2AggId(lnode, 0);
                               const LO myNumDofs = Kokkos::atomic_fetch_add(&numDofs(myAgg), 1);
                               aggToRowMap(aggStart(myAgg) + myNumDofs) = lnode;
                             }
                           });

    } else {
      Kokkos::parallel_for("",
                           Kokkos::RangePolicy<execution_space>(0, size),
                           KOKKOS_LAMBDA(const LO lnode) {
                             if (procWinner(lnode, 0) == myPid) {
                               LO myAgg   = vertex2AggId(lnode, 0);
                               // GO nodeGID = nodeLocalMap.getGlobalElement(lnode);

                               for (LO k = 0; k < stridedblocksize_; k++) {
                                 // GO GID = ComputeGlobalDOF(nodeGID, k);
                                 // if (columnMap_->isNodeGlobalElement(GID)) {
                                 const LO myNumDofs = Kokkos::atomic_fetch_add(&numDofs(myAgg), 1);
                                 aggToRowMap(aggStart(myAgg) + myNumDofs) = lnode*stridedblocksize_ + k;
                                 // }
                               }
                             }
                           });
    }
    // todo plausibility check: entry numDofs[k] == aggToRowMap[k].size()

  } //UnamalgamateAggregates

  /////////////////////////////////////////////////////////////////////////////

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > AmalgamationInfo_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  ComputeUnamalgamatedImportDofMap(const Aggregates_kokkos& aggregates) const {

    RCP<const Map> nodeMap = aggregates.GetMap();

    RCP<std::vector<GO> > myDofGids = Teuchos::rcp(new std::vector<GO>);
    ArrayView<const GO> gEltList = nodeMap->getNodeElementList();
    LO nodeElements = Teuchos::as<LO>(nodeMap->getNodeNumElements());
    if (stridedblocksize_ == 1) {
      for (LO n = 0; n<nodeElements; n++) {
        GO gDofIndex = ComputeGlobalDOF(gEltList[n]);
        myDofGids->push_back(gDofIndex);
      }
    } else {
      for (LO n = 0; n<nodeElements; n++) {
        for (LocalOrdinal k = 0; k < stridedblocksize_; k++) {
          GO gDofIndex = ComputeGlobalDOF(gEltList[n],k);
          if (columnMap_->isNodeGlobalElement(gDofIndex))
            myDofGids->push_back(gDofIndex);
        }
      }
    }

    ArrayRCP<GO> arr_myDofGids = Teuchos::arcp( myDofGids );
    RCP<Map> importDofMap = MapFactory::Build(aggregates.GetMap()->lib(),
                                              Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                              arr_myDofGids(),
                                              aggregates.GetMap()->getIndexBase(),
                                              aggregates.GetMap()->getComm());
    return importDofMap;
  }

  /////////////////////////////////////////////////////////////////////////////

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  KOKKOS_INLINE_FUNCTION
  GlobalOrdinal AmalgamationInfo_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  ComputeGlobalDOF(GlobalOrdinal const &gNodeID, LocalOrdinal const &k) const {
    // here, the assumption is, that the node map has the same indexBase as the dof map
    GlobalOrdinal gDofIndex = offset_
      + (gNodeID-indexBase_)*fullblocksize_  // this is the node map index base
      + nStridedOffset_ + k + indexBase_;    // this is the dof map index base
    return gDofIndex;
  }

} //namespace


#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif /* MUELU_AMALGAMATIONINFO_KOKKOS_DEF_HPP_ */
