// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2014 Sandia Corporation
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

#ifndef MUELUTESTS_CONTACT_HELPERS_HPP_
#define MUELUTESTS_CONTACT_HELPERS_HPP_

// MueLu
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_SegregatedAFactory.hpp"

// Xpetra
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_StridedMapFactory.hpp>


namespace MueLuTests {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> createRedundantMaps(
        Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> localDropMap) {
#include "MueLu_UseShortNames.hpp"

  Teuchos::RCP<const Teuchos::Comm<int>> comm = localDropMap->getComm();

  Teuchos::Array<GO> localDropMapGIDList = localDropMap->getLocalElementList();
  const int GIDListSize = localDropMap->getMaxAllGlobalIndex()+1;

  //  Create a list of GID with only an incomplete/partial set of GIDs, which can then be completed by reduceAll
  Teuchos::Array<GO> partialDropMapGIDList(GIDListSize, -Teuchos::ScalarTraits<GlobalOrdinal>::one());
  Teuchos::Array<GO> redundantDropMapGIDList(GIDListSize, -Teuchos::ScalarTraits<GlobalOrdinal>::one());

  for(GO gid : localDropMapGIDList){
    partialDropMapGIDList[gid] = gid;
  }

  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, GIDListSize, &partialDropMapGIDList[0],
                     &redundantDropMapGIDList[0]);
  redundantDropMapGIDList.erase(std::remove(redundantDropMapGIDList.begin(), redundantDropMapGIDList.end(), -1),
                                redundantDropMapGIDList.end());
  Teuchos::RCP<const Map> redundantDropMap = MapFactory::Build(
          localDropMap->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), redundantDropMapGIDList, 0, comm);

  return redundantDropMap;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void checkAggregatesBlockmap(Teuchos::RCP<MueLu::Aggregates<LocalOrdinal, GlobalOrdinal, Node>> aggregates,
                             std::vector<size_t> &stridingInfo,
                             Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> blockmap) {
#include "MueLu_UseShortNames.hpp"

  // Check if Graph filtering worked
  typename Aggregates::LO_view aggPtr;
  typename Aggregates::LO_view aggNodes;
  typename Aggregates::LO_view unaggregated;

  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  aggregates->ComputeNodesInAggregate(aggPtr, aggNodes, unaggregated);

  // Create redundant maps to be absolutely sure no entry is missed during the check
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> blockmap_final;
  blockmap_final = MueLuTests::createRedundantMaps<SC,LO,GO,NO>(blockmap);

  // Loop over all aggregates
  for (size_t aggregate = 0; aggregate < aggPtr.extent(0) - 1; aggregate++) {
    bool aggContainsNodesInBlockmap = false;
    bool aggContainsNodesNotInBlockmap = false;
    // Loop over all nodes contained in the aggregate
    for (LO nodeIter = aggPtr[aggregate]; nodeIter < aggPtr[aggregate + 1]; nodeIter++) {
      LO node = aggNodes[nodeIter];
      // Loop over all node dofs
      for (size_t i = 0; i < stridingInfo[0]; i++) {
        // Check if dof is contained in rowmap (inner + interface dofs) of one of the bodies in contact
        if (blockmap_final->isNodeGlobalElement(node * stridingInfo[0] + i)) {
          aggContainsNodesInBlockmap = true;
        } else {
          aggContainsNodesNotInBlockmap = true;
        }
      }
    }
    // Aggregates must not cross from one solid body to next and respect contact interface
    TEUCHOS_TEST_FOR_EXCEPTION(aggContainsNodesInBlockmap == aggContainsNodesNotInBlockmap,
                               MueLu::Exceptions::RuntimeError,
                               "Aggregate " << aggregate <<
                                            " crosses contact interface! Not allowed. Error in segregated aggregation procedure. \n")
  }
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void checkAggregatesMapPair(Teuchos::RCP<MueLu::Aggregates<LocalOrdinal, GlobalOrdinal, Node>> aggregates,
                            std::vector<size_t> &stridingInfo,
                            Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> slavemap,
                            Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> mastermap) {
#include "MueLu_UseShortNames.hpp"

  // Check if Graph filtering worked
  typename Aggregates::LO_view aggPtr;
  typename Aggregates::LO_view aggNodes;
  typename Aggregates::LO_view unaggregated;

  aggregates->ComputeNodesInAggregate(aggPtr, aggNodes, unaggregated);

  // Create redundant maps to be absolutely sure no entry is missed during the check
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> slavemap_final;
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> mastermap_final;
  slavemap_final = MueLuTests::createRedundantMaps<SC,LO,GO,NO>(slavemap);
  mastermap_final = MueLuTests::createRedundantMaps<SC,LO,GO,NO>(mastermap);

  // Loop over all aggregates
  for (size_t aggregate = 0; aggregate < aggPtr.extent(0) - 1; aggregate++) {
    bool aggContainsNodesInSlavemap = false;
    bool aggContainsNodesInMastermap = false;
    // Loop over all nodes contained in the aggregate
    for (LO nodeIter = aggPtr[aggregate]; nodeIter < aggPtr[aggregate + 1]; nodeIter++) {
      LO node = aggNodes[nodeIter];
      // Loop over all node dofs
      for (size_t i = 0; i < stridingInfo[0]; i++) {
        // Check if dof is contained in interface mapping of one of the bodies in contact
        if (slavemap_final->isNodeGlobalElement(node * stridingInfo[0] + i)) {
          aggContainsNodesInSlavemap = true;
        }
        if (mastermap_final->isNodeGlobalElement(node * stridingInfo[0] + i)) {
          aggContainsNodesInMastermap = true;
        }
      }
    }
    // Aggregates must not cross from one solid body to next and respect contact interface
    TEUCHOS_TEST_FOR_EXCEPTION((aggContainsNodesInSlavemap == true) && (aggContainsNodesInMastermap == true),
                               MueLu::Exceptions::RuntimeError,
                               "Aggregate " << aggregate <<
                               " crosses contact interface! Not allowed. Error in segregated aggregation procedure. \n")
  }
}
}

#endif