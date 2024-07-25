// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  const int GIDListSize                  = localDropMap->getMaxAllGlobalIndex() + 1;

  //  Create a list of GID with only an incomplete/partial set of GIDs, which can then be completed by reduceAll
  Teuchos::Array<GO> partialDropMapGIDList(GIDListSize, -Teuchos::ScalarTraits<GlobalOrdinal>::one());
  Teuchos::Array<GO> redundantDropMapGIDList(GIDListSize, -Teuchos::ScalarTraits<GlobalOrdinal>::one());

  for (GO gid : localDropMapGIDList) {
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

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void checkAggregatesBlockmap(Teuchos::RCP<MueLu::Aggregates<LocalOrdinal, GlobalOrdinal, Node>> aggregates,
                             std::vector<size_t> &stridingInfo,
                             Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> blockmap) {
#include "MueLu_UseShortNames.hpp"

  // Lists of nodes in each aggregate
  struct {
    // GH: For now, copy everything to host until we properly set this factory to run device code
    // Instead, we'll copy data into HostMirrors and run the algorithms on host, saving optimization for later.
    typename Aggregates::LO_view ptr, nodes, unaggregated;
    typename Aggregates::LO_view::HostMirror ptr_h, nodes_h, unaggregated_h;
  } nodesInAgg;
  aggregates->ComputeNodesInAggregate(nodesInAgg.ptr, nodesInAgg.nodes, nodesInAgg.unaggregated);
  nodesInAgg.ptr_h          = Kokkos::create_mirror_view(nodesInAgg.ptr);
  nodesInAgg.nodes_h        = Kokkos::create_mirror_view(nodesInAgg.nodes);
  nodesInAgg.unaggregated_h = Kokkos::create_mirror_view(nodesInAgg.unaggregated);
  Kokkos::deep_copy(nodesInAgg.ptr_h, nodesInAgg.ptr);
  Kokkos::deep_copy(nodesInAgg.nodes_h, nodesInAgg.nodes);
  Kokkos::deep_copy(nodesInAgg.unaggregated_h, nodesInAgg.unaggregated);

  // Create redundant maps to be absolutely sure no entry is missed during the check
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> blockmap_final;
  blockmap_final = MueLuTests::createRedundantMaps<SC, LO, GO, NO>(blockmap);

  // Loop over all aggregates
  for (LO aggregate = 0; aggregate < (LO)nodesInAgg.ptr_h.extent(0) - 1; aggregate++) {
    bool aggContainsNodesInBlockmap    = false;
    bool aggContainsNodesNotInBlockmap = false;
    // Loop over all nodes contained in the aggregate
    for (LO nodeIter = nodesInAgg.ptr_h[aggregate]; nodeIter < (LO)nodesInAgg.ptr_h[aggregate + 1]; nodeIter++) {
      LO node = nodesInAgg.nodes_h[nodeIter];
      // Loop over all node dofs
      for (LO i = 0; i < (LO)stridingInfo[0]; i++) {
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
                               "Aggregate " << aggregate << " crosses contact interface! Not allowed. Error in segregated aggregation procedure. \n")
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void checkAggregatesMapPair(Teuchos::RCP<MueLu::Aggregates<LocalOrdinal, GlobalOrdinal, Node>> aggregates,
                            std::vector<size_t> &stridingInfo,
                            Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> integration_side_map,
                            Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> projection_side_map) {
#include "MueLu_UseShortNames.hpp"

  // Lists of nodes in each aggregate
  struct {
    // GH: For now, copy everything to host until we properly set this factory to run device code
    // Instead, we'll copy data into HostMirrors and run the algorithms on host, saving optimization for later.
    typename Aggregates::LO_view ptr, nodes, unaggregated;
    typename Aggregates::LO_view::HostMirror ptr_h, nodes_h, unaggregated_h;
  } nodesInAgg;
  aggregates->ComputeNodesInAggregate(nodesInAgg.ptr, nodesInAgg.nodes, nodesInAgg.unaggregated);
  nodesInAgg.ptr_h          = Kokkos::create_mirror_view(nodesInAgg.ptr);
  nodesInAgg.nodes_h        = Kokkos::create_mirror_view(nodesInAgg.nodes);
  nodesInAgg.unaggregated_h = Kokkos::create_mirror_view(nodesInAgg.unaggregated);
  Kokkos::deep_copy(nodesInAgg.ptr_h, nodesInAgg.ptr);
  Kokkos::deep_copy(nodesInAgg.nodes_h, nodesInAgg.nodes);
  Kokkos::deep_copy(nodesInAgg.unaggregated_h, nodesInAgg.unaggregated);

  // Create redundant maps to be absolutely sure no entry is missed during the check
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> integration_side_map_final;
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> projection_side_map_final;
  integration_side_map_final = MueLuTests::createRedundantMaps<SC, LO, GO, NO>(integration_side_map);
  projection_side_map_final  = MueLuTests::createRedundantMaps<SC, LO, GO, NO>(projection_side_map);

  // Loop over all aggregates
  for (LO aggregate = 0; aggregate < (LO)nodesInAgg.ptr_h.extent(0) - 1; aggregate++) {
    bool aggContainsNodesInIntegrationSideMap = false;
    bool aggContainsNodesInProjectionSideMap  = false;
    // Loop over all nodes contained in the aggregate
    for (LO nodeIter = nodesInAgg.ptr_h[aggregate]; nodeIter < (LO)nodesInAgg.ptr_h[aggregate + 1]; nodeIter++) {
      LO node = nodesInAgg.nodes_h[nodeIter];
      // Loop over all node dofs
      for (LO i = 0; i < (LO)stridingInfo[0]; i++) {
        // Check if dof is contained in interface mapping of one of the bodies in contact
        if (integration_side_map_final->isNodeGlobalElement(node * stridingInfo[0] + i)) {
          aggContainsNodesInIntegrationSideMap = true;
        }
        if (projection_side_map_final->isNodeGlobalElement(node * stridingInfo[0] + i)) {
          aggContainsNodesInProjectionSideMap = true;
        }
      }
    }
    // Aggregates must not cross from one solid body to next and respect contact interface
    TEUCHOS_TEST_FOR_EXCEPTION((aggContainsNodesInIntegrationSideMap == true) && (aggContainsNodesInProjectionSideMap == true),
                               MueLu::Exceptions::RuntimeError,
                               "Aggregate " << aggregate << " crosses contact interface! Not allowed. Error in segregated aggregation procedure. \n")
  }
}
}  // namespace MueLuTests

#endif