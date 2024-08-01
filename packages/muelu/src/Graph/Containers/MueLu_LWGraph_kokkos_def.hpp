// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_LWGRAPH_KOKKOS_DEF_HPP
#define MUELU_LWGRAPH_KOKKOS_DEF_HPP

#include "MueLu_LWGraph.hpp"
#include "MueLu_LWGraph_kokkos_decl.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::LWGraph<LocalOrdinal, GlobalOrdinal, Node> > MueLu::LWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Node>::copyToHost() {
  auto graph = this->getGraph();

  auto row_map_h = Kokkos::create_mirror_view(graph.row_map);
  auto entries_h = Kokkos::create_mirror_view(graph.entries);
  Kokkos::deep_copy(row_map_h, graph.row_map);
  Kokkos::deep_copy(entries_h, graph.entries);

  using local_graph_type_host = typename MueLu::LWGraphBase<LocalOrdinal, GlobalOrdinal, Node, true>::local_graph_type;
  auto graph_h                = local_graph_type_host(entries_h, row_map_h);

  auto lw_h = rcp(new MueLu::LWGraph<LocalOrdinal, GlobalOrdinal, Node>(graph_h, this->GetDomainMap(), this->GetImportMap(), this->getObjectLabel()));

  using bndry_nodes_type = typename MueLu::LWGraphBase<LocalOrdinal, GlobalOrdinal, Node, true>::boundary_nodes_type;

  auto bndry   = this->GetBoundaryNodeMap();
  auto bndry_h = bndry_nodes_type("boundary_nodes", bndry.extent(0));
  Kokkos::deep_copy(bndry_h, bndry);
  lw_h->SetBoundaryNodeMap(bndry_h);

  return lw_h;
}

}  // namespace MueLu

#endif  // MUELU_LWGRAPH_KOKKOS_DEF_HPP
