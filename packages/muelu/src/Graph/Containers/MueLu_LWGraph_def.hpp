// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_LWGRAPH_DEF_HPP
#define MUELU_LWGRAPH_DEF_HPP

#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_LWGraph_decl.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::LWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Node> > MueLu::LWGraph<LocalOrdinal, GlobalOrdinal, Node>::copyToDevice() {
  // This could be improved to skip copies for UVM.

  auto graph = this->getGraph();

  using dev_crs_graph_type = typename MueLu::LWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Node>::crs_graph_type;
  auto rows                = typename dev_crs_graph_type::local_graph_device_type::row_map_type::non_const_type("rows", graph.numRows() + 1);
  auto entries             = typename dev_crs_graph_type::local_graph_device_type::entries_type::non_const_type("columns", graph.entries.extent(0));
  Kokkos::deep_copy(rows, graph.row_map);
  Kokkos::deep_copy(entries, graph.entries);

  using local_graph_type_device = typename MueLu::LWGraphBase<LocalOrdinal, GlobalOrdinal, Node, false>::local_graph_type;
  auto graph_d                  = local_graph_type_device(entries, rows);

  auto lw_d = rcp(new MueLu::LWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Node>(graph_d, this->GetDomainMap(), this->GetImportMap(), this->getObjectLabel()));

  using bndry_nodes_type = typename MueLu::LWGraphBase<LocalOrdinal, GlobalOrdinal, Node, false>::boundary_nodes_type;

  auto bndry   = this->GetBoundaryNodeMap();
  auto bndry_d = bndry_nodes_type("boundary_nodes", bndry.extent(0));
  Kokkos::deep_copy(bndry_d, bndry);
  lw_d->SetBoundaryNodeMap(bndry_d);

  return lw_d;
}

}  // namespace MueLu

#endif  // MUELU_LWGRAPH_DEF_HPP
