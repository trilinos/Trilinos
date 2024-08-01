// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_LWGRAPH_KOKKOS_DECL_HPP
#define MUELU_LWGRAPH_KOKKOS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include <Xpetra_ConfigDefs.hpp>  // global_size_t

#include "MueLu_LWGraphBase.hpp"
#include "MueLu_LWGraph_kokkos_fwd.hpp"
#include "MueLu_LWGraph_fwd.hpp"

namespace MueLu {

/*!
  @class LWGraph_kokkos
  @brief Lightweight MueLu representation of a compressed row storage graph

  This class is lightweight in the sense that it holds to local graph
  information. These were built without using fillComplete.
 */
template <class LocalOrdinal, class GlobalOrdinal, class Node>
class LWGraph_kokkos : public MueLu::LWGraphBase<LocalOrdinal, GlobalOrdinal, Node, false> {
  using LWGraphBase<LocalOrdinal, GlobalOrdinal, Node, false>::LWGraphBase;

 public:
  RCP<MueLu::LWGraph<LocalOrdinal, GlobalOrdinal, Node> > copyToHost();
};

}  // namespace MueLu

#define MUELU_LWGRAPH_KOKKOS_SHORT
#endif  // MUELU_LWGRAPH_KOKKOS_DECL_HPP
