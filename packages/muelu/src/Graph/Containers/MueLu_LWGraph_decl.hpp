// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_LWGRAPH_DECL_HPP
#define MUELU_LWGRAPH_DECL_HPP

#include <Xpetra_ConfigDefs.hpp>  // global_size_t
#include "MueLu_ConfigDefs.hpp"

#include "MueLu_LWGraph_fwd.hpp"
#include "MueLu_LWGraphBase.hpp"

namespace MueLu {

/*!
   @class LWGraph
   @brief Lightweight MueLu representation of a compressed row storage graph.

   This class is lightweight in the sense that it holds to local graph information.  These were built without using
   fillComplete.
   TODO handle systems
*/
template <class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class LWGraph : public MueLu::LWGraphBase<LocalOrdinal, GlobalOrdinal, Node, true> {
 public:
  using LWGraphBase<LocalOrdinal, GlobalOrdinal, Node, true>::LWGraphBase;
};

}  // namespace MueLu

#define MUELU_LWGRAPH_SHORT
#endif  // MUELU_LWGRAPH_DECL_HPP
