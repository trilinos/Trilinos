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
#ifndef MUELU_LOCALLWGRAPH_KOKKOS_DECL_HPP
#define MUELU_LOCALLWGRAPH_KOKKOS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include <Kokkos_StaticCrsGraph.hpp>
#include <Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

#include <Xpetra_ConfigDefs.hpp>  // global_size_t
#include <Xpetra_Map_fwd.hpp>

#include "MueLu_VerbosityLevel.hpp"
#include "MueLu_LocalLWGraph_kokkos_fwd.hpp"

#include "MueLu_Exceptions.hpp"

namespace MueLu {

/*!
  @class LocalLWGraph_kokkos
  @brief Lightweight MueLu representation of a compressed row storage graph

  This class is lightweight in the sense that it holds to local graph
  information. These were built without using fillComplete.
 */
template <class LocalOrdinal, class GlobalOrdinal, class Node>
class LocalLWGraph_kokkos;

// Partial specialization for DeviceType
template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
class LocalLWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosDeviceWrapperNode<DeviceType>> {
 public:
  using local_ordinal_type  = LocalOrdinal;
  using global_ordinal_type = GlobalOrdinal;
  using execution_space     = typename DeviceType::execution_space;
  using memory_space        = typename DeviceType::memory_space;
  using device_type         = Kokkos::Device<execution_space, memory_space>;
  using range_type          = Kokkos::RangePolicy<local_ordinal_type, execution_space>;
  using node_type           = Tpetra::KokkosCompat::KokkosDeviceWrapperNode<DeviceType>;
  using size_type           = size_t;

  using local_graph_type    = Kokkos::StaticCrsGraph<LocalOrdinal,
                                                  Kokkos::LayoutLeft,
                                                  device_type, void, size_t>;
  using boundary_nodes_type = Kokkos::View<const bool*, memory_space>;
  using row_type            = Kokkos::View<const LocalOrdinal*, memory_space>;
  using map_type            = Xpetra::Map<LocalOrdinal, GlobalOrdinal, node_type>;

 private:
  // For compatibility
  typedef node_type Node;
#undef MUELU_LOCALLWGRAPH_KOKKOS_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! LocalLWGraph constructor
  //
  // @param[in] graph: local graph of type Kokkos::StaticCrsGraph containing CRS data
  LocalLWGraph_kokkos(const local_graph_type& graph,
                      const RCP<const map_type>& domainMap);

  ~LocalLWGraph_kokkos() = default;
  //@}

  //! Return number of graph vertices
  KOKKOS_INLINE_FUNCTION size_type GetNodeNumVertices() const {
    return graph_.numRows();
  }
  //! Return number of graph edges
  KOKKOS_INLINE_FUNCTION size_type GetNodeNumEdges() const {
    return graph_.row_map(GetNodeNumVertices());
  }

  //! Returns the maximum number of entries across all rows/columns on this node
  KOKKOS_INLINE_FUNCTION size_type getLocalMaxNumRowEntries() const {
    return maxNumRowEntries_;
  }

  //! Return the row pointers of the local graph
  KOKKOS_INLINE_FUNCTION typename local_graph_type::row_map_type getRowPtrs() const {
    return graph_.row_map;
  }

  //! Return the list entries in the local graph
  KOKKOS_INLINE_FUNCTION typename local_graph_type::entries_type getEntries() const {
    return graph_.entries;
  }

  //! Return the list of vertices adjacent to the vertex 'v'.
  // Unfortunately, C++11 does not support the following:
  //    auto getNeighborVertices(LO i) const -> decltype(rowView)
  // auto return with decltype was only introduced in C++14
  KOKKOS_INLINE_FUNCTION
  Kokkos::GraphRowViewConst<local_graph_type> getNeighborVertices(LO i) const {
    auto rowView = graph_.rowConst(i);

    return rowView;
  }

  //! Return true if vertex with local id 'v' is on current process.
  KOKKOS_INLINE_FUNCTION bool isLocalNeighborVertex(LO i) const {
    return i >= minLocalIndex_ && i <= maxLocalIndex_;
  }

  //! Set boolean array indicating which rows correspond to Dirichlet boundaries.
  KOKKOS_INLINE_FUNCTION void SetBoundaryNodeMap(const boundary_nodes_type bndry) {
    dirichletBoundaries_ = bndry;
  }

  //! Returns map with global ids of boundary nodes.
  KOKKOS_INLINE_FUNCTION const boundary_nodes_type GetBoundaryNodeMap() const {
    return dirichletBoundaries_;
  }

  const local_graph_type& getGraph() const {
    return graph_;
  }

 private:
  //! Underlying graph (with label)
  const local_graph_type graph_;

  //! Boolean array marking Dirichlet rows.
  boundary_nodes_type dirichletBoundaries_;

  //! Local index boundaries (cached from domain map)
  LO minLocalIndex_, maxLocalIndex_;
  size_type maxNumRowEntries_;
};

}  // namespace MueLu

#define MUELU_LOCALLWGRAPH_KOKKOS_SHORT
#endif  // MUELU_LOCALLWGRAPH_KOKKOS_DECL_HPP
