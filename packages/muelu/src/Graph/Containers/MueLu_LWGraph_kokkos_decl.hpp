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
#ifndef MUELU_LWGRAPH_KOKKOS_DECL_HPP
#define MUELU_LWGRAPH_KOKKOS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include <Kokkos_StaticCrsGraph.hpp>
#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

#include <Xpetra_ConfigDefs.hpp>   // global_size_t
#include <Xpetra_Map.hpp>

#include "MueLu_VerbosityLevel.hpp"
#include "MueLu_LWGraph_kokkos_fwd.hpp"
#include <MueLu_LocalLWGraph_kokkos_decl.hpp>

#include "MueLu_Exceptions.hpp"

namespace MueLu {

  /*!
    @class LWGraph_kokkos
    @brief Lightweight MueLu representation of a compressed row storage graph

    This class is lightweight in the sense that it holds to local graph
    information. These were built without using fillComplete.
   */
  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  class LWGraph_kokkos;

  // Partial specialization for DeviceType
  template<class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  class LWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>> {
  public:
    using local_ordinal_type  = LocalOrdinal;
    using global_ordinal_type = GlobalOrdinal;
    using execution_space     = typename DeviceType::execution_space;
    using memory_space        = typename DeviceType::memory_space;
    using device_type         = Kokkos::Device<execution_space, memory_space>;
    using node_type           = Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>;
    using local_lw_graph_type = MueLu::LocalLWGraph_kokkos<LocalOrdinal, GlobalOrdinal, node_type>;
    using size_type           = size_t;

    using map_type            = Xpetra::Map<LocalOrdinal, GlobalOrdinal, node_type>;
    using local_graph_type    = typename local_lw_graph_type::local_graph_type;
    using boundary_nodes_type = typename local_lw_graph_type::boundary_nodes_type;

  private:
    // For compatibility
    typedef node_type                                           Node;
#undef MUELU_LWGRAPH_KOKKOS_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

  public:

    //! @name Constructors/Destructors.
    //@{

    //! LWGraph constructor
    //
    // @param[in] graph: local graph of type Kokkos::StaticCrsGraph containing CRS data
    // @param[in] domainMap: non-overlapping (domain) map for graph. Usually provided by AmalgamationFactory stored in UnAmalgamationInfo container
    // @param[in] importMap: overlapping map for graph. Usually provided by AmalgamationFactory stored in UnAmalgamationInfo container
    // @param[in] objectLabel: label string
    LWGraph_kokkos(const local_graph_type&    graph,
                   const RCP<const map_type>& domainMap,
                   const RCP<const map_type>& importMap,
                   const std::string&         objectLabel = "")
      : lclLWGraph_(graph, domainMap), domainMap_(domainMap), importMap_(importMap), objectLabel_(objectLabel) { }

    ~LWGraph_kokkos() = default;
    //@}

    const RCP<const Teuchos::Comm<int> > GetComm() const {
      return domainMap_->getComm();
    }
    const RCP<const Map> GetDomainMap() const {
      return domainMap_;
    }
    //! Return overlapping import map (nodes).
    const RCP<const Map> GetImportMap() const {
      return importMap_;
    }

    //! Return number of graph vertices
    KOKKOS_INLINE_FUNCTION size_type GetNodeNumVertices() const {
      return lclLWGraph_.GetNodeNumVertices();
    }
    //! Return number of graph edges
    KOKKOS_INLINE_FUNCTION size_type GetNodeNumEdges() const {
      return lclLWGraph_.GetNodeNumEdges();
    }

    //! Returns the maximum number of entries across all rows/columns on this node
    KOKKOS_INLINE_FUNCTION size_type getLocalMaxNumRowEntries () const {
      return lclLWGraph_.getLocalMaxNumRowEntries();
    }

    /// Return a simple one-line description of the Graph.
    std::string description() const {
      return "LWGraph (" + objectLabel_ + ")";
    }

    //! Print the Graph with some verbosity level to an FancyOStream object.
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

    local_lw_graph_type& getLocalLWGraph() const {
      return lclLWGraph_;
    }

  private:

    //! Underlying graph (with label)
    mutable local_lw_graph_type   lclLWGraph_;

    //! Graph maps
    const RCP<const map_type>   domainMap_;
    const RCP<const map_type>   importMap_;

    //! Name of this graph.
    const std::string objectLabel_;
  };

}

#define MUELU_LWGRAPH_KOKKOS_SHORT
#endif // MUELU_LWGRAPH_KOKKOS_DECL_HPP
