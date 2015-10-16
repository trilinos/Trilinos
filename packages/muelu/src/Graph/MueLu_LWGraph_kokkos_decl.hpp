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
#ifdef HAVE_MUELU_KOKKOS_REFACTOR

#include <Kokkos_StaticCrsGraph.hpp>
#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

#include <Xpetra_ConfigDefs.hpp>   // global_size_t
#include <Xpetra_Map.hpp>

#include "MueLu_LWGraph_kokkos_fwd.hpp"

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
    typedef LocalOrdinal                                        local_ordinal_type;
    typedef GlobalOrdinal                                       global_ordinal_type;
    typedef typename DeviceType::execution_space                execution_space;
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> node_type;

  private:
    typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, node_type> map_type;

    typedef Kokkos::StaticCrsGraph<LocalOrdinal, Kokkos::LayoutLeft, execution_space> local_graph_type;

    // For compatibility
    typedef node_type                                           Node;
#undef MUELU_LWGRAPH_KOKKOS_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

  public:

    //! @name Constructors/Destructors.
    //@{
    LWGraph_kokkos(const local_graph_type&    graph,
                   const RCP<const map_type>& domainMap,
                   const RCP<const map_type>& rangeMap,
                   const std::string&         objectLabel = "");

    ~LWGraph_kokkos() { }
    //@}

    size_t GetNodeNumVertices() const                           { return graph_.row_map.dimension_0()-1; }
    size_t GetNodeNumEdges()    const                           { return graph_.row_map(GetNodeNumVertices()); }

    const RCP<const Teuchos::Comm<int> > GetComm()      const   { return domainMap_->getComm(); }
    const RCP<const Map>                 GetDomainMap() const   { return domainMap_; }
    //! Returns overlapping import map (nodes).
    const RCP<const Map>                 GetImportMap() const   { return importMap_; }

    //! Return the list of vertices adjacent to the vertex 'v'.
    ArrayView<const LO> getNeighborVertices(LO i) const;

    //! Return true if vertex with local id 'v' is on current process.
    bool isLocalNeighborVertex(LO i) const                      { return i >= minLocalIndex_ && i <= maxLocalIndex_; }

    //! Set boolean array indicating which rows correspond to Dirichlet boundaries.
    void SetBoundaryNodeMap(const Kokkos::View<const bool*, typename NO::device_type> bndry)  { dirichletBoundaries_ = bndry; }

    //! Returns the maximum number of entries across all rows/columns on this node
    size_t getNodeMaxNumRowEntries () const                     { return maxNumRowEntries_; }

    //! Returns map with global ids of boundary nodes.
    const Kokkos::View<const bool*, typename NO::device_type> GetBoundaryNodeMap() const        { return dirichletBoundaries_; }

    /// Return a simple one-line description of the Graph.
    std::string description() const                             { return "MueLu.description()"; } //FIXME use object's label

    //! Print the Graph with some verbosity level to an FancyOStream object.
    // void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

  private:

    //! Underlying graph (with label)
    const local_graph_type      graph_;

    //! Graph maps
    const RCP<const map_type>   domainMap_;
    const RCP<const map_type>   importMap_;

    //! Boolean array marking Dirichlet rows.
    Kokkos::View<const bool*, typename NO::device_type>         dirichletBoundaries_;

    //! Local index boundaries (cached from domain map)
    LO     minLocalIndex_, maxLocalIndex_;
    size_t maxNumRowEntries_;

    //! Name of this graph.
    const std::string & objectLabel_;
  };

}

#define MUELU_LWGRAPH_KOKKOS_SHORT
#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_LWGRAPH_KOKKOS_DECL_HPP
