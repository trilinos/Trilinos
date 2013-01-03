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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_GRAPH_DECL_HPP
#define MUELU_GRAPH_DECL_HPP

#include <Xpetra_ConfigDefs.hpp>   // global_size_t
#include <Xpetra_CrsGraph.hpp>     // inline functions requires class declaration
#include <Xpetra_Map_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_BaseClass.hpp"
#include "MueLu_Graph_fwd.hpp"

namespace MueLu {

/*!
   @class Graph
   @brief MueLu representation of a graph.

   This class holds an underlying Xpetra_CrsGraph.
   This class can be considered a facade, as MueLu needs only limited functionality for aggregation.
*/
  template <class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class Graph
    : public BaseClass {
#undef MUELU_GRAPH_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

  public:

    //! @name Constructors/Destructors.
    //@{
    Graph(const RCP<const CrsGraph> & graph, const std::string & objectLabel="") : graph_(graph), gBoundaryNodeMap_(Teuchos::null) {
    }

    virtual ~Graph() {}
    //@}

    size_t GetNodeNumVertices() const { return graph_->getNodeNumRows(); }
    size_t GetNodeNumEdges()    const { return graph_->getNodeNumEntries(); }

    Xpetra::global_size_t GetGlobalNumEdges() const { return graph_->getGlobalNumEntries(); }

    const RCP<const Teuchos::Comm<int> > GetComm() const { return graph_->getComm(); }
    const RCP<const Map> GetDomainMap() const { return graph_->getDomainMap(); }

    //! Returns overlapping import map (nodes).
    const RCP<const Map> GetImportMap() const { return graph_->getColMap();    }

    //! Set map with global ids of boundary nodes.
    void SetBoundaryNodeMap(const RCP<const Map> & gBoundaryNodeMap) { gBoundaryNodeMap_ = gBoundaryNodeMap; }

    //! Returns map with global ids of boundary nodes.
    const RCP<const Map> GetBoundaryNodeMap() const { return gBoundaryNodeMap_; }

    //! Return the list of vertices adjacent to the vertex 'v'.
    Teuchos::ArrayView<const LocalOrdinal> getNeighborVertices(LocalOrdinal v) const;

    //! Return true if vertex with local id 'v' is on current process.
    bool isLocalNeighborVertex(LocalOrdinal v) const;

#ifdef MUELU_UNUSED
    size_t GetNodeNumGhost() const;
#endif

    /// Return a simple one-line description of the Graph.
    std::string description() const;

    //! Print the Graph with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    //void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;;
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

  private:

    RCP<const CrsGraph> graph_;

    //RCP<std::map<GlobalOrdinal,bool> > gBoundaryNodes_;
    //! Vector of global boundary node IDs on current process.
    RCP<const Map> gBoundaryNodeMap_;

  };

} // namespace MueLu

#define MUELU_GRAPH_SHORT
#endif // MUELU_GRAPH_DECL_HPP
