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
#ifndef MUELU_LWGRAPH_DECL_HPP
#define MUELU_LWGRAPH_DECL_HPP

#include <Xpetra_ConfigDefs.hpp>   // global_size_t
#include <Xpetra_CrsGraph.hpp>     // inline functions requires class declaration
#include <Xpetra_Map_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_LWGraph_fwd.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

/*!
   @class LWGraph
   @brief Lightweight MueLu representation of a compressed row storage graph.

   This class is lightweight in the sense that it holds to local graph information.  These were built without using
   fillComplete.
   TODO handle systems
*/
  template <class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class LWGraph
    : public MueLu::GraphBase<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> { //FIXME  shortnames isn't working
#undef MUELU_LWGRAPH_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

  public:

    //! @name Constructors/Destructors.
    //@{
    LWGraph(const ArrayRCP<const LocalOrdinal> & rowPtrs, const ArrayRCP<const LocalOrdinal> & colPtrs,
            const RCP<const Map>& domainMap, const RCP<const Map>& rangeMap, std::string const & objectLabel="")
            : rows_(rowPtrs), columns_(colPtrs), domainMap_(domainMap), importMap_(rangeMap), objectLabel_(objectLabel) {}

    virtual ~LWGraph() {}
    //@}

    size_t GetNodeNumVertices() const { return rows_.size()-1; }
    size_t GetNodeNumEdges()    const { return rows_[rows_.size()-1]; }

    // TODO: do we really need this function
    // It is being called from CoupledAggregation, but do we need it there?
    Xpetra::global_size_t GetGlobalNumEdges() const {
      Xpetra::global_size_t in = GetNodeNumEdges(), out;
      Teuchos::reduceAll(*domainMap_->getComm(), Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));
      return out;
    }

    const RCP<const Teuchos::Comm<int> > GetComm() const { return domainMap_->getComm(); }
    const RCP<const Map> GetDomainMap() const            { return domainMap_; }
    //! Returns overlapping import map (nodes).
    const RCP<const Map> GetImportMap() const            { return importMap_; }

    void SetBoundaryNodeMap(RCP<const Map> const &map)   {throw(Exceptions::NotImplemented("LWGraph: Boundary node map not implemented."));}

    //! Return the list of vertices adjacent to the vertex 'v'.
    Teuchos::ArrayView<const LocalOrdinal> getNeighborVertices(LocalOrdinal v) const;

    //! Return true if vertex with local id 'v' is on current process.
    bool isLocalNeighborVertex(LocalOrdinal v) const;

    //! Set boolean array indicating which rows correspond to Dirichlet boundaries.
    void SetBoundaryNodeMap(const ArrayRCP<const bool>& bndry) { dirichletBoundaries_ = bndry; }

    //! Returns map with global ids of boundary nodes.
    const ArrayRCP<const bool> GetBoundaryNodeMap() const { return dirichletBoundaries_; }


    /// Return a simple one-line description of the Graph.
    std::string description() const;

    //! Print the Graph with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    //void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;;
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

  private:

    //! Indices into columns_ array.  Part of local graph information.
    const ArrayRCP<const LocalOrdinal> rows_;
    //! Columns corresponding to connections.  Part of local graph information.
    const ArrayRCP<const LocalOrdinal> columns_;
    //! Graph maps
    const RCP<const Map> domainMap_, importMap_;
    //! Name of this graph.
    const std::string & objectLabel_;
    //! Boolean array marking Dirichlet rows.
    ArrayRCP<const bool> dirichletBoundaries_;

  };

} // namespace MueLu

#define MUELU_LWGRAPH_SHORT
#endif // MUELU_LWGRAPH_DECL_HPP
