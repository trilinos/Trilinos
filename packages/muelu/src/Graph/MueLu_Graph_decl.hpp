#ifndef MUELU_GRAPH_HPP
#define MUELU_GRAPH_HPP

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_CrsGraph.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_Exceptions.hpp"

/******************************************************************************
   MueLu representation of a graph.
******************************************************************************/

namespace MueLu {

  template <class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class Graph 
    : public BaseClass {

#include "MueLu_UseShortNamesOrdinal.hpp"

  public:

    Graph(const RCP<const CrsGraph> & graph, const std::string & objectLabel="") : graph_(graph) ;
    virtual ~Graph() ;
    
    inline size_t GetNodeNumVertices() const ;
    inline size_t GetNodeNumEdges()    const ;
    
    inline Xpetra::global_size_t GetGlobalNumEdges() const ;

    inline const RCP<const Teuchos::Comm<int> > GetComm() const ;
    inline const RCP<const Map> GetDomainMap() const ;
    inline const RCP<const Map> GetImportMap() const ;

    //! Return the list of vertices adjacent to the vertex 'v'
    inline Teuchos::ArrayView<const LocalOrdinal> getNeighborVertices(LocalOrdinal v) const ;

#ifdef MUELU_UNUSED
    size_t GetNodeNumGhost() const ;
#endif

    /// Return a simple one-line description of this object.
    std::string description() const ;

    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    //void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const ;;

}

#define MUELU_GRAPH_SHORT
#endif
