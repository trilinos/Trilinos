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

    Graph(const RCP<const CrsGraph> & graph, const std::string & objectLabel="") : graph_(graph) { 
      //setObjectLabel(objectLabel); 
    }
    virtual ~Graph() {}
    
    inline size_t GetNodeNumVertices() const { return graph_->getNodeNumRows(); }
    inline size_t GetNodeNumEdges()    const { return graph_->getNodeNumEntries(); }
    
    inline Xpetra::global_size_t GetGlobalNumEdges() const { return graph_->getGlobalNumEntries(); }

    inline const RCP<const Teuchos::Comm<int> > GetComm() const { return graph_->getComm(); }
    inline const RCP<const Map> GetDomainMap() const { return graph_->getDomainMap(); }
    inline const RCP<const Map> GetImportMap() const { return graph_->getColMap(); }

    //! Return the list of vertices adjacent to the vertex 'v'
    inline Teuchos::ArrayView<const LocalOrdinal> getNeighborVertices(LocalOrdinal v) const { 
      Teuchos::ArrayView<const LocalOrdinal> neighborVertices;
      graph_->getLocalRowView(v, neighborVertices); 
      return neighborVertices;
    }

#ifdef MUELU_UNUSED
    size_t GetNodeNumGhost() const { 
      /*
        Ray's comments about nGhost:
        Graph->NGhost == graph_->RowMatrixColMap()->NumMyElements() - graph_->OperatorDomainMap()->NumMyElements()
        is basically right. But we've had some issues about how epetra handles empty columns.
        Probably worth discussing this with Jonathan and Chris to see if this is ALWAYS right. 
      */
      size_t nGhost = graph_->getColMap()->getNodeNumElements() - graph_->getDomainMap()->getNodeNumElements();
      if (nGhost < 0) nGhost = 0; // FIXME: size_t is unsigned.
      
      return nGhost;
    }
#endif

    /// Return a simple one-line description of this object.
    std::string description() const {
      return "MueLu.description()";
    }

    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    //void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const {
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const {
      MUELU_DESCRIBE;

      if (verbLevel & Parameters0) {
        //out0 << "Prec. type: " << type_ << std::endl;
      }

      if (verbLevel & Parameters1) {
        //out0 << "Linear Algebra: " << toString(lib_) << std::endl;
        //out0 << "PrecType: " << type_ << std::endl;
        //out0 << "Parameter list: " << std::endl; { Teuchos::OSTab tab2(out); out << paramList_; }
        //out0 << "Overlap: " << overlap_ << std::endl;
      }

      if (verbLevel & Debug) {
        graph_->describe(out0, Teuchos::VERB_EXTREME);
      }
    }

  private:

    RCP<const CrsGraph> graph_;

  };

}

#define MUELU_GRAPH_SHORT
#endif
