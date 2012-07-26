#ifndef MUELU_GRAPH_DEF_HPP
#define MUELU_GRAPH_DEF_HPP

#include <Xpetra_MapFactory.hpp>  // TODO: can go away?

#include "MueLu_Graph_decl.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  Teuchos::ArrayView<const LocalOrdinal> Graph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::getNeighborVertices(LocalOrdinal v) const { 
    Teuchos::ArrayView<const LocalOrdinal> neighborVertices;
    graph_->getLocalRowView(v, neighborVertices); 
    return neighborVertices;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool Graph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::isLocalNeighborVertex(LocalOrdinal v) const {
    return graph_->getDomainMap()->isNodeLocalElement(v);
  }

#ifdef MUELU_UNUSED
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  size_t Graph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetNodeNumGhost() const { 
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
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  std::string Graph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    return "MueLu.description()";
  }

  //! Print the object with some verbosity level to an FancyOStream object.
  //using MueLu::Describable::describe; // overloading, not hiding
  //void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const {
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void Graph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
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

}

#endif // MUELU_GRAPH_DEF_HPP
