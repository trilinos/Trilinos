#ifndef MUELU_GRAPH_HPP
#define MUELU_GRAPH_HPP

#include <Teuchos_ArrayView.hpp>
#include <Epetra_CrsGraph.h>

#include "MueLu_Exceptions.hpp"

/******************************************************************************
   MueLu representation of a graph.
******************************************************************************/

namespace MueLu {
  
  class Graph : public Teuchos::Describable {
    
  public:

    Graph(const Epetra_CrsGraph & graph, const std::string & objectLabel="") : graph_(graph) { setObjectLabel(objectLabel); }
    ~Graph() {}
    
    inline int GetNodeNumVertices() { return graph_.NumMyRows(); }
    inline int GetNodeNumEdges()    { return graph_.NumMyNonzeros();  }
    
    inline int GetGlobalNumEdges()  { return graph_.NumGlobalNonzeros(); }

    inline const Epetra_Comm & GetComm() { return graph_.Comm(); }
    inline const Epetra_BlockMap & GetDomainMap() { return graph_.DomainMap(); }
    inline const Epetra_BlockMap & GetImportMap() { return graph_.ImportMap(); }

    //! Return the list of vertices adjacent to the vertex 'v'
    Teuchos::ArrayView<int> getNeighborVertices(int v) {
      int* startPtr;
      int  size;
      
      graph_.ExtractMyRowView(v, size, startPtr);
      return Teuchos::ArrayView<int>(startPtr, size);
    }

    int GetNodeNumGhost()    { 
      /*
        Ray comments about nGhost:
        Graph.NGhost == graph_.RowMatrixColMap().NumMyElements() - graph_.OperatorDomainMap().NumMyElements()
        is basically right. But we've had some issues about how epetra handles empty columns.
        Probably worth discussing this with Jonathan and Chris to see if this is ALWAYS right. 
      */
      int nGhost = graph_.ColMap().NumMyElements() - graph_.DomainMap().NumMyElements();
      if (nGhost < 0) nGhost = 0;
      
      return nGhost;
    }

  private:

    const Epetra_CrsGraph & graph_;

  };

}

#endif
