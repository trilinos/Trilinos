#ifndef MUELU_GRAPH_HPP
#define MUELU_GRAPH_HPP

#include <Teuchos_ArrayView.hpp>
#include <Epetra_CrsGraph.h>

#include "MueLu_Exceptions.hpp"

/******************************************************************************
   MueLu representation of a graph. Some of this is redundant with an 
   Epetra_CrsGraph so we might want to clean up. In particular, things
   like VertexNeighbors, NVertices, NEdges, etc. are available somewhere 
   in EGraph (though perhaps protected).              
******************************************************************************/

namespace MueLu {
  
  class Graph : public Teuchos::Describable {
    
  public:

    Graph(const Epetra_CrsMatrix & A, const std::string & objectLabel="");
    ~Graph();

    inline int GetNodeNumVertices() { return nVertices_; }
    inline int GetNodeNumEdges()    { return nEdges_;    }
    inline int GetNodeNumGhost()    { return nGhost_;    }

    inline int GetGlobalNumEdges() { 
      TEST_FOR_EXCEPTION(eGraph_ == NULL, MueLu::Exceptions::RuntimeError, "Graph: eGraph_ == NULL. Cannot call this method");
      return eGraph_->NumGlobalNonzeros(); 
    }

    inline const Epetra_Comm & GetComm() { 
      TEST_FOR_EXCEPTION(eGraph_ == NULL, MueLu::Exceptions::RuntimeError, "Graph: eGraph_ == NULL. Cannot call this method");
      return eGraph_->Comm();
    }

    inline const Epetra_BlockMap & GetDomainMap() { 
      TEST_FOR_EXCEPTION(eGraph_ == NULL, MueLu::Exceptions::RuntimeError, "Graph: eGraph_ == NULL. Cannot call this method");
      return eGraph_->DomainMap();
    }

    inline const Epetra_BlockMap & GetImportMap() { 
      TEST_FOR_EXCEPTION(eGraph_ == NULL, MueLu::Exceptions::RuntimeError, "Graph: eGraph_ == NULL. Cannot call this method");
      return eGraph_->ImportMap();
    }

    //! Return the list of vertices adjacent to the vertex 'v'
    Teuchos::ArrayView<int> getNeighborVertices(int v) {
      int* startPtr = &(vertexNeighbors_[vertexNeighborsPtr_[v]]);
      int  size     = vertexNeighborsPtr_[v+1] - vertexNeighborsPtr_[v];

      return Teuchos::ArrayView<int>(startPtr, size);
    }
    // TODO: use directly a GetRow or something like that.

  private:

    // inline const Epetra_CrsGraph * GetCrsGraph() { return eGraph_; } // not used.

    // inline void GetCrsDataPointers(int** vertexNeighborsPtr, int** vertexNeighbors) { 
    //   *vertexNeighborsPtr = vertexNeighborsPtr_;
    //   *vertexNeighbors    = vertexNeighbors_;
    // }

    int  nVertices_;
    int  nEdges_;
    int  nGhost_;

    int  *vertexNeighborsPtr_;  /* VertexNeighbors[VertexNeighborsPtr[i]:     */
    int  *vertexNeighbors_;     /*                 VertexNeighborsPtr[i+1]-1] */
                                /* corresponds to vertices Adjacent to vertex */
                                /* i in graph                                 */ 

    const Epetra_CrsGraph * eGraph_;

  };

  /**********************************************************************************/
  /* Function to build Graph (take an Epetra_CrsMatrix and                          */
  /* extract out the Epetra_CrsGraph).                                              */
  /**********************************************************************************/
  Graph::Graph(const Epetra_CrsMatrix & A, const std::string & objectLabel)
  {

    setObjectLabel(objectLabel);

    nVertices_          = 0;
    nEdges_             = 0;
    nGhost_             = 0;
    vertexNeighborsPtr_ = NULL;
    vertexNeighbors_    = NULL;
    eGraph_             = NULL;

    if (A.NumMyRows() != 0) {
      nVertices_ = A.NumMyRows();

      {
        double *unusedArg = NULL;
        A.ExtractCrsDataPointers(vertexNeighborsPtr_, vertexNeighbors_, unusedArg);
      }

      nEdges_    = vertexNeighborsPtr_[nVertices_]; // TODO: simplify
      eGraph_    = &(A.Graph());
    }

    /*
      Ray comments about nGhost:
      Graph.NGhost == A.RowMatrixColMap().NumMyElements() - A.OperatorDomainMap().NumMyElements()
      is basically right. But we've had some issues about how epetra handles empty columns.
      Probably worth discussing this with Jonathan and Chris to see if this is ALWAYS right. 
    */
    if (eGraph_ != NULL) {
      nGhost_ = A.RowMatrixColMap().NumMyElements() - A.OperatorDomainMap().NumMyElements();
      if (nGhost_ < 0) nGhost_ = 0;
    }

  }

  Graph::~Graph()
  {

  }

}

#endif

// TODO: use directly a CrsGraph (not a CrsMatrix) ?
