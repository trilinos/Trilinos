#include "Xpetra_EpetraCrsGraph.hpp"

#include "Xpetra_Utils.hpp"

namespace Xpetra {

  EpetraCrsGraph::EpetraCrsGraph(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype)
    : graph_(Teuchos::rcp(new Epetra_CrsGraph(Copy, toEpetra(rowMap), maxNumEntriesPerRow, toEpetra(pftype)))) { }
  
  // TODO: convert array size_t to int
  //   EpetraCrsGraph::EpetraCrsGraph(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype)
  //     : graph_(Teuchos::rcp(new Epetra_CrsGraph(Copy, toEpetra(rowMap), NumEntriesPerRowToAlloc.getRawPtr(), toEpetra(pftype)))) { }
  
  EpetraCrsGraph::EpetraCrsGraph(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype)
    : graph_(Teuchos::rcp(new Epetra_CrsGraph(Copy, toEpetra(rowMap), toEpetra(colMap), maxNumEntriesPerRow, toEpetra(pftype)))) { }

  // TODO: convert array size_t to int
  //   EpetraCrsGraph::EpetraCrsGraph(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype)
  //     : graph_(Teuchos::rcp(new Epetra_CrsGraph(Copy, toEpetra(rowMap), toEpetra(colMap), NumEntriesPerRowToAlloc.getRawPtr(), toEpetra(pftype)))) { }

  void EpetraCrsGraph::insertGlobalIndices(int globalRow, const ArrayView<const int> &indices) { 
    int* indices_rawPtr = const_cast<int*>(indices.getRawPtr()); // there is no const in the Epetra interface :(
    XPETRA_ERR_CHECK(graph_->InsertGlobalIndices(globalRow, indices.size(), indices_rawPtr)); 
  }

  void EpetraCrsGraph::insertLocalIndices(int localRow, const ArrayView<const int> &indices) { 
    int* indices_rawPtr = const_cast<int*>(indices.getRawPtr()); // there is no const in the Epetra interface :(
    XPETRA_ERR_CHECK(graph_->InsertMyIndices(localRow, indices.size(), indices_rawPtr)); 
  }

  void EpetraCrsGraph::getGlobalRowView(int GlobalRow, ArrayView<const int> &Indices) const { 
    int      numEntries;
    int    * eIndices;
      
    XPETRA_ERR_CHECK(graph_->ExtractGlobalRowView(GlobalRow, numEntries, eIndices));
    if (numEntries == 0) { eIndices = NULL; } // Cf. TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

    Indices = ArrayView<const int>(eIndices, numEntries);
  }

  void EpetraCrsGraph::getLocalRowView(int LocalRow, ArrayView<const int> &indices) const {
    int      numEntries;
    int    * eIndices;
      
    XPETRA_ERR_CHECK(graph_->ExtractMyRowView(LocalRow, numEntries, eIndices));
    if (numEntries == 0) { eIndices = NULL; } // Cf. TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

    indices = ArrayView<const int>(eIndices, numEntries);
  }

  void EpetraCrsGraph::fillComplete(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap, OptimizeOption os){ 
    graph_->FillComplete(toEpetra(domainMap), toEpetra(rangeMap)); 
    if (os == DoOptimizeStorage) graph_->OptimizeStorage();
  }
  
  void EpetraCrsGraph::fillComplete(OptimizeOption os){ 
    graph_->FillComplete();
    if (os == DoOptimizeStorage) graph_->OptimizeStorage();
  }

  std::string EpetraCrsGraph::description() const { return "NotImplemented"; }
  
  void EpetraCrsGraph::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const { } //TODO: throw exception or implements it

  // TODO: move that elsewhere
  RCP< const CrsGraph<int, int> > toXpetra(const Epetra_CrsGraph &g) {
    RCP<const Epetra_CrsGraph> const_graph = rcp(new Epetra_CrsGraph(g));
    
    RCP<Epetra_CrsGraph> graph = Teuchos::rcp_const_cast<Epetra_CrsGraph>(const_graph); //TODO: can I avoid the const_cast ?
    return rcp ( new Xpetra::EpetraCrsGraph(graph) );
  }
  //

} // namespace Xpetra
