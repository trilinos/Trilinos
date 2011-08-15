#ifndef CTHULHU_EPETRACRSGRAPH_HPP
#define CTHULHU_EPETRACRSGRAPH_HPP

#include "Cthulhu_EpetraConfigDefs.hpp"

#include <Teuchos_Describable.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_CompileTimeAssert.hpp>

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>
#include <Kokkos_CrsGraph.hpp>
#include <Kokkos_NodeHelpers.hpp>

#include <Epetra_CrsGraph.h>
#include <Epetra_Import.h>
#include <Cthulhu_EpetraMap.hpp>
#include <Cthulhu_EpetraImport.hpp>
#include <Cthulhu_EpetraBlockMap.hpp>//TMP?

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_CrsGraph.hpp"
#include "Cthulhu_Comm.hpp"

namespace Cthulhu {
  
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration
  template <class S, class LO, class GO, class N, class SpMatOps>
  class CrsMatrix;
#endif
  
  class EpetraCrsGraph
    : public CrsGraph<int,int>
  {
                   
  public: 
    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructor with fixed number of indices per row.
    EpetraCrsGraph(const RCP<const Map<int,int> > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);

    //! Constructor with variable number of indices per row.
    EpetraCrsGraph(const RCP<const Map<int,int> > &rowMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile);

    //! Constructor with fixed number of indices per row and specified column map.
    EpetraCrsGraph(const RCP<const Map<int,int> > &rowMap, const RCP<const Map<int,int> > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);

    //! Constructor with variable number of indices per row and specified column map.
    EpetraCrsGraph(const RCP<const Map<int,int> > &rowMap, const RCP<const Map<int,int> > &colMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile);

    EpetraCrsGraph(const Teuchos::RCP<Epetra_CrsGraph> &graph) : graph_(graph) { }

    // !Destructor.
    virtual ~EpetraCrsGraph() { }

    //@}

    //! @name Insertion/Removal Methods
    //@{ 

    //! Insert graph indices, using global IDs.
    void insertGlobalIndices(int globalRow, const ArrayView<const int> &indices);

    //! Insert graph indices, using local IDs.
    void insertLocalIndices(int localRow, const ArrayView<const int> &indices);

    //! Remove all graph indices from the specified local row.
    void removeLocalIndices(int localRow) {  graph_->RemoveMyIndices(localRow); }

    //@}

    //! @name Transformational Methods
    //@{ 

    void fillComplete(const RCP<const Map<int,int> > &domainMap, const RCP<const Map<int,int> > &rangeMap, OptimizeOption os = DoOptimizeStorage) { 
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, domainMap, tDomainMap, "Cthulhu::TpetraCrsMatrix::fillComplete() only accept Cthulhu::TpetraMap as input arguments.");
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, rangeMap,  tRangeMap,  "Cthulhu::TpetraCrsMatrix::fillComplete() only accept Cthulhu::TpetraMap as input arguments.");
      graph_->FillComplete(tDomainMap->getEpetra_BlockMap(), tRangeMap->getEpetra_BlockMap());       //TODO: os
    }

    void fillComplete(OptimizeOption os = DoOptimizeStorage) {  //TODO: os
      graph_->FillComplete(); 
    }

    //@}

    //! @name Methods implementing RowGraph.
    //@{ 

    //! Returns the communicator.
    const RCP<const Comm<int> > getComm() const { 
       
      RCP<const Epetra_Comm> rcpComm = rcpFromRef(graph_->Comm());
      return Epetra2Teuchos_Comm(rcpComm);
    }

    //! Returns the Map that describes the row distribution in this graph.
    const RCP<const Map<int,int> > getRowMap() const { 
       
    
      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(graph_->RowMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
      
      // TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "get*Map() of EpetraCrsGraph()");
      // return Teuchos::null;
    }

    //! Returns the Map that describes the column distribution in this graph.
    const RCP<const Map<int,int> > getColMap() const {  //TODO TODO TODO BlockMap vs Map
       
    
      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(graph_->ColMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }

    //! Returns the Map associated with the domain of this graph.
    const RCP<const Map<int,int> > getDomainMap() const { 
       
    
      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(graph_->DomainMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }

    //! Returns the Map associated with the domain of this graph.
    const RCP<const Map<int,int> > getRangeMap() const { 
      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(graph_->RangeMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }

    //! Returns the importer associated with this graph.
    RCP<const Import<int,int> > getImporter() const { 
      RCP<const Epetra_Import> imp = rcp(new Epetra_Import(*graph_->Importer())); //NOTE: non consitent: return pointer, take ref
      return rcp ( new Cthulhu::EpetraImport(imp) );
    }

    //! Returns the number of global rows in the graph.
    global_size_t getGlobalNumRows() const {  return graph_->NumGlobalRows(); }

    //! Returns the number of global columns in the graph.
    global_size_t getGlobalNumCols() const {  return graph_->NumGlobalCols(); }

    //! Returns the number of graph rows owned on the calling node.
    size_t getNodeNumRows() const {  return graph_->NumMyRows(); }

    //! Returns the number of columns connected to the locally owned rows of this graph.
    size_t getNodeNumCols() const {  return graph_->NumMyCols(); }

    //! Returns the index base for global indices for this graph. 
    int getIndexBase() const {  return graph_->IndexBase(); }

    //! Returns the global number of entries in the graph.
    global_size_t getGlobalNumEntries() const {  return graph_->NumGlobalEntries(); }

    //! Returns the local number of entries in the graph.
    size_t getNodeNumEntries() const {  return graph_->NumMyEntries(); }

    //! Returns the current number of entries on this node in the specified global row.
    size_t getNumEntriesInGlobalRow(int globalRow) const {  return graph_->NumGlobalIndices(globalRow); }

    //! Returns the current number of entries on this node in the specified local row.
    size_t getNumEntriesInLocalRow(int localRow) const {  return graph_->NumMyIndices(localRow); }

    //! Returns the current number of allocated entries for this node in the specified global row .
    size_t getNumAllocatedEntriesInGlobalRow(int globalRow) const {  return graph_->NumAllocatedGlobalIndices(globalRow); }

    //! Returns the current number of allocated entries on this node in the specified local row.
    size_t getNumAllocatedEntriesInLocalRow(int localRow) const {  return graph_->NumAllocatedMyIndices(localRow); }

    //! Returns the number of global diagonal entries, based on global row/column index comparisons. 
    global_size_t getGlobalNumDiags() const {  return graph_->NumGlobalDiagonals(); }

    //! Returns the number of local diagonal entries, based on global row/column index comparisons. 
    size_t getNodeNumDiags() const {  return graph_->NumMyDiagonals(); }

    //! Returns the maximum number of entries across all rows/columns on all nodes. 
    size_t getGlobalMaxNumRowEntries() const {  return graph_->GlobalMaxNumIndices(); }

    //! Returns the maximum number of entries across all rows/columns on this node. 
    size_t getNodeMaxNumRowEntries() const {  return graph_->MaxNumIndices(); } //Note: why is it not *My*MaxNumIndices ?

    //! Indicates whether the graph has a well-defined column map. 
    bool hasColMap() const {  return graph_->HaveColMap(); } 

    //! Indicates whether the graph is lower triangular.
    bool isLowerTriangular() const {  return graph_->LowerTriangular(); }

    //! Indicates whether the graph is upper triangular.
    bool isUpperTriangular() const {  return graph_->UpperTriangular(); }

    //! If graph indices are in the local range, this function returns true. Otherwise, this function returns false. */
    bool isLocallyIndexed() const {  return graph_->IndicesAreLocal(); }

    //! If graph indices are in the global range, this function returns true. Otherwise, this function returns false. */
    bool isGloballyIndexed() const {  return graph_->IndicesAreGlobal(); }

    //! Returns \c true if fillComplete() has been called and the graph is in compute mode.
    bool isFillComplete() const {  return graph_->Filled(); }

    //! Returns \c true if storage has been optimized.
    bool isStorageOptimized() const {  return graph_->StorageOptimized(); }

    //! Extract a const, non-persisting view of global indices in a specified row of the graph.
    void getGlobalRowView(int GlobalRow, ArrayView<const int> &Indices) const;

    //! Extract a const, non-persisting view of local indices in a specified row of the graph.
    void getLocalRowView(int LocalRow, ArrayView<const int> &indices) const;

    //@}

    //! @name Overridden from Teuchos::Describable 
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const {  return "TODO"; }

    //! Print the object with some verbosity level to an FancyOStream object.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const {  } //TODO

    //@}

    RCP< const Epetra_CrsGraph> getEpetra_CrsGraph() const {  return graph_; }
    
  private:
    
    RCP<Epetra_CrsGraph> graph_;
    
  }; // class EpetraCrsGraph

} // namespace Cthulhu

#endif
