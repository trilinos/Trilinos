#ifndef CTHULHU_EPETRACRSGRAPH_HPP
#define CTHULHU_EPETRACRSGRAPH_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifndef HAVE_CTHULHU_EPETRA
#error This file should be included only if HAVE_CTHULHU_EPETRA is defined.
#endif

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
//#include "Cthulhu_RowGraph.hpp"
//#include "Cthulhu_DistObject.hpp"
//#include "Cthulhu_Util.hpp"

namespace Cthulhu {
  
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration
  template <class S, class LO, class GO, class N, class SpMatOps>
  class CrsMatrix;
#endif
  
  //! \brief A class for constructing and using sparse compressed graphs with row access.
  /*!
   * This class allows the construction of sparse graphs with row-access. 
   * 
   * <b>Local vs. Global</b>
   * 
   * Graph entries can be added using either local or global coordinates for the indices. The 
   * accessors isGloballyIndexed() and isLocallyIndexed() indicate whether the indices are currently
   * stored as global or local indices. Many of the class methods are divided into global and local 
   * versions, which differ only in whether they accept/return indices in the global or local coordinate
   * space. Some of these methods may only be used if the graph coordinates are in the appropriate coordinates.
   * For example, getGlobalRowView() returns a View to the indices in global coordinates; if the indices are 
   * not in global coordinates, then no such View can be created.
   * 
   * The global/local distinction does distinguish between operation on the global/local graph. Almost all methods 
   * operate on the local graph, i.e., the rows of the graph associated with the local node, per the distribution specified
   * by the row map. Access to non-local rows requires performing an explicit communication via the import/export capabilities of the
   * CrsGraph object; see DistObject. However, the method insertGlobalValues() is an exception to this rule, as non-local rows are 
   * allowed to be added via the local graph. These rows are stored in the local graph and communicated to the appropriate node 
   * on the next call to globalAssemble() or fillComplete() (the latter calls the former).
   * 
   */
  class EpetraCrsGraph : public CrsGraph<int,int> {
                   
   // template <class S, class LO, class GO, class N, class SpMatOps>
   // friend class CrsMatrix;

  public: 
    //! @name Constructor/Destructor Methods
    //@{ 

  //! Constructor with fixed number of indices per row.
  EpetraCrsGraph(const RCP<const Map<int,int> > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);

  //! Constructor with variable number of indices per row.
  EpetraCrsGraph(const RCP<const Map<int,int> > &rowMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile);

  //! Constructor with fixed number of indices per row and specified column map.
  /** The column map will be used to filter any graph indices inserted using insertLocalIndices() or insertGlobalIndices().
   */
  EpetraCrsGraph(const RCP<const Map<int,int> > &rowMap, const RCP<const Map<int,int> > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);

  //! Constructor with variable number of indices per row and specified column map.
  /** The column map will be used to filter any graph indices inserted using insertLocalIndices() or insertGlobalIndices().
   */
  EpetraCrsGraph(const RCP<const Map<int,int> > &rowMap, const RCP<const Map<int,int> > &colMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile);

    EpetraCrsGraph(const Teuchos::RCP<Epetra_CrsGraph> &graph) : graph_(graph) { CTHULHU_DEBUG_ME; }

    // !Destructor.
    inline ~EpetraCrsGraph() { CTHULHU_DEBUG_ME; };

    //@}

    //! @name Insertion/Removal Methods
    //@{ 

    //! Insert graph indices, using global IDs.
    /** All index values must be in the global space. 
        \pre \c globalRow exists as an ID in the global row map
        \pre <tt>isLocallyIndexed() == false</tt>
        \pre <tt>isStorageOptimized() == false</tt>

        \post <tt>indicesAreAllocated() == true</tt>
        \post <tt>isGloballyIndexed() == true</tt>

        \note If \c globalRow does not belong to the graph on this node, then it will be communicated to the appropriate node when globalAssemble() is called (which will, at the latest, occur during the next call to fillComplete().) Otherwise, the entries will be inserted in the local graph. 
        \note If the graph row already contains entries at the indices corresponding to values in \c indices, then the redundant indices will be eliminated; this may happen at insertion or during the next call to fillComplete().
    */
    inline void insertGlobalIndices(int globalRow, const ArrayView<const int> &indices) { 
      CTHULHU_DEBUG_ME;

      int* indices_rawPtr = const_cast<int*>(indices.getRawPtr()); // there is no const in the Epetra interface :(
      CTHULHU_ERR_CHECK(graph_->InsertGlobalIndices(globalRow, indices.size(), indices_rawPtr)); 
    };

    //! Insert graph indices, using local IDs.
    /**
       \pre \c localRow is a local row belonging to the graph on this node
       \pre <tt>isGloballyIndexed() == false</tt>
       \pre <tt>isStorageOptimized() == false</tt>
       \pre <tt>hasColMap() == true</tt>

       \post <tt>indicesAreAllocated() == true</tt>
       \post <tt>isLocallyIndexed() == true</tt>

       \note If the graph row already contains entries at the indices corresponding to values in \c indices, then the redundant indices will be eliminated; this may happen at insertion or during the next call to fillComplete().
    */
    inline void insertLocalIndices(int localRow, const ArrayView<const int> &indices) { 
      CTHULHU_DEBUG_ME;

      int* indices_rawPtr = const_cast<int*>(indices.getRawPtr()); // there is no const in the Epetra interface :(
      CTHULHU_ERR_CHECK(graph_->InsertMyIndices(localRow, indices.size(), indices_rawPtr)); 
    }

    //! Remove all graph indices from the specified local row.
    /**
       \pre \c localRow is a local row of this graph.
       \pre <tt>isGloballyIndexed() == false</tt>
       \pre <tt>isStorageOptimized() == false</tt>

       \post <tt>getNumEntriesInLocalRow(localRow) == 0</tt>
       \post <tt>indicesAreAllocated() == true</tt>
       \post <tt>isLocallyIndexed() == true</tt>
    */
    inline void removeLocalIndices(int localRow) { CTHULHU_DEBUG_ME; graph_->RemoveMyIndices(localRow); };

    //@}

    //! @name Transformational Methods
    /** 
        Each of the methods in this group is a global collective. It is
        necessary to call these mehtods on all nodes participating in the
        communicator associated with this graph.
    */
    //@{ 

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! \brief Communicate non-local contributions to other nodes.
    inline void globalAssemble() { CTHULHU_DEBUG_ME; graph_->globalAssemble(); };
#endif

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    /*! Resume fill operations.
      After calling fillComplete(), resumeFill() must be called before initiating any changes to the graph.

      resumeFill() may be called repeatedly. 

      \post  <tt>isFillActive() == true<tt>
      \post  <tt>isFillComplete() == false<tt>
    */
    inline void resumeFill() { CTHULHU_DEBUG_ME; graph_->resumeFill(); };
#endif

    /*! \brief Signal that data entry is complete, specifying domain and range maps. 

    Off-node indices are distributed (via globalAssemble()), indices are sorted, redundant indices are eliminated, and global indices are transformed to local indices.

    \pre  <tt>isFillActive() == true<tt>
    \pre <tt>isFillComplete()() == false<tt>

    \post <tt>isFillActive() == false<tt>
    \post <tt>isFillComplete() == true<tt>
    \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
    */ 
    inline void fillComplete(const RCP<const Map<int,int> > &domainMap, const RCP<const Map<int,int> > &rangeMap, OptimizeOption os = DoOptimizeStorage) { 
      CTHULHU_DEBUG_ME; 
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, domainMap, tDomainMap, "Cthulhu::TpetraCrsMatrix::fillComplete() only accept Cthulhu::TpetraMap as input arguments.");
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, rangeMap,  tRangeMap,  "Cthulhu::TpetraCrsMatrix::fillComplete() only accept Cthulhu::TpetraMap as input arguments.");
      graph_->FillComplete(tDomainMap->getEpetra_BlockMap(), tRangeMap->getEpetra_BlockMap());       //TODO: os
    };

    /*! \brief Signal that data entry is complete. 

    Off-node entries are distributed (via globalAssemble()), repeated entries are summed, and global indices are transformed to local indices.

    \note This method calls fillComplete( getRowMap(), getRowMap(), os ).

    \pre  <tt>isFillActive() == true<tt>
    \pre <tt>isFillComplete()() == false<tt>

    \post <tt>isFillActive() == false<tt>
    \post <tt>isFillComplete() == true<tt>
    \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
    */
    inline void fillComplete(OptimizeOption os = DoOptimizeStorage) {  //TODO: os
      CTHULHU_DEBUG_ME; 
      graph_->FillComplete(); 
    };

    //@}

    //! @name Methods implementing RowGraph.
    //@{ 

    //! Returns the communicator.
    inline const RCP<const Comm<int> > getComm() const { 
      CTHULHU_DEBUG_ME; 
      RCP<const Epetra_Comm> rcpComm = rcpFromRef(graph_->Comm());
      return Epetra2Teuchos_Comm(rcpComm);
    };

#ifdef NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns the underlying node.
    inline RCP<Node> getNode() const { CTHULHU_DEBUG_ME; return graph_->getNode(); };
#endif

    //! Returns the Map that describes the row distribution in this graph.
    inline const RCP<const Map<int,int> > getRowMap() const { 
      CTHULHU_DEBUG_ME; 
    
      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(graph_->RowMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
      
      // TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "get*Map() of EpetraCrsGraph()");
      // return Teuchos::null;
    };

    //! \brief Returns the Map that describes the column distribution in this graph.
    inline const RCP<const Map<int,int> > getColMap() const {  //TODO TODO TODO BlockMap vs Map
      CTHULHU_DEBUG_ME; 
    
      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(graph_->ColMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    };

    //! Returns the Map associated with the domain of this graph.
    inline const RCP<const Map<int,int> > getDomainMap() const { 
      CTHULHU_DEBUG_ME; 
    
      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(graph_->DomainMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    };

    //! Returns the Map associated with the domain of this graph.
    inline const RCP<const Map<int,int> > getRangeMap() const { 
      CTHULHU_DEBUG_ME; 

      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(graph_->RangeMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    };

    //! Returns the importer associated with this graph.
    inline RCP<const Import<int,int> > getImporter() const { 
      CTHULHU_DEBUG_ME; 
      
      RCP<const Epetra_Import> imp = rcp(new Epetra_Import(*graph_->Importer())); //NOTE: non consitent: return pointer, take ref
      return rcp ( new Cthulhu::EpetraImport(imp) );

    };

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Returns the exporter associated with this graph.
    inline RCP<const Export<int,int> > getExporter() const { CTHULHU_DEBUG_ME; return graph_->getExporter(); };
#endif // CTHULHU_NOT_IMPLEMENTED

    //! Returns the number of global rows in the graph.
    /** Undefined if isFillActive().
     */
    inline global_size_t getGlobalNumRows() const { CTHULHU_DEBUG_ME; return graph_->NumGlobalRows(); };

    //! \brief Returns the number of global columns in the graph.
    /** Undefined if isFillActive().
     */
    inline global_size_t getGlobalNumCols() const { CTHULHU_DEBUG_ME; return graph_->NumGlobalCols(); };

    //! Returns the number of graph rows owned on the calling node.
    inline size_t getNodeNumRows() const { CTHULHU_DEBUG_ME; return graph_->NumMyRows(); };

    //! Returns the number of columns connected to the locally owned rows of this graph.
    /** Throws std::runtime_error if <tt>hasColMap() == false</tt>
     */
    inline size_t getNodeNumCols() const { CTHULHU_DEBUG_ME; return graph_->NumMyCols(); };

    //! Returns the index base for global indices for this graph. 
    inline int getIndexBase() const { CTHULHU_DEBUG_ME; return graph_->IndexBase(); };

    //! Returns the global number of entries in the graph.
    /** Undefined if isFillActive().
     */
    inline global_size_t getGlobalNumEntries() const { CTHULHU_DEBUG_ME; return graph_->NumGlobalEntries(); };

    //! Returns the local number of entries in the graph.
    inline size_t getNodeNumEntries() const { CTHULHU_DEBUG_ME; return graph_->NumMyEntries(); };

    //! \brief Returns the current number of entries on this node in the specified global row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this graph. */
    inline size_t getNumEntriesInGlobalRow(int globalRow) const { CTHULHU_DEBUG_ME; return graph_->NumGlobalIndices(globalRow); };

    //! Returns the current number of entries on this node in the specified local row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this graph. */
    inline size_t getNumEntriesInLocalRow(int localRow) const { CTHULHU_DEBUG_ME; return graph_->NumMyIndices(localRow); };

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! \brief Returns the total number of indices allocated for the graph, across all rows on this node.
    /*! This is the allocation available to the user. Actual allocation may be larger, for example, after 
      calling fillComplete(), and thus this does not necessarily reflect the memory consumption of the 
      this graph.  

      This quantity is computed during the actual allocation. Therefore, if <tt>indicesAreAllocated() == false</tt>,
      this method returns <tt>OrdinalTraits<size_t>::invalid()</tt>.
    */
    inline size_t getNodeAllocationSize() const { CTHULHU_DEBUG_ME; return graph_->(); };
#endif

    //! \brief Returns the current number of allocated entries for this node in the specified global row .
    /** Throws exception std::runtime_error if the specified global row does not belong to this node. */
    inline size_t getNumAllocatedEntriesInGlobalRow(int globalRow) const { CTHULHU_DEBUG_ME; return graph_->NumAllocatedGlobalIndices(globalRow); };

    //! Returns the current number of allocated entries on this node in the specified local row.
    /** Throws exception std::runtime_error if the specified local row is not valid for this node. */
    inline size_t getNumAllocatedEntriesInLocalRow(int localRow) const { CTHULHU_DEBUG_ME; return graph_->NumAllocatedMyIndices(localRow); };

    //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
    /** Undefined if isFillActive().
     */
    inline global_size_t getGlobalNumDiags() const { CTHULHU_DEBUG_ME; return graph_->NumGlobalDiagonals(); };

    //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
    /** Undefined if isFillActive().
     */
    inline size_t getNodeNumDiags() const { CTHULHU_DEBUG_ME; return graph_->NumMyDiagonals(); };

    //! \brief Returns the maximum number of entries across all rows/columns on all nodes. 
    /** Undefined if isFillActive().
     */
    inline size_t getGlobalMaxNumRowEntries() const { CTHULHU_DEBUG_ME; return graph_->GlobalMaxNumIndices(); };

    //! \brief Returns the maximum number of entries across all rows/columns on this node. 
    /** Undefined if isFillActive().
     */
    inline size_t getNodeMaxNumRowEntries() const { CTHULHU_DEBUG_ME; return graph_->MaxNumIndices(); }; //Note: why is it not *My*MaxNumIndices ?

    //! \brief Indicates whether the graph has a well-defined column map. 
    inline bool hasColMap() const { CTHULHU_DEBUG_ME; return graph_->HaveColMap(); }; 

    //! \brief Indicates whether the graph is lower triangular.
    /** Undefined if isFillActive().
     */
    inline bool isLowerTriangular() const { CTHULHU_DEBUG_ME; return graph_->LowerTriangular(); };

    //! \brief Indicates whether the graph is upper triangular.
    /** Undefined if isFillActive().
     */
    inline bool isUpperTriangular() const { CTHULHU_DEBUG_ME; return graph_->UpperTriangular(); };

    //! \brief If graph indices are in the local range, this function returns true. Otherwise, this function returns false. */
    inline bool isLocallyIndexed() const { CTHULHU_DEBUG_ME; return graph_->IndicesAreLocal(); };

    //! \brief If graph indices are in the global range, this function returns true. Otherwise, this function returns false. */
    inline bool isGloballyIndexed() const { CTHULHU_DEBUG_ME; return graph_->IndicesAreGlobal(); };

    //! Returns \c true if fillComplete() has been called and the graph is in compute mode.
    inline bool isFillComplete() const { CTHULHU_DEBUG_ME; return graph_->Filled(); };

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns \c true if resumeFill() has been called and the graph is in edit mode.
    inline bool isFillActive() const { CTHULHU_DEBUG_ME; return graph_->isFillActive(); };
#endif

    //! \brief Returns \c true if storage has been optimized.
    /**
       Optimized storage means that the allocation of each row is equal to the
       number of entries. The effect is that a pass through the matrix, i.e.,
       during a mat-vec, requires minimal memory traffic. One limitation of
       optimized storage is that no new indices can be added to the graph.
    */
    inline bool isStorageOptimized() const { CTHULHU_DEBUG_ME; return graph_->StorageOptimized(); };

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Returns \c true if the graph was allocated with static data structures.
    inline ProfileType getProfileType() const { CTHULHU_DEBUG_ME; return graph_->getProfileType(); };
#endif

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Extract a list of elements in a specified global row of the graph. Put into pre-allocated storage.
    /*!
      \param LocalRow - (In) Global row number for which indices are desired.
      \param Indices - (Out) Global column indices corresponding to values.
      \param NumIndices - (Out) Number of indices.

      Note: A std::runtime_error exception is thrown indices is not large enough to hold the column indices associated
      with row \c GlobalRow. If \c GlobalRow does not belong to this node, then \c indices is unchanged and \c NumIndices is 
      returned as OrdinalTraits<size_t>::invalid().
    */
    inline void getGlobalRowCopy(int GlobalRow, 
                                 const ArrayView<int> &Indices, 
                                 size_t &NumIndices
                                 ) const { 
      CTHULHU_DEBUG_ME; 

      graph_->getGlobalRowCopy(GlobalRow, Indices, NumIndices); 


    };
#endif

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Extract a list of elements in a specified local row of the graph. Put into storage allocated by calling routine.
    /*!
      \param LocalRow - (In) Local row number for which indices are desired.
      \param Indices - (Out) Local column indices corresponding to values.
      \param NumIndices - (Out) Number of indices.

      Note: A std::runtime_error exception is thrown indices is not large enough to hold the column indices associated
      with row \c LocalRow. If \c LocalRow is not valid for this node, then \c indices is unchanged and \c NumIndices is 
      returned as OrdinalTraits<size_t>::invalid().

      \pre <tt>isLocallyIndexed()==true</tt> or <tt>hasColMap() == true</tt>
    */
    inline void getLocalRowCopy(int LocalRow, 
                                const ArrayView<int> &indices, 
                                size_t &NumIndices
                                ) const { CTHULHU_DEBUG_ME; graph_->getLocalRowCopy(LocalRow, indices, NumIndices); };
#endif

    //! Extract a const, non-persisting view of global indices in a specified row of the graph.
    /*!
      \param GlobalRow - (In) Global row number for which indices are desired.
      \param Indices   - (Out) Global column indices corresponding to values.
      \pre <tt>isLocallyIndexed() == false</tt>
      \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

      Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
    */
    inline void getGlobalRowView(int GlobalRow, ArrayView<const int> &Indices) const { 
      CTHULHU_DEBUG_ME; 
    
      int      numEntries;
      int    * eIndices;
      
      CTHULHU_ERR_CHECK(graph_->ExtractGlobalRowView(GlobalRow, numEntries, eIndices));
      if (numEntries == 0) { eIndices = NULL; } // Cf. TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

      Indices = ArrayView<const int>(eIndices, numEntries);
    };

    //! Extract a const, non-persisting view of local indices in a specified row of the graph.
    /*!
      \param LocalRow - (In) Local row number for which indices are desired.
      \param Indices  - (Out) Global column indices corresponding to values.
      \pre <tt>isGloballyIndexed() == false</tt>
      \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

      Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
    */
    inline void getLocalRowView(int LocalRow, ArrayView<const int> &indices) const {
      CTHULHU_DEBUG_ME; 
      
      int      numEntries;
      int    * eIndices;
      
      CTHULHU_ERR_CHECK(graph_->ExtractMyRowView(LocalRow, numEntries, eIndices));
      if (numEntries == 0) { eIndices = NULL; } // Cf. TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

      indices = ArrayView<const int>(eIndices, numEntries);
    }

    //@}

    //! @name Overridden from Teuchos::Describable 
    //@{

    /** \brief Return a simple one-line description of this object. */
    inline std::string description() const { CTHULHU_DEBUG_ME; return "TODO"; };

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    inline void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const { CTHULHU_DEBUG_ME; }; //TODO

    //@}

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! @name Methods implementing Cthulhu::DistObject
    //@{

    inline bool checkSizes(const DistObject<int,int,int>& source) { CTHULHU_DEBUG_ME; return graph_->checkSizes(source); };

    inline void copyAndPermute(const DistObject<int,int,int> & source,
                               size_t numSameIDs,
                               const ArrayView<const int> &permuteToLIDs,
                               const ArrayView<const int> &permuteFromLIDs) { CTHULHU_DEBUG_ME; graph_->copyAndPermute(sourcem numSameIDs, permuteToLIDs, permuteFromLIDs); };

    inline void packAndPrepare(const DistObject<int,int,int> & source,
                               const ArrayView<const int> &exportLIDs,
                               Array<int> &exports,
                               const ArrayView<size_t> & numPacketsPerLID,
                               size_t& constantNumPackets,
                               Distributor &distor) { CTHULHU_DEBUG_ME; graph_->packAndPrepare(source, exportLIDs, exports, numPacketsPerLID, constantNumPackets, distor); };

    inline void unpackAndCombine(const ArrayView<const int> &importLIDs,
                                 const ArrayView<const int> &imports,
                                 const ArrayView<size_t> &numPacketsPerLID,
                                 size_t constantNumPackets,
                                 Distributor &distor,
                                 CombineMode CM) { CTHULHU_DEBUG_ME; graph_->unpackAndCombine(importLIDs, imports, numPacketsPerLID, constantNumPackets, distor, CM); };
    //@}
#endif // CTHULHU_NOT_IMPLEMENTED

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! \name Advanced methods, at increased risk of deprecation.
    //@{

    //! Get an ArrayRCP of the row-offsets.
    /*! Returns null if optimizeStorage() hasn't been called.
      The returned buffer exists in host-memory.
    */
    inline ArrayRCP<const size_t>       getNodeRowBegs() const { CTHULHU_DEBUG_ME; return graph_->getNodeRowBegs(); };

    //! Get an ArrayRCP of the packed column-indices.
    /*! Returns null if optimizeStorage() hasn't been called.
      The returned buffer exists in host-memory.
    */
    inline ArrayRCP<const int> getNodePackedIndices() const { CTHULHU_DEBUG_ME; return graph_->getNodePackedIndices(); };

    //@}
#endif

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! \name Deprecated methods; will be removed at some point in the near future.
    //@{

    /** \brief Re-allocate the data into contiguous storage.

    This method is deprecated and will be removed in a future version of Cthulhu, as 
    the implementation of storage optimization has been below Cthulhu to Kokkos.

    Currently, the implementation simply calls resumeFill() and then fillComplete(OptimizeStorage). As such, it is 
    required to be called by all nodes that participate in the associated communicator.
    */
    CTHULHU_DEPRECATED inline void optimizeStorage() { CTHULHU_DEBUG_ME; graph_->optimizeStorage(); };

    //! Deprecated. Get a persisting const view of the elements in a specified global row of the graph.
    CTHULHU_DEPRECATED inline ArrayRCP<const int> getGlobalRowView(int GlobalRow) const { CTHULHU_DEBUG_ME; return graph_->getGlobalRowView(GlobalRow); };

    //! Deprecated. Get a persisting const view of the elements in a specified local row of the graph.
    CTHULHU_DEPRECATED inline ArrayRCP<const int> getLocalRowView(int LocalRow) const { CTHULHU_DEBUG_ME; return graph_->getLocalRowView(LocalRow); };

    //@}
#endif

    RCP< const Epetra_CrsGraph> getEpetra_CrsGraph() const { CTHULHU_DEBUG_ME; return graph_; }
    
  private:
    
    RCP<Epetra_CrsGraph> graph_;
    
  }; // class EpetraCrsGraph

} // namespace Cthulhu

#endif
