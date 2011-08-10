#ifndef CTHULHU_TPETRACRSGRAPH_HPP
#define CTHULHU_TPETRACRSGRAPH_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifndef HAVE_CTHULHU_TPETRA
#error This file should be included only if HAVE_CTHULHU_TPETRA is defined.
#endif

#include <Teuchos_Describable.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_CompileTimeAssert.hpp>

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>
#include <Kokkos_CrsGraph.hpp>
#include <Kokkos_NodeHelpers.hpp>

#include <Tpetra_CrsGraph.hpp>

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_CrsGraph.hpp"
#include "Cthulhu_TpetraImport.hpp"
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
  template <class LocalOrdinal, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node = Kokkos::DefaultNode::DefaultNodeType,
            class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps >
  class TpetraCrsGraph : public CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> {
                   
   // template <class S, class LO, class GO, class N, class SpMatOps>
   // friend class CrsMatrix;

    // The following typedef are used by the CTHULHU_DYNAMIC_CAST() macro.
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMapClass;
    
  public: 
    //! @name Constructor/Destructor Methods
    //@{ 

  //! Constructor with fixed number of indices per row.
  TpetraCrsGraph(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile)
    { 
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, rowMap, tRowMap, "Cthulhu::TpetraCrsGraph constructors only accept Cthulhu::TpetraMap as input arguments.");
      graph_ = rcp(new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node, LocalMatOps>(tRowMap->getTpetra_Map(), maxNumEntriesPerRow)); //TODO: convert, pftype));
    }

  //! Constructor with variable number of indices per row.
  TpetraCrsGraph(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile)
    { 
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, rowMap, tRowMap, "Cthulhu::TpetraCrsGraph constructors only accept Cthulhu::TpetraMap as input arguments.");
      graph_ = rcp(new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node, LocalMatOps>(tRowMap->getTpetra_Map(), NumEntriesPerRowToAlloc)); //TODO: convert, pftype));
    }
    
  //! Constructor with fixed number of indices per row and specified column map.
  /** The column map will be used to filter any graph indices inserted using insertLocalIndices() or insertGlobalIndices().
   */
  TpetraCrsGraph(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile)
    { 
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, rowMap, tRowMap, "Cthulhu::TpetraCrsGraph constructors only accept Cthulhu::TpetraMap as input arguments.");
      graph_ = rcp(new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node, LocalMatOps>(tRowMap->getTpetra_Map(), colMap, maxNumEntriesPerRow)); //TODO: convert, pftype));
    }

  //! Constructor with variable number of indices per row and specified column map.
  /** The column map will be used to filter any graph indices inserted using insertLocalIndices() or insertGlobalIndices().
   */
  TpetraCrsGraph(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile)
    { 
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, rowMap, tRowMap, "Cthulhu::TpetraCrsGraph constructors only accept Cthulhu::TpetraMap as input arguments.");
      graph_ = rcp(new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node, LocalMatOps>(tRowMap->getTpetra_Map(), colMap, NumEntriesPerRowToAlloc)); //TODO: convert, pftype));
    }
    
    TpetraCrsGraph(const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > &graph) : graph_(graph) {       
       
    }

    // !Destructor.
    inline ~TpetraCrsGraph() {  };

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
    inline void insertGlobalIndices(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &indices) {  graph_->insertGlobalIndices(globalRow, indices); };

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
    inline void insertLocalIndices(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &indices) {  graph_->insertLocalIndices(localRow, indices); };

    //! Remove all graph indices from the specified local row.
    /**
       \pre \c localRow is a local row of this graph.
       \pre <tt>isGloballyIndexed() == false</tt>
       \pre <tt>isStorageOptimized() == false</tt>

       \post <tt>getNumEntriesInLocalRow(localRow) == 0</tt>
       \post <tt>indicesAreAllocated() == true</tt>
       \post <tt>isLocallyIndexed() == true</tt>
    */
    inline void removeLocalIndices(LocalOrdinal localRow) {  graph_->removeLocalIndices(localRow); };

    //@}

    //! @name Transformational Methods
    /** 
        Each of the methods in this group is a global collective. It is
        necessary to call these mehtods on all nodes participating in the
        communicator associated with this graph.
    */
    //@{ 

    //! \brief Communicate non-local contributions to other nodes.
    inline void globalAssemble() {  graph_->globalAssemble(); };

    /*! Resume fill operations.
      After calling fillComplete(), resumeFill() must be called before initiating any changes to the graph.

      resumeFill() may be called repeatedly. 

      \post  <tt>isFillActive() == true<tt>
      \post  <tt>isFillComplete() == false<tt>
    */
    inline void resumeFill() {  graph_->resumeFill(); };

    /*! \brief Signal that data entry is complete, specifying domain and range maps. 

    Off-node indices are distributed (via globalAssemble()), indices are sorted, redundant indices are eliminated, and global indices are transformed to local indices.

    \pre  <tt>isFillActive() == true<tt>
    \pre <tt>isFillComplete()() == false<tt>

    \post <tt>isFillActive() == false<tt>
    \post <tt>isFillComplete() == true<tt>
    \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
    */ 
    inline void fillComplete(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap, OptimizeOption os = DoOptimizeStorage) {  //TODO: os
       
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, domainMap, tDomainMap, "Cthulhu::TpetraCrsGraph::fillComplete() only accept Cthulhu::TpetraMap as input arguments.");
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, rangeMap,  tRangeMap,  "Cthulhu::TpetraCrsGraph::fillComplete() only accept Cthulhu::TpetraMap as input arguments.");
      graph_->fillComplete(tDomainMap->getTpetra_Map(), tRangeMap->getTpetra_Map()); 
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
    inline void fillComplete(OptimizeOption os = DoOptimizeStorage) {  graph_->fillComplete(); }; //TODO: os

    //@}

    //! @name Methods implementing RowGraph.
    //@{ 

    //! Returns the communicator.
    inline const RCP<const Comm<int> > getComm() const {  return graph_->getComm(); };

    //! Returns the underlying node.
    inline RCP<Node> getNode() const {  return graph_->getNode(); };

    //! Returns the Map that describes the row distribution in this graph.
    inline const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getRowMap() const {  return rcp( new TpetraMapClass(graph_->getRowMap())); };

    //! \brief Returns the Map that describes the column distribution in this graph.
    inline const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getColMap() const {  return rcp( new TpetraMapClass(graph_->getColMap())); };

    //! Returns the Map associated with the domain of this graph.
    inline const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const {  return rcp( new TpetraMapClass(graph_->getDomainMap())); };

    //! Returns the Map associated with the domain of this graph.
    inline const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const {  return rcp( new TpetraMapClass(graph_->getRangeMap())); };

    //! Returns the importer associated with this graph.
    inline RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > getImporter() const {  return rcp( new TpetraImport<LocalOrdinal, GlobalOrdinal, Node>(graph_->getImporter())); };

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Returns the exporter associated with this graph.
    inline RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > getExporter() const {  return graph_->getExporter(); };
#endif // CTHULHU_NOT_IMPLEMENTED

    //! Returns the number of global rows in the graph.
    /** Undefined if isFillActive().
     */
    inline global_size_t getGlobalNumRows() const {  return graph_->getGlobalNumRows(); };

    //! \brief Returns the number of global columns in the graph.
    /** Undefined if isFillActive().
     */
    inline global_size_t getGlobalNumCols() const {  return graph_->getGlobalNumCols(); };

    //! Returns the number of graph rows owned on the calling node.
    inline size_t getNodeNumRows() const {  return graph_->getNodeNumRows(); };

    //! Returns the number of columns connected to the locally owned rows of this graph.
    /** Throws std::runtime_error if <tt>hasColMap() == false</tt>
     */
    inline size_t getNodeNumCols() const {  return graph_->getNodeNumCols(); };

    //! Returns the index base for global indices for this graph. 
    inline GlobalOrdinal getIndexBase() const {  return graph_->getIndexBase(); };

    //! Returns the global number of entries in the graph.
    /** Undefined if isFillActive().
     */
    inline global_size_t getGlobalNumEntries() const {  return graph_->getGlobalNumEntries(); };

    //! Returns the local number of entries in the graph.
    inline size_t getNodeNumEntries() const {  return graph_->getNodeNumEntries(); };

    //! \brief Returns the current number of entries on this node in the specified global row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this graph. */
    inline size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const {  return graph_->getNumEntriesInGlobalRow(globalRow); };

    //! Returns the current number of entries on this node in the specified local row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this graph. */
    inline size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const {  return graph_->getNumEntriesInLocalRow(localRow); };

    //! \brief Returns the total number of indices allocated for the graph, across all rows on this node.
    /*! This is the allocation available to the user. Actual allocation may be larger, for example, after 
      calling fillComplete(), and thus this does not necessarily reflect the memory consumption of the 
      this graph.  

      This quantity is computed during the actual allocation. Therefore, if <tt>indicesAreAllocated() == false</tt>,
      this method returns <tt>OrdinalTraits<size_t>::invalid()</tt>.
    */
    inline size_t getNodeAllocationSize() const {  return graph_->getNodeAllocationSize(); };

    //! \brief Returns the current number of allocated entries for this node in the specified global row .
    /** Throws exception std::runtime_error if the specified global row does not belong to this node. */
    inline size_t getNumAllocatedEntriesInGlobalRow(GlobalOrdinal globalRow) const {  return graph_->getNumAllocatedEntriesInGlobalRow(globalRow); };

    //! Returns the current number of allocated entries on this node in the specified local row.
    /** Throws exception std::runtime_error if the specified local row is not valid for this node. */
    inline size_t getNumAllocatedEntriesInLocalRow(LocalOrdinal localRow) const {  return graph_->getNumAllocatedEntriesInLocalRow(localRow); };

    //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
    /** Undefined if isFillActive().
     */
    inline global_size_t getGlobalNumDiags() const {  return graph_->getGlobalNumDiags(); };

    //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
    /** Undefined if isFillActive().
     */
    inline size_t getNodeNumDiags() const {  return graph_->getNodeNumDiags(); };

    //! \brief Returns the maximum number of entries across all rows/columns on all nodes. 
    /** Undefined if isFillActive().
     */
    inline size_t getGlobalMaxNumRowEntries() const {  return graph_->getGlobalMaxNumRowEntries(); };

    //! \brief Returns the maximum number of entries across all rows/columns on this node. 
    /** Undefined if isFillActive().
     */
    inline size_t getNodeMaxNumRowEntries() const {  return graph_->getNodeMaxNumRowEntries(); };

    //! \brief Indicates whether the graph has a well-defined column map. 
    inline bool hasColMap() const {  return graph_->hasColMap(); }; 

    //! \brief Indicates whether the graph is lower triangular.
    /** Undefined if isFillActive().
     */
    inline bool isLowerTriangular() const {  return graph_->isLowerTriangular(); };

    //! \brief Indicates whether the graph is upper triangular.
    /** Undefined if isFillActive().
     */
    inline bool isUpperTriangular() const {  return graph_->isUpperTriangular(); };

    //! \brief If graph indices are in the local range, this function returns true. Otherwise, this function returns false. */
    inline bool isLocallyIndexed() const {  return graph_->isLocallyIndexed(); };

    //! \brief If graph indices are in the global range, this function returns true. Otherwise, this function returns false. */
    inline bool isGloballyIndexed() const {  return graph_->isGloballyIndexed(); };

    //! Returns \c true if fillComplete() has been called and the graph is in compute mode.
    inline bool isFillComplete() const {  return graph_->isFillComplete(); };

    //! Returns \c true if resumeFill() has been called and the graph is in edit mode.
    inline bool isFillActive() const {  return graph_->isFillActive(); };

    //! \brief Returns \c true if storage has been optimized.
    /**
       Optimized storage means that the allocation of each row is equal to the
       number of entries. The effect is that a pass through the matrix, i.e.,
       during a mat-vec, requires minimal memory traffic. One limitation of
       optimized storage is that no new indices can be added to the graph.
    */
    inline bool isStorageOptimized() const {  return graph_->isStorageOptimized(); };

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Returns \c true if the graph was allocated with static data structures.
    inline ProfileType getProfileType() const {  return graph_->getProfileType(); };
#endif

    //! Extract a list of elements in a specified global row of the graph. Put into pre-allocated storage.
    /*!
      \param LocalRow - (In) Global row number for which indices are desired.
      \param Indices - (Out) Global column indices corresponding to values.
      \param NumIndices - (Out) Number of indices.

      Note: A std::runtime_error exception is thrown indices is not large enough to hold the column indices associated
      with row \c GlobalRow. If \c GlobalRow does not belong to this node, then \c indices is unchanged and \c NumIndices is 
      returned as OrdinalTraits<size_t>::invalid().
    */
    inline void getGlobalRowCopy(GlobalOrdinal GlobalRow, 
                                 const ArrayView<GlobalOrdinal> &Indices, 
                                 size_t &NumIndices
                                 ) const {  graph_->getGlobalRowCopy(GlobalRow, Indices, NumIndices); };

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
    inline void getLocalRowCopy(LocalOrdinal LocalRow, 
                                const ArrayView<LocalOrdinal> &indices, 
                                size_t &NumIndices
                                ) const {  graph_->getLocalRowCopy(LocalRow, indices, NumIndices); };

    //! Extract a const, non-persisting view of global indices in a specified row of the graph.
    /*!
      \param GlobalRow - (In) Global row number for which indices are desired.
      \param Indices   - (Out) Global column indices corresponding to values.
      \pre <tt>isLocallyIndexed() == false</tt>
      \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

      Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
    */
    inline void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &Indices) const {  graph_->getGlobalRowView(GlobalRow, Indices); };

    //! Extract a const, non-persisting view of local indices in a specified row of the graph.
    /*!
      \param LocalRow - (In) Local row number for which indices are desired.
      \param Indices  - (Out) Global column indices corresponding to values.
      \pre <tt>isGloballyIndexed() == false</tt>
      \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

      Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
    */
    inline void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices) const {  graph_->getLocalRowView(LocalRow, indices); };

    //@}

    //! @name Overridden from Teuchos::Describable 
    //@{

    /** \brief Return a simple one-line description of this object. */
    inline std::string description() const {  return graph_->description(); };

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    inline void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const {  graph_->describe(out, verbLevel); };

    //@}

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! @name Methods implementing Cthulhu::DistObject
    //@{

    inline bool checkSizes(const DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node>& source) {  return graph_->checkSizes(source); };

    inline void copyAndPermute(const DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node> & source,
                               size_t numSameIDs,
                               const ArrayView<const LocalOrdinal> &permuteToLIDs,
                               const ArrayView<const LocalOrdinal> &permuteFromLIDs) {  graph_->copyAndPermute(sourcem numSameIDs, permuteToLIDs, permuteFromLIDs); };

    inline void packAndPrepare(const DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node> & source,
                               const ArrayView<const LocalOrdinal> &exportLIDs,
                               Array<GlobalOrdinal> &exports,
                               const ArrayView<size_t> & numPacketsPerLID,
                               size_t& constantNumPackets,
                               Distributor &distor) {  graph_->packAndPrepare(source, exportLIDs, exports, numPacketsPerLID, constantNumPackets, distor); };

    inline void unpackAndCombine(const ArrayView<const LocalOrdinal> &importLIDs,
                                 const ArrayView<const GlobalOrdinal> &imports,
                                 const ArrayView<size_t> &numPacketsPerLID,
                                 size_t constantNumPackets,
                                 Distributor &distor,
                                 CombineMode CM) {  graph_->unpackAndCombine(importLIDs, imports, numPacketsPerLID, constantNumPackets, distor, CM); };
    //@}
#endif // CTHULHU_NOT_IMPLEMENTED

    //! \name Advanced methods, at increased risk of deprecation.
    //@{

    //! Get an ArrayRCP of the row-offsets.
    /*! Returns null if optimizeStorage() hasn't been called.
      The returned buffer exists in host-memory.
    */
    inline ArrayRCP<const size_t>       getNodeRowBegs() const {  return graph_->getNodeRowBegs(); };

    //! Get an ArrayRCP of the packed column-indices.
    /*! Returns null if optimizeStorage() hasn't been called.
      The returned buffer exists in host-memory.
    */
    inline ArrayRCP<const LocalOrdinal> getNodePackedIndices() const {  return graph_->getNodePackedIndices(); };

    //@}

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! \name Deprecated methods; will be removed at some point in the near future.
    //@{

    /** \brief Re-allocate the data into contiguous storage.

    This method is deprecated and will be removed in a future version of Cthulhu, as 
    the implementation of storage optimization has been below Cthulhu to Kokkos.

    Currently, the implementation simply calls resumeFill() and then fillComplete(OptimizeStorage). As such, it is 
    required to be called by all nodes that participate in the associated communicator.
    */
    CTHULHU_DEPRECATED inline void optimizeStorage() {  graph_->optimizeStorage(); };

    //! Deprecated. Get a persisting const view of the elements in a specified global row of the graph.
    CTHULHU_DEPRECATED inline ArrayRCP<const GlobalOrdinal> getGlobalRowView(GlobalOrdinal GlobalRow) const {  return graph_->getGlobalRowView(GlobalRow); };

    //! Deprecated. Get a persisting const view of the elements in a specified local row of the graph.
    CTHULHU_DEPRECATED inline ArrayRCP<const LocalOrdinal> getLocalRowView(LocalOrdinal LocalRow) const {  return graph_->getLocalRowView(LocalRow); };

    //@}
#endif

    RCP< const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > getTpetra_CrsGraph() const {  return graph_; }
    
  private:
    
    RCP< Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > graph_;
    
  }; // class TpetraCrsGraph

} // namespace Cthulhu

#define CTHULHU_TPETRACRSGRAPH_SHORT
#endif
