#ifndef CTHULHU_CRSGRAPH_HPP
#define CTHULHU_CRSGRAPH_HPP

#include <Teuchos_Describable.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_CompileTimeAssert.hpp>

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>
#include <Kokkos_CrsGraph.hpp>
#include <Kokkos_NodeHelpers.hpp>

#include "Cthulhu_ConfigDefs.hpp"
// #include "Cthulhu_RowGraph.hpp"
// #include "Cthulhu_DistObject.hpp"
// #include "Cthulhu_Util.hpp"
#include "Cthulhu_Import.hpp"

namespace Cthulhu {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration
  template <class S, class LO, class GO, class N, class SpMatOps>
  class CrsMatrix;
#endif

  struct RowInfo {
    size_t localRow;
    size_t allocSize;
    size_t numEntries;
    size_t offset1D;
  };

  enum ELocalGlobal {
    LocalIndices,
    GlobalIndices
  };

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
  class CrsGraph { // : public RowGraph<LocalOrdinal,GlobalOrdinal,Node>,
//                    public DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node> { TODO

    // template <class S, class LO, class GO, class N, class SpMatOps>
    // friend class CrsMatrix;

  public: 
    //! @name Constructor/Destructor Methods
    //@{ 

    // !Destructor.
    virtual ~CrsGraph() {}

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
    virtual void insertGlobalIndices(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &indices) = 0;

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
    virtual void insertLocalIndices(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &indices) = 0;

    //! Remove all graph indices from the specified local row.
    /**
       \pre \c localRow is a local row of this graph.
       \pre <tt>isGloballyIndexed() == false</tt>
       \pre <tt>isStorageOptimized() == false</tt>

       \post <tt>getNumEntriesInLocalRow(localRow) == 0</tt>
       \post <tt>indicesAreAllocated() == true</tt>
       \post <tt>isLocallyIndexed() == true</tt>
    */
    virtual void removeLocalIndices(LocalOrdinal localRow) = 0;

    //@}

    //! @name Transformational Methods
    /** 
        Each of the methods in this group is a global collective. It is
        necessary to call these mehtods on all nodes participating in the
        communicator associated with this graph.
    */
    //@{ 

    /*! \brief Signal that data entry is complete, specifying domain and range maps. 

    Off-node indices are distributed (via globalAssemble()), indices are sorted, redundant indices are eliminated, and global indices are transformed to local indices.

    \pre  <tt>isFillActive() == true<tt>
    \pre <tt>isFillComplete()() == false<tt>

    \post <tt>isFillActive() == false<tt>
    \post <tt>isFillComplete() == true<tt>
    \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
    */ 
    virtual void fillComplete(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap, OptimizeOption os = DoOptimizeStorage) = 0;

    /*! \brief Signal that data entry is complete. 

    Off-node entries are distributed (via globalAssemble()), repeated entries are summed, and global indices are transformed to local indices.

    \note This method calls fillComplete( getRowMap(), getRowMap(), os ).

    \pre  <tt>isFillActive() == true<tt>
    \pre <tt>isFillComplete()() == false<tt>

    \post <tt>isFillActive() == false<tt>
    \post <tt>isFillComplete() == true<tt>
    \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
    */
    virtual void fillComplete(OptimizeOption os = DoOptimizeStorage) = 0;

    //@}

    //! @name Methods implementing RowGraph.
    //@{ 

    //! Returns the communicator.
    virtual const RCP<const Comm<int> > getComm() const = 0;

    //! Returns the Map that describes the row distribution in this graph.
    virtual const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getRowMap() const = 0;

    //! \brief Returns the Map that describes the column distribution in this graph.
    virtual const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getColMap() const = 0;

    //! Returns the Map associated with the domain of this graph.
    virtual const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const = 0;

    //! Returns the Map associated with the domain of this graph.
    virtual const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const = 0;

    //! Returns the importer associated with this graph.
    virtual RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > getImporter() const = 0;

    //! Returns the number of global rows in the graph.
    /** Undefined if isFillActive().
     */
    virtual global_size_t getGlobalNumRows() const = 0;

    //! \brief Returns the number of global columns in the graph.
    /** Undefined if isFillActive().
     */
    virtual global_size_t getGlobalNumCols() const = 0;

    //! Returns the number of graph rows owned on the calling node.
    virtual size_t getNodeNumRows() const = 0;

    //! Returns the number of columns connected to the locally owned rows of this graph.
    /** Throws std::runtime_error if <tt>hasColMap() == false</tt>
     */
    virtual size_t getNodeNumCols() const = 0;

    //! Returns the index base for global indices for this graph. 
    virtual GlobalOrdinal getIndexBase() const = 0;

    //! Returns the global number of entries in the graph.
    /** Undefined if isFillActive().
     */
    virtual global_size_t getGlobalNumEntries() const = 0;

    //! Returns the local number of entries in the graph.
    virtual size_t getNodeNumEntries() const = 0;

    //! \brief Returns the current number of entries on this node in the specified global row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this graph. */
    virtual size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const = 0;

    //! Returns the current number of entries on this node in the specified local row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this graph. */
    virtual size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const = 0;

    //! \brief Returns the current number of allocated entries for this node in the specified global row .
    /** Throws exception std::runtime_error if the specified global row does not belong to this node. */
    virtual size_t getNumAllocatedEntriesInGlobalRow(GlobalOrdinal globalRow) const = 0;

    //! Returns the current number of allocated entries on this node in the specified local row.
    /** Throws exception std::runtime_error if the specified local row is not valid for this node. */
    virtual size_t getNumAllocatedEntriesInLocalRow(LocalOrdinal localRow) const = 0;

    //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
    /** Undefined if isFillActive().
     */
    virtual global_size_t getGlobalNumDiags() const = 0;

    //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
    /** Undefined if isFillActive().
     */
    virtual size_t getNodeNumDiags() const = 0;

    //! \brief Returns the maximum number of entries across all rows/columns on all nodes. 
    /** Undefined if isFillActive().
     */
    virtual size_t getGlobalMaxNumRowEntries() const = 0;

    //! \brief Returns the maximum number of entries across all rows/columns on this node. 
    /** Undefined if isFillActive().
     */
    virtual size_t getNodeMaxNumRowEntries() const = 0;

    //! \brief Indicates whether the graph has a well-defined column map. 
    virtual bool hasColMap() const = 0; 

    //! \brief Indicates whether the graph is lower triangular.
    /** Undefined if isFillActive().
     */
    virtual bool isLowerTriangular() const = 0;

    //! \brief Indicates whether the graph is upper triangular.
    /** Undefined if isFillActive().
     */
    virtual bool isUpperTriangular() const = 0;

    //! \brief If graph indices are in the local range, this function returns true. Otherwise, this function returns false. */
    virtual bool isLocallyIndexed() const = 0;

    //! \brief If graph indices are in the global range, this function returns true. Otherwise, this function returns false. */
    virtual bool isGloballyIndexed() const = 0;

    //! Returns \c true if fillComplete() has been called and the graph is in compute mode.
    virtual bool isFillComplete() const = 0;

    //! \brief Returns \c true if storage has been optimized.
    /**
       Optimized storage means that the allocation of each row is equal to the
       number of entries. The effect is that a pass through the matrix, i.e.,
       during a mat-vec, requires minimal memory traffic. One limitation of
       optimized storage is that no new indices can be added to the graph.
    */
    virtual bool isStorageOptimized() const = 0;

    //! Extract a const, non-persisting view of global indices in a specified row of the graph.
    /*!
      \param GlobalRow - (In) Global row number for which indices are desired.
      \param Indices   - (Out) Global column indices corresponding to values.
      \pre <tt>isLocallyIndexed() == false</tt>
      \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

      Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
    */
    virtual void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &Indices) const = 0;

    //! Extract a const, non-persisting view of local indices in a specified row of the graph.
    /*!
      \param LocalRow - (In) Local row number for which indices are desired.
      \param Indices  - (Out) Global column indices corresponding to values.
      \pre <tt>isGloballyIndexed() == false</tt>
      \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

      Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
    */
    virtual void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices) const = 0;

    //@}

    //! @name Overridden from Teuchos::Describable 
    //@{

    /** \brief Return a simple one-line description of this object. */
    virtual std::string description() const = 0;

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const = 0;

    //@}

  }; // class CrsGraph

} // namespace Cthulhu

#define CTHULHU_CRSGRAPH_SHORT
#endif
