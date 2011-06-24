#ifndef CTHULHU_CRSMATRIX_HPP
#define CTHULHU_CRSMATRIX_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>
#include <Kokkos_CrsMatrix.hpp>

#include "Cthulhu_ConfigDefs.hpp"
//#include "Cthulhu_RowMatrix.hpp"
//TODO #include "Cthulhu_DistObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Cthulhu_CrsGraph.hpp"
//TODO ??? #include "Cthulhu_CrsMatrixMultiplyOp_decl.hpp"

#include "Cthulhu_Map.hpp"
#include "Cthulhu_Vector.hpp"
#include "Cthulhu_MultiVector.hpp"

//#include "Tpetra_Map.hpp" //TODO TMP

#include "Cthulhu_Debug.hpp"

namespace Cthulhu {

  //! \brief A class for constructing and using sparse compressed matrices with row access.
  /*!
   * This class allows the construction of sparse matrices with row-access. 
   * 
   * <b>Local vs. Global</b>
   * 
   * Matrix entries can be added using either local or global coordinates for the indices. The 
   * accessors isGloballyIndexed() and isLocallyIndexed() indicate whether the indices are currently
   * stored as global or local indices. Many of the class methods are divided into global and local 
   * versions, which differ only in whether they accept/return indices in the global or local coordinate
   * space. Some of these methods may only be used if the matrix coordinates are in the appropriate coordinates.
   * For example, getGlobalRowView() returns a View to the indices in global coordinates; if the indices are 
   * not in global coordinates, then no such View can be created.
   * 
   * The global/local distinction does distinguish between operation on the global/local matrix. Almost all methods 
   * operate on the local matrix, i.e., the rows of the matrix associated with the local node, per the distribution specified
   * by the row map. Access to non-local rows requires performing an explicit communication via the import/export capabilities of the
   * CrsMatrix object; see DistObject. However, the method insertGlobalValues() is an exception to this rule, as non-local rows are 
   * allowed to be added via the local matrix. These rows are stored in the local matrix and communicated to the appropriate node 
   * on the next call to globalAssemble() or fillComplete() (the latter calls the former).
   * 
   */
  template <class Scalar, 
            class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType, 
            class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps >
  class CrsMatrix // : public RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>//, public DistObject<char, LocalOrdinal,GlobalOrdinal,Node> { //TODO
    : public Teuchos::LabeledObject { // TODO: Describable

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    // !Destructor.
    virtual ~CrsMatrix() { CTHULHU_DEBUG_ME; }

    //@}

    //! @name Insertion/Removal Methods
    //@{ 

    //! Insert matrix entries, using global IDs.
    /** All index values must be in the global space. 
        \pre \c globalRow exists as an ID in the global row map
        \pre <tt>isLocallyIndexed() == false</tt>
        \pre <tt>isStorageOptimized() == false</tt>

        \post <tt>isGloballyIndexed() == true</tt>

        \note If \c globalRow does not belong to the matrix on this node, then it will be communicated to the appropriate node when globalAssemble() is called (which will, at the latest, occur during the next call to fillComplete().) Otherwise, the entries will be inserted in the local matrix. 
        \note If the matrix row already contains values at the indices corresponding to values in \c cols, then the new values will be summed with the old values; this may happen at insertion or during the next call to fillComplete().
        \note If <tt>hasColMap() == true</tt>, only (cols[i],vals[i]) where cols[i] belongs to the column map on this node will be inserted into the matrix.
    */
    virtual void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &cols, const ArrayView<const Scalar> &vals) =0;

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Insert matrix entries, using local IDs.
    /**
       \pre \c localRow is a local row belonging to the matrix on this node
       \pre <tt>isGloballyIndexed() == false</tt>
       \pre <tt>isStorageOptimized() == false</tt>
       \pre <tt>hasColMap() == true</tt>

       \post <tt>isLocallyIndexed() == true</tt>

       \note If the matrix row already contains entries at the indices corresponding to values in \c cols, then the new values will be summed with the old values; this may happen at insertion or during the next call to fillComplete().
       \note If <tt>hasColMap() == true</tt>, only (cols[i],vals[i]) where cols[i] belongs to the column map on this node will be inserted into the matrix.
    */
    virtual void insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &cols, const ArrayView<const Scalar> &vals) =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! \brief Replace matrix entries, using global IDs.
    /** All index values must be in the global space. 

    \pre \c globalRow is a global row belonging to the matrix on this node.

    \note If (globalRow,cols[i]) corresponds to an entry that is duplicated in this matrix row (likely because it was inserted more than once and fillComplete() has not been called in the interim), the behavior of this function is not defined. */
    virtual void replaceGlobalValues(GlobalOrdinal globalRow, 
                                     const ArrayView<const GlobalOrdinal> &cols,
                                     const ArrayView<const Scalar>        &vals) =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Replace matrix entries, using local IDs.
    /** All index values must be in the local space. 
     */
    virtual void replaceLocalValues(LocalOrdinal localRow, 
                                    const ArrayView<const LocalOrdinal> &cols,
                                    const ArrayView<const Scalar>       &vals) =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Sum into multiple entries, using global IDs.
    /** All index values must be in the global space. 

    \pre \c globalRow is a global row belonging to the matrix on this node.

    */
    virtual void sumIntoGlobalValues(GlobalOrdinal globalRow, 
                                     const ArrayView<const GlobalOrdinal> &cols,
                                     const ArrayView<const Scalar>        &vals) =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Sum into multiple entries, using local IDs.
    /** All index values must be in the local space. 

    \pre \c localRow is a local row belonging to the matrix on this node.

    */
    virtual void sumIntoLocalValues(LocalOrdinal globalRow, 
                                    const ArrayView<const LocalOrdinal>  &cols,
                                    const ArrayView<const Scalar>        &vals) =0; 
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Set all matrix entries equal to scalarThis.
    virtual void setAllToScalar(const Scalar &alpha) =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Scale the current values of a matrix, this = alpha*this. 
    virtual void scale(const Scalar &alpha) =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //@}

    //! @name Transformational Methods
    //@{ 

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! \brief Communicate non-local contributions to other nodes.
    virtual void globalAssemble() =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    /*! Resume fill operations.
      After calling fillComplete(), resumeFill() must be called before initiating any changes to the matrix.

      resumeFill() may be called repeatedly. 

      \post  <tt>isFillActive() == true<tt>
      \post  <tt>isFillComplete() == false<tt>
    */
    virtual void resumeFill() =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    /*! \brief Signal that data entry is complete, specifying domain and range maps.

    Off-node indices are distributed (via globalAssemble()), indices are sorted, redundant indices are eliminated, and global indices are transformed to local indices.

    \pre  <tt>isFillActive() == true<tt>
    \pre <tt>isFillComplete()() == false<tt>

    \post <tt>isFillActive() == false<tt>
    \post <tt>isFillComplete() == true<tt>
    \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
    */ 
    virtual void fillComplete(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap, OptimizeOption os = DoOptimizeStorage) =0;

    /*! \brief Signal that data entry is complete. 

    Off-node entries are distributed (via globalAssemble()), repeated entries are summed, and global indices are transformed to local indices.

    \note This method calls fillComplete( getRowMap(), getRowMap(), os ).

    \pre  <tt>isFillActive() == true<tt>
    \pre <tt>isFillComplete()() == false<tt>

    \post <tt>isFillActive() == false<tt>
    \post <tt>isFillComplete() == true<tt>
    \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
    */
    //TODO : Get ride of "Tpetra"::OptimizeOption
    virtual void fillComplete(Cthulhu::OptimizeOption os = Cthulhu::DoOptimizeStorage) =0;

    //@}

    //! @name Methods implementing RowMatrix
    //@{ 

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns the communicator.
    virtual const RCP<const Comm<int> > getComm() const =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns the underlying node.
    virtual RCP<Node> getNode() const =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //! Returns the Map that describes the row distribution in this matrix.
    virtual const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getRowMap() const =0;

    //! \brief Returns the Map that describes the column distribution in this matrix.
    virtual const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getColMap() const =0;

#ifdef CTHULHU_NOT_IMPLEMENTED    
    //! Returns the RowGraph associated with this matrix. 
    virtual RCP<const RowGraph<LocalOrdinal,GlobalOrdinal,Node> > getGraph() const =0;
#endif // CTHULHU_NOT_IMPLEMENTED

    //! Returns the CrsGraph associated with this matrix. 
    virtual RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > getCrsGraph() const =0;

    //! Returns the number of global rows in this matrix.
    /** Undefined if isFillActive().
     */
    virtual global_size_t getGlobalNumRows() const =0;

    //! \brief Returns the number of global columns in the matrix.
    /** Undefined if isFillActive().
     */
    virtual global_size_t getGlobalNumCols() const =0;

    //! Returns the number of matrix rows owned on the calling node.
    virtual size_t getNodeNumRows() const =0;

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns the number of columns connected to the locally owned rows of this matrix.
    /** Throws std::runtime_error if <tt>hasColMap() == false</tt>
     */
    virtual size_t getNodeNumCols() const =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns the index base for global indices for this matrix. 
    virtual GlobalOrdinal getIndexBase() const =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //! Returns the global number of entries in this matrix.
    virtual global_size_t getGlobalNumEntries() const =0;

    //! Returns the local number of entries in this matrix.
    virtual size_t getNodeNumEntries() const =0;

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! \brief Returns the current number of entries on this node in the specified global row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this matrix. */
    virtual size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //! Returns the current number of entries on this node in the specified local row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
    virtual size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const =0;

    //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
    /** Undefined if isFillActive().
     */
    virtual global_size_t getGlobalNumDiags() const =0;

    //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
    /** Undefined if isFillActive().
     */
    virtual size_t getNodeNumDiags() const =0;

    //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
    /** Undefined if isFillActive().
     */
    virtual size_t getGlobalMaxNumRowEntries() const =0;

    //! \brief Returns the maximum number of entries across all rows/columns on this node.
    /** Undefined if isFillActive().
     */
    virtual size_t getNodeMaxNumRowEntries() const =0;

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! \brief Indicates whether the matrix has a well-defined column map. 
    virtual bool hasColMap() const =0; 
#endif //CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! \brief Indicates whether the matrix is lower triangular.
    /** Undefined if isFillActive().
     */
    virtual bool isLowerTriangular() const =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! \brief Indicates whether the matrix is upper triangular.
    /** Undefined if isFillActive().
     */
    virtual bool isUpperTriangular() const =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
    virtual bool isLocallyIndexed() const =0;

    //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
    virtual bool isGloballyIndexed() const =0;

    //! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
    virtual bool isFillComplete() const =0;

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns \c true if resumeFill() has been called and the matrix is in edit mode.
    virtual bool isFillActive() const =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! \brief Returns \c true if storage has been optimized.
    /**
       Optimized storage means that the allocation of each row is equal to the
       number of entries. The effect is that a pass through the matrix, i.e.,
       during a mat-vec, requires minimal memory traffic. One limitation of
       optimized storage is that no new indices can be added to the matrix.
    */
    virtual bool isStorageOptimized() const =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns \c true if the matrix was allocated with static data structures.
    virtual Tpetra::ProfileType getProfileType() const =0; //TODO Tpetra::ProfileType
#endif

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Indicates that the graph is static, so that new entries cannot be added to this matrix. */
    virtual bool isStaticGraph() const =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Extract a list of entries in a specified global row of this matrix. Put into pre-allocated storage.
    /*!
      \param LocalRow - (In) Global row number for which indices are desired.
      \param Indices - (Out) Global column indices corresponding to values.
      \param Values - (Out) Matrix values.
      \param NumEntries - (Out) Number of indices.

      Note: A std::runtime_error exception is thrown if either \c Indices or \c Values is not large enough to hold the data associated
      with row \c GlobalRow. If \c GlobalRow does not belong to this node, then \c Indices and \c Values are unchanged and \c NumIndices is 
      returned as OrdinalTraits<size_t>::invalid().
    */
    virtual void getGlobalRowCopy(GlobalOrdinal GlobalRow,
                                  const ArrayView<GlobalOrdinal> &Indices,
                                  const ArrayView<Scalar> &Values,
                                  size_t &NumEntries
                                  ) const =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //! Extract a list of entries in a specified local row of the matrix. Put into storage allocated by calling routine.
    /*!
      \param LocalRow - (In) Local row number for which indices are desired.
      \param Indices - (Out) Local column indices corresponding to values.
      \param Values - (Out) Matrix values.
      \param NumIndices - (Out) Number of indices.

      Note: A std::runtime_error exception is thrown if either \c Indices or \c Values is not large enough to hold the data associated
      with row \c LocalRow. If \c LocalRow is not valid for this node, then \c Indices and \c Values are unchanged and \c NumIndices is 
      returned as OrdinalTraits<size_t>::invalid().

      \pre <tt>isLocallyIndexed()==true</tt> or <tt>hasColMap() == true</tt>
    */
    virtual void getLocalRowCopy(LocalOrdinal LocalRow, 
                                 const ArrayView<LocalOrdinal> &Indices, 
                                 const ArrayView<Scalar> &Values,
                                 size_t &NumEntries
                                 ) const =0;

    //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
    /*!
      \param GlobalRow - (In) Global row number for which indices are desired.
      \param Indices   - (Out) Global column indices corresponding to values.
      \param Values    - (Out) Row values
      \pre <tt>isLocallyIndexed() == false</tt>
      \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

      Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
    */
    virtual void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &indices, ArrayView<const Scalar> &values) const =0;

    //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
    /*!
      \param LocalRow - (In) Local row number for which indices are desired.
      \param Indices  - (Out) Global column indices corresponding to values.
      \param Values   - (Out) Row values
      \pre <tt>isGloballyIndexed() == false</tt>
      \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

      Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
    */
    virtual void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices, ArrayView<const Scalar> &values) const =0;

    //! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
    /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the 
      the zero and non-zero diagonals owned by this node. */
    virtual void getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const =0;

    //@}

    //! @name Advanced Matrix-vector multiplication and solve methods
    //@{

    //! Multiplies this matrix by a MultiVector.
    /*! \c X is required to be post-imported, i.e., described by the column map of the matrix. \c Y is required to be pre-exported, i.e., described by the row map of the matrix.

    Both are required to have constant stride, and they are not permitted to ocupy overlapping space. No runtime checking will be performed in a non-debug build.

    This method is templated on the scalar type of MultiVector objects, allowing this method to be applied to MultiVector objects of arbitrary type. However, it is recommended that multiply() not be called directly; instead, use the CrsMatrixMultiplyOp, as it will handle the import/exprt operations required to apply a matrix with non-trivial communication needs.
          
    If \c beta is equal to zero, the operation will enjoy overwrite semantics (\c Y will be overwritten with the result of the multiplication). Otherwise, the result of the multiplication
    will be accumulated into \c Y.
    */
    //TODO virtual=0 // TODO: Add default parameters ?
#ifdef CTHULHU_NOT_IMPLEMENTED
    template <class DomainScalar, class RangeScalar>
    void multiply(const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans, RangeScalar alpha, RangeScalar beta) const;
#endif // CTHULHU_NOT_IMPLEMENTED

    //    virtual void multiply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans, Scalar alpha, Scalar beta) const=0;

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Solves a linear system when the underlying matrix is triangular.
    /*! \c X is required to be post-imported, i.e., described by the column map of the matrix. \c Y is required to be pre-exported, i.e., described by the row map of the matrix.

    This method is templated on the scalar type of MultiVector objects, allowing this method to be applied to MultiVector objects of arbitrary type. However, it is recommended that solve() not be called directly; instead, use the CrsMatrixSolveOp, as it will handle the import/exprt operations required to apply a matrix with non-trivial communication needs.
          
    Both are required to have constant stride. However, unlike multiply(), it is permissible for <tt>&X == &Y</tt>. No runtime checking will be performed in a non-debug build.
    */
    template <class DomainScalar, class RangeScalar>
    void solve(const MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> & Y, MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X, Teuchos::ETransp trans) const;
#endif // CTHULHU_NOT_IMPLEMENTED
          
    //@}

    //! @name Methods implementing Operator
    //@{ 

    //! \brief Computes the sparse matrix-multivector multiplication.
    /*! Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exceptions:
      - if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
    */
    virtual void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                       Teuchos::ETransp mode = Teuchos::NO_TRANS,
                       Scalar alpha = ScalarTraits<Scalar>::one(),
                       Scalar beta = ScalarTraits<Scalar>::zero()) const =0;

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Indicates whether this operator supports applying the adjoint operator.
    virtual bool hasTransposeApply() const =0;
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //! \brief Returns the Map associated with the domain of this operator.
    //! This will be <tt>null</tt> until fillComplete() is called.
    virtual const RCP<const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const =0;

    //! Returns the Map associated with the domain of this operator.
    //! This will be <tt>null</tt> until fillComplete() is called.
    virtual const RCP<const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const =0;

    //@}

    //! @name Overridden from Teuchos::Describable 
    //@{

    /** \brief Return a simple one-line description of this object. */
    virtual std::string description() const =0;

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const =0;

    //@}

    //! @name Methods implementing Cthulhu::DistObject
    //@{
#ifdef CTHULHU_NOT_IMPLEMENTED
    // TODO: DistObject
    virtual bool checkSizes(const DistObject<char, LocalOrdinal,GlobalOrdinal,Node>& source) =0;
#endif // CTHULHU_NOT_IMPLEMENTED

#ifdef CTHULHU_NOT_IMPLEMENTED
    virtual void copyAndPermute(const DistObject<char, LocalOrdinal,GlobalOrdinal,Node>& source,
                                size_t numSameIDs,
                                const ArrayView<const LocalOrdinal> &permuteToLIDs,
                                const ArrayView<const LocalOrdinal> &permuteFromLIDs) =0;
#endif // CTHULHU_NOT_IMPLEMENTED

#ifdef CTHULHU_NOT_IMPLEMENTED

    virtual void packAndPrepare(const DistObject<char, LocalOrdinal,GlobalOrdinal,Node>& source,
                                const ArrayView<const LocalOrdinal> &exportLIDs,
                                Array<char> &exports,
                                const ArrayView<size_t> & numPacketsPerLID,
                                size_t& constantNumPackets,
                                Distributor &distor) =0;
#endif // CTHULHU_NOT_IMPLEMENTED

#ifdef CTHULHU_NOT_IMPLEMENTED
    virtual void unpackAndCombine(const ArrayView<const LocalOrdinal> &importLIDs,
                                  const ArrayView<const char> &imports,
                                  const ArrayView<size_t> &numPacketsPerLID,
                                  size_t constantNumPackets,
                                  Distributor &distor,
                                  CombineMode CM) =0;
#endif // CTHULHU_NOT_IMPLEMENTED
    //@}

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! \name Deprecated routines to be removed at some point in the future.
    //@{

    /** \brief Deprecated. Re-allocate the data into contiguous storage.

    This method is deprecated and will be removed in a future version of Cthulhu, as 
    the implementation of storage optimization has been below Cthulhu to Kokkos.

    Currently, the implementation simply calls resumeFill() and then fillComplete(OptimizeStorage). As such, it is 
    required to be called by all nodes that participate in the associated communicator.
    */
    virtual CTHULHU_DEPRECATED void optimizeStorage();

   //! Deprecated. Get a persisting const view of the entries in a specified global row of this matrix.
    virtual CTHULHU_DEPRECATED void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayRCP<const GlobalOrdinal> &indices, ArrayRCP<const Scalar> &values) const;

   //! Deprecated. Get a persisting const view of the entries in a specified local row of this matrix.
    virtual CTHULHU_DEPRECATED void getLocalRowView(LocalOrdinal LocalRow, ArrayRCP<const LocalOrdinal> &indices, ArrayRCP<const Scalar> &values) const;

    //@}
#endif //CTHULHU_NOT_IMPLEMENTED

  }; // class CrsMatrix

} // namespace Cthulhu

#define CTHULHU_CRSMATRIX_SHORT
#endif
