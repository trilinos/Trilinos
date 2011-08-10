#ifndef CTHULHU_TPETRACRSMATRIX_HPP
#define CTHULHU_TPETRACRSMATRIX_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifndef HAVE_CTHULHU_TPETRA
#error This file should be included only if HAVE_CTHULHU_TPETRA is defined.
#endif

#include "Tpetra_CrsMatrix.hpp"

#include "Cthulhu_CrsMatrix.hpp"

#include "Cthulhu_TpetraMap.hpp"
#include "Cthulhu_TpetraMultiVector.hpp"
#include "Cthulhu_TpetraVector.hpp"
#include "Cthulhu_TpetraCrsGraph.hpp"

#include "Cthulhu_Exceptions.hpp"

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
            class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps> 
  class TpetraCrsMatrix : public CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> {
    
    // The following typedef are used by the CTHULHU_DYNAMIC_CAST() macro.
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMapClass;
    typedef TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TpetraMultiVectorClass;
    typedef TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TpetraVectorClass;
    typedef TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> TpetraCrsMatrixClass;
    typedef TpetraImport<LocalOrdinal, GlobalOrdinal, Node> TpetraImportClass;
    typedef TpetraExport<LocalOrdinal, GlobalOrdinal, Node> TpetraExportClass;    

  public:
    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructor specifying the number of non-zeros for all rows.
    TpetraCrsMatrix(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, size_t maxNumEntriesPerRow, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile) 
    { 
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, rowMap, tRowMap, "Cthulhu::TpetraCrsMatrix constructors only accept Cthulhu::TpetraMap as input arguments.");
      mtx_ = rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps>(tRowMap->getTpetra_Map(), maxNumEntriesPerRow)); //TODO: convert, pftype));
    }

    //! Constructor specifying the number of non-zeros for each row.
    TpetraCrsMatrix(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile)
    { 
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, rowMap, tRowMap, "Cthulhu::TpetraCrsMatrix constructors only accept Cthulhu::TpetraMap as input arguments.");
      mtx_ = rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>(tRowMap->getTpetra_Map(), NumEntriesPerRowToAlloc)); // TODO convert, pftype));
    }

    //! Constructor specifying a column map and the number of non-zeros for all rows.
    /** The column map will be used to filter any matrix entries inserted using insertLocalValues() or insertGlobalValues().
     */
    TpetraCrsMatrix(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, size_t maxNumEntriesPerRow, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile)
    { 
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, rowMap, tRowMap, "Cthulhu::TpetraCrsMatrix constructors only accept Cthulhu::TpetraMap as input arguments.");
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, colMap, tColMap, "Cthulhu::TpetraCrsMatrix constructors only accept Cthulhu::TpetraMap as input arguments.");
      mtx_ = rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>(tRowMap->getTpetra_Map(), tColMap->getTpetra_Map(), maxNumEntriesPerRow));// TODO convert, pftype));
    }
    
    //! Constructor specifying a column map and the number of non-zeros for each row.
    /** The column map will be used to filter any matrix entries inserted using insertLocalValues() or insertGlobalValues().
     */
    TpetraCrsMatrix(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile)
    { 
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, rowMap, tRowMap, "Cthulhu::TpetraCrsMatrix constructors only accept Cthulhu::TpetraMap as input arguments.");
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, colMap, tColMap, "Cthulhu::TpetraCrsMatrix constructors only accept Cthulhu::TpetraMap as input arguments.");
      mtx_ = rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>(tRowMap->getTpetra_Map(), tColMap->getTpetra_Map(), NumEntriesPerRowToAlloc));// TODO convert, pftype));
    } 

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Constructor specifying a pre-constructed graph.
    // TODO: need a CrsGraph
    explicit TpetraCrsMatrix(const RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &graph) {  }
#endif // CTHULHU_NOT_IMPLEMENTED

    TpetraCrsMatrix(const Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > &mtx) : mtx_(mtx) {  }

    // !Destructor.
    virtual ~TpetraCrsMatrix() { }

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
    inline void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &cols, const ArrayView<const Scalar> &vals) {  mtx_->insertGlobalValues(globalRow, cols, vals); }

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
     inline void insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &cols, const ArrayView<const Scalar> &vals) {  mtx_->insertLocalValues(localRow, cols, vals); }

     //! \brief Replace matrix entries, using global IDs.
     /** All index values must be in the global space. 

     \pre \c globalRow is a global row belonging to the matrix on this node.

     \note If (globalRow,cols[i]) corresponds to an entry that is duplicated in this matrix row (likely because it was inserted more than once and fillComplete() has not been called in the interim), the behavior of this function is not defined. */
     inline void replaceGlobalValues(GlobalOrdinal globalRow, 
                                     const ArrayView<const GlobalOrdinal> &cols,
                                     const ArrayView<const Scalar>        &vals) {  mtx_->replaceGlobalValues(globalRow, cols, vals); }

     //! Replace matrix entries, using local IDs.
     /** All index values must be in the local space. 
      */
     inline void replaceLocalValues(LocalOrdinal localRow, 
                                    const ArrayView<const LocalOrdinal> &cols,
                                    const ArrayView<const Scalar>       &vals) {  mtx_->replaceLocalValues(localRow, cols, vals); }

     //! Sum into multiple entries, using global IDs.
     /** All index values must be in the global space. 

     \pre \c globalRow is a global row belonging to the matrix on this node.

     */
     inline void sumIntoGlobalValues(GlobalOrdinal globalRow, 
                                     const ArrayView<const GlobalOrdinal> &cols,
                                     const ArrayView<const Scalar>        &vals) {  mtx_->sumIntoGlobalValues(globalRow, cols, vals); }


     //! Sum into multiple entries, using local IDs.
     /** All index values must be in the local space. 

     \pre \c localRow is a local row belonging to the matrix on this node.

     */
     inline void sumIntoLocalValues(LocalOrdinal globalRow, 
                                    const ArrayView<const LocalOrdinal>  &cols,
                                    const ArrayView<const Scalar>        &vals) {  mtx_->sumIntoLocalValues(globalRow, cols, vals); } 

     //! Set all matrix entries equal to scalarThis.
     inline void setAllToScalar(const Scalar &alpha) {  mtx_->setAllToScalar(alpha); }

     //! Scale the current values of a matrix, this = alpha*this. 
     inline void scale(const Scalar &alpha) {  mtx_->scale(alpha); }

     //@}

     //! @name Transformational Methods
     //@{ 

     //! \brief Communicate non-local contributions to other nodes.
     inline void globalAssemble() {  mtx_->globalAssemble(); }

     /*! Resume fill operations.
       After calling fillComplete(), resumeFill() must be called before initiating any changes to the matrix.

       resumeFill() may be called repeatedly. 

       \post  <tt>isFillActive() == true<tt>
       \post  <tt>isFillComplete() == false<tt>
     */
     inline void resumeFill() {  mtx_->resumeFill(); }

     /*! \brief Signal that data entry is complete, specifying domain and range maps.

     Off-node indices are distributed (via globalAssemble()), indices are sorted, redundant indices are eliminated, and global indices are transformed to local indices.

     \pre  <tt>isFillActive() == true<tt>
     \pre <tt>isFillComplete()() == false<tt>

     \post <tt>isFillActive() == false<tt>
     \post <tt>isFillComplete() == true<tt>
     \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
     */ 
    inline void fillComplete(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap, OptimizeOption os = DoOptimizeStorage) { 
       
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, domainMap, tDomainMap, "Cthulhu::TpetraCrsMatrix::fillComplete() only accept Cthulhu::TpetraMap as input arguments.");
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, rangeMap,  tRangeMap,  "Cthulhu::TpetraCrsMatrix:fillComplete() only accept Cthulhu::TpetraMap as input arguments.");
      mtx_->fillComplete(tDomainMap->getTpetra_Map(), tRangeMap->getTpetra_Map()); // TODO: os 
    }

     /*! \brief Signal that data entry is complete. 

     Off-node entries are distributed (via globalAssemble()), repeated entries are summed, and global indices are transformed to local indices.

     \note This method calls fillComplete( getRowMap(), getRowMap(), os ).

     \pre  <tt>isFillActive() == true<tt>
     \pre <tt>isFillComplete()() == false<tt>

     \post <tt>isFillActive() == false<tt>
     \post <tt>isFillComplete() == true<tt>
     \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
     */
    // TODO: Tpetra::OptimizeOption
    inline void fillComplete(Cthulhu::OptimizeOption os = Cthulhu::DoOptimizeStorage) { 
       
      if (os == Cthulhu::DoOptimizeStorage)
        mtx_->fillComplete(Tpetra::DoOptimizeStorage); 
      else if (os == Cthulhu::DoNotOptimizeStorage)
        mtx_->fillComplete(Tpetra::DoNotOptimizeStorage); 
      else TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::BadCast, "Cannot convert 'Cthulhu::OptimizeOption' to a 'Tpetra::OptimizeOption'. Cthulhu::OptimizeOption 'os' have a bad value.");
    }


    //! leftScale
    inline void leftScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x){
       
      CTHULHU_DYNAMIC_CAST(const TpetraVectorClass, x, tX, "Cthulhu::TpetraCrsMatrix->leftScale() only accept Cthulhu::TpetraVector as input arguments.");
      mtx_->leftScale(tX);
    }


     //@}

     //! @name Methods implementing RowMatrix
     //@{ 

     //! Returns the communicator.
    inline const RCP<const Comm<int> > getComm() const {  return mtx_->getComm(); } // removed &

     //! Returns the underlying node.
     inline RCP<Node> getNode() const {  return mtx_->getNode(); }

     //! Returns the Map that describes the row distribution in this matrix.
    inline const RCP<const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> > getRowMap() const { 
       
      return rcp( new TpetraMapClass(mtx_->getRowMap()) );
    }
     
     //! \brief Returns the Map that describes the column distribution in this matrix.
    inline const RCP<const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> > getColMap() const { 
       
      return rcp( new TpetraMapClass(mtx_->getColMap()) );
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
     //! Returns the RowGraph associated with this matrix. 
    inline RCP<const RowGraph<LocalOrdinal,GlobalOrdinal,Node> > getGraph() const {  return mtx_->getGraph(); } // wrapped by a Cthulhu object
#endif // CTHULHU_NOT_IMPLEMENTED

     //! Returns the CrsGraph associated with this matrix. 
    inline RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > getCrsGraph() const { 
       
      
      RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > graph = Teuchos::rcp_const_cast<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >(mtx_->getCrsGraph()); //TODO: can I avoid the const_cast ?
      return rcp ( new Cthulhu::TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>(graph) );
    }

     //! Returns the number of global rows in this matrix.
     /** Undefined if isFillActive().
      */
     inline global_size_t getGlobalNumRows() const {  return mtx_->getGlobalNumRows(); }

     //! \brief Returns the number of global columns in the matrix.
     /** Undefined if isFillActive().
      */
     inline global_size_t getGlobalNumCols() const {  return mtx_->getGlobalNumCols(); }

     //! Returns the number of matrix rows owned on the calling node.
     inline size_t getNodeNumRows() const {  return mtx_->getNodeNumRows(); }

     //! Returns the number of columns connected to the locally owned rows of this matrix.
     /** Throws std::runtime_error if <tt>hasColMap() == false</tt>
      */
     inline size_t getNodeNumCols() const {  return mtx_->getNodeNumCols(); }

     //! Returns the index base for global indices for this matrix. 
     inline GlobalOrdinal getIndexBase() const {  return mtx_->getIndexBase(); }

     //! Returns the global number of entries in this matrix.
     inline global_size_t getGlobalNumEntries() const {  return mtx_->getGlobalNumEntries(); }

     //! Returns the local number of entries in this matrix.
     inline size_t getNodeNumEntries() const {  return mtx_->getNodeNumEntries(); }

     //! \brief Returns the current number of entries on this node in the specified global row.
     /*! Returns OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this matrix. */
     inline size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const {  return mtx_->getNumEntriesInGlobalRow(globalRow); }

     //! Returns the current number of entries on this node in the specified local row.
     /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
     inline size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const {  return mtx_->getNumEntriesInLocalRow(localRow); }

     //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
     /** Undefined if isFillActive().
      */
     inline global_size_t getGlobalNumDiags() const {  return mtx_->getGlobalNumDiags(); }

     //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
     /** Undefined if isFillActive().
      */
     inline size_t getNodeNumDiags() const {  return mtx_->getNodeNumDiags(); }

     //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
     /** Undefined if isFillActive().
      */
     inline size_t getGlobalMaxNumRowEntries() const {  return mtx_->getGlobalMaxNumRowEntries(); }

     //! \brief Returns the maximum number of entries across all rows/columns on this node.
     /** Undefined if isFillActive().
      */
     inline size_t getNodeMaxNumRowEntries() const {  return mtx_->getNodeMaxNumRowEntries(); }

     //! \brief Indicates whether the matrix has a well-defined column map. 
     inline bool hasColMap() const {  return mtx_->hasColMap(); } 

     //! \brief Indicates whether the matrix is lower triangular.
     /** Undefined if isFillActive().
      */
     inline bool isLowerTriangular() const {  return mtx_->isLowerTriangular(); }

     //! \brief Indicates whether the matrix is upper triangular.
     /** Undefined if isFillActive().
      */
     inline bool isUpperTriangular() const {  return mtx_->isUpperTriangular(); }

     //! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
     inline bool isLocallyIndexed() const {  return mtx_->isLocallyIndexed(); }

     //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
     inline bool isGloballyIndexed() const {  return mtx_->isGloballyIndexed(); }

     //! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
     inline bool isFillComplete() const {  return mtx_->isFillComplete(); }

     //! Returns \c true if resumeFill() has been called and the matrix is in edit mode.
     inline bool isFillActive() const {  return mtx_->isFillActive(); }

     //! \brief Returns \c true if storage has been optimized.
     /**
        Optimized storage means that the allocation of each row is equal to the
        number of entries. The effect is that a pass through the matrix, i.e.,
        during a mat-vec, requires minimal memory traffic. One limitation of
        optimized storage is that no new indices can be added to the matrix.
     */
     inline bool isStorageOptimized() const {  return mtx_->isStorageOptimized(); }

     //! Returns \c true if the matrix was allocated with static data structures.
    inline Cthulhu::ProfileType getProfileType() const {  return mtx_->getProfileType(); } // TODO Tpetra::ProfileType

     //! Indicates that the graph is static, so that new entries cannot be added to this matrix. */
     inline bool isStaticGraph() const {  return mtx_->isStaticGraph(); }

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
     inline void getGlobalRowCopy(GlobalOrdinal GlobalRow,
                                  const ArrayView<GlobalOrdinal> &Indices,
                                  const ArrayView<Scalar> &Values,
                                  size_t &NumEntries
                                  ) const {  mtx_->getGlobalRowCopy(GlobalRow, Indices, Values, NumEntries); }

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
     inline void getLocalRowCopy(LocalOrdinal LocalRow, 
                                 const ArrayView<LocalOrdinal> &Indices, 
                                 const ArrayView<Scalar> &Values,
                                 size_t &NumEntries
                                 ) const {  mtx_->getLocalRowCopy(LocalRow, Indices, Values, NumEntries); }

     //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
     /*!
       \param GlobalRow - (In) Global row number for which indices are desired.
       \param Indices   - (Out) Global column indices corresponding to values.
       \param Values    - (Out) Row values
       \pre <tt>isLocallyIndexed() == false</tt>
       \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

       Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
     */
     inline void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &indices, ArrayView<const Scalar> &values) const {  mtx_->getGlobalRowView(GlobalRow, indices, values); }

     //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
     /*!
       \param LocalRow - (In) Local row number for which indices are desired.
       \param Indices  - (Out) Global column indices corresponding to values.
       \param Values   - (Out) Row values
       \pre <tt>isGloballyIndexed() == false</tt>
       \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

       Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
     */
     inline void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices, ArrayView<const Scalar> &values) const {  mtx_->getLocalRowView(LocalRow, indices, values); }

     //! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
     /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the 
       the zero and non-zero diagonals owned by this node. */
    inline void getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const { 
       
      
      CTHULHU_DYNAMIC_CAST(TpetraVectorClass, diag, tDiag, "Cthulhu::TpetraCrsMatrix.getLocalDiagCopy() only accept Cthulhu::TpetraVector as input arguments.");
      mtx_->getLocalDiagCopy(*tDiag.getTpetra_Vector()); 
    }
    
     //@}

     //! @name Advanced Matrix-vector multiplication and solve methods
     //@{

#ifdef CTHULHU_NOT_IMPLEMENTED
     //! Multiplies this matrix by a MultiVector.
     /*! \c X is required to be post-imported, i.e., described by the column map of the matrix. \c Y is required to be pre-exported, i.e., described by the row map of the matrix.

     Both are required to have constant stride, and they are not permitted to ocupy overlapping space. No runtime checking will be performed in a non-debug build.

     This method is templated on the scalar type of MultiVector objects, allowing this method to be applied to MultiVector objects of arbitrary type. However, it is recommended that multiply() not be called directly; instead, use the TpetraCrsMatrixMultiplyOp, as it will handle the import/exprt operations required to apply a matrix with non-trivial communication needs.
          
     If \c beta is equal to zero, the operation will enjoy overwrite semantics (\c Y will be overwritten with the result of the multiplication). Otherwise, the result of the multiplication
     will be accumulated into \c Y.
     */
    template <class DomainScalar, class RangeScalar>
    inline void multiply(const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans, RangeScalar alpha, RangeScalar beta) const {  mtx_->multiply(X, Y, trans, alpha, beta); }
#endif // CTHULHU_NOT_IMPLEMENTED

    // TODO Note: Do we need to use a Tpetra::CrsMatrixMultiplyOp ?? 
    //            (Tpetra Doc of multiply: it is recommended that multiply() not be called directly; instead, use the CrsMatrixMultiplyOp, as it will handle the import/exprt operations required to apply a matrix with non-trivial communication needs.)

    // TODO : templated type

//     inline void multiply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans, Scalar alpha, Scalar beta) const {
//        

//       CTHULHU_DYNAMIC_CAST(const TpetraMultiVectorClass, X, tX, "Cthulhu::TpetraCrsMatrix->multiply() only accept Cthulhu::TpetraMultiVector as input arguments.");
//       CTHULHU_DYNAMIC_CAST(      TpetraMultiVectorClass, Y, tY, "Cthulhu::TpetraCrsMatrix->multiply() only accept Cthulhu::TpetraMultiVector as input arguments.");
//       mtx_->template multiply<Scalar,Scalar>(*tX.getTpetra_MultiVector(), *tY.getTpetra_MultiVector(), trans, alpha, beta);
//     }

#ifdef CTHULHU_NOT_IMPLEMENTED
     //! Solves a linear system when the underlying matrix is triangular.
     /*! \c X is required to be post-imported, i.e., described by the column map of the matrix. \c Y is required to be pre-exported, i.e., described by the row map of the matrix.

     This method is templated on the scalar type of MultiVector objects, allowing this method to be applied to MultiVector objects of arbitrary type. However, it is recommended that solve() not be called directly; instead, use the TpetraCrsMatrixSolveOp, as it will handle the import/exprt operations required to apply a matrix with non-trivial communication needs.
          
     Both are required to have constant stride. However, unlike multiply(), it is permissible for <tt>&X == &Y</tt>. No runtime checking will be performed in a non-debug build.
     */
     template <class DomainScalar, class RangeScalar>
     inline void solve(const MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> & Y, MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X, Teuchos::ETransp trans) const {  mtx_->solve(Y, X, trans); }
#endif // CTHULHU_NOT_IMPLEMENTED          
     //@}

     //! @name Methods implementing Operator
     //@{ 

     //! \brief Computes the sparse matrix-multivector multiplication.
     /*! Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exceptions:
       - if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
     */
     inline void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                       Teuchos::ETransp mode = Teuchos::NO_TRANS,
                       Scalar alpha = ScalarTraits<Scalar>::one(),
                       Scalar beta = ScalarTraits<Scalar>::zero()) const { 
        

      CTHULHU_DYNAMIC_CAST(const TpetraMultiVectorClass, X, tX, "Cthulhu::TpetraCrsMatrix->multiply() only accept Cthulhu::TpetraMultiVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(      TpetraMultiVectorClass, Y, tY, "Cthulhu::TpetraCrsMatrix->multiply() only accept Cthulhu::TpetraMultiVector as input arguments.");
      mtx_->apply(*tX.getTpetra_MultiVector(), *tY.getTpetra_MultiVector(), mode, alpha, beta);
     }

     //! Indicates whether this operator supports applying the adjoint operator.
     inline bool hasTransposeApply() const {  return mtx_->hasTransposeApply(); }

     //! \brief Returns the Map associated with the domain of this operator.
     //! This will be <tt>null</tt> until fillComplete() is called.
    inline const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const { 
       
      return rcp( new TpetraMapClass(mtx_->getDomainMap()) );
    }
    
     //! Returns the Map associated with the domain of this operator.
     //! This will be <tt>null</tt> until fillComplete() is called.
    inline const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const { 
       
      return rcp( new TpetraMapClass(mtx_->getRangeMap()) );
    }

     //@}

     //! @name Overridden from Teuchos::Describable 
     //@{

     /** \brief Return a simple one-line description of this object. */
    inline std::string description() const {  return mtx_->description(); }
    
     /** \brief Print the object with some verbosity level to an FancyOStream object. */
    inline void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const {  mtx_->describe(out,verbLevel); }
    
    //@}

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! @name Methods implementing Cthulhu::DistObject
    //@{
    inline bool checkSizes(const DistObject<char, LocalOrdinal,GlobalOrdinal,Node>& source) {  return mtx_->checkSizes(); }

 inline void copyAndPermute(const DistObject<char, LocalOrdinal,GlobalOrdinal,Node>& source,
                            size_t numSameIDs,
                            const ArrayView<const LocalOrdinal> &permuteToLIDs,
                            const ArrayView<const LocalOrdinal> &permuteFromLIDs) {  mtx_->copyAndPermute(); }

 inline void packAndPrepare(const DistObject<char, LocalOrdinal,GlobalOrdinal,Node>& source,
                            const ArrayView<const LocalOrdinal> &exportLIDs,
                            Array<char> &exports,
                            const ArrayView<size_t> & numPacketsPerLID,
                            size_t& constantNumPackets,
                            Distributor &distor) {  mtx_->packAndPrepare(); }

 inline void unpackAndCombine(const ArrayView<const LocalOrdinal> &importLIDs,
                              const ArrayView<const char> &imports,
                              const ArrayView<size_t> &numPacketsPerLID,
                              size_t constantNumPackets,
                              Distributor &distor,
                              CombineMode CM) {  mtx_->unpackAndCombine(); }
    //@}
#endif // CTHULHU_NOT_IMPLEMENTED

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! \name Deprecated routines to be removed at some point in the future.
    //@{

    /** \brief Deprecated. Re-allocate the data into contiguous storage.

    This method is deprecated and will be removed in a future version of Cthulhu, as 
    the implementation of storage optimization has been below Cthulhu to Kokkos.

    Currently, the implementation simply calls resumeFill() and then fillComplete(OptimizeStorage). As such, it is 
    required to be called by all nodes that participate in the associated communicator.
    */
    inline CTHULHU_DEPRECATED void optimizeStorage() {  }
    
    //! Deprecated. Get a persisting const view of the entries in a specified global row of this matrix.
    inline CTHULHU_DEPRECATED void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayRCP<const GlobalOrdinal> &indices, ArrayRCP<const Scalar> &values) const {  } 
    
    //! Deprecated. Get a persisting const view of the entries in a specified local row of this matrix.
    inline CTHULHU_DEPRECATED void getLocalRowView(LocalOrdinal LocalRow, ArrayRCP<const LocalOrdinal> &indices, ArrayRCP<const Scalar> &values) const {  } 
    
    //@}
#endif // CTHULHU_NOT_IMPLEMENTED

    RCP< const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > getTpetra_CrsMatrix() const {  return mtx_; }
    RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > getTpetra_CrsMatrixNonConst() const {  return mtx_; }

    /** TODO : interface of Teuchos_LabeledObject.hpp **/
    void setObjectLabel (const std::string &objectLabel) {  mtx_->setObjectLabel(objectLabel); }

    //{@
    // Implements DistObject interface
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getMap() const { 
       
      return rcp( new TpetraMapClass(mtx_->getMap()) );
    }

    inline void doImport(const DistObject<char, LocalOrdinal,GlobalOrdinal,Node> &source, 
                         const Import<LocalOrdinal,GlobalOrdinal,Node> &importer, CombineMode CM) { 
      
      CTHULHU_DYNAMIC_CAST(const TpetraCrsMatrixClass, source, tSource, "Cthulhu::TpetraCrsMatrix::doImport only accept Cthulhu::TpetraCrsMatrix as input arguments.");
      CTHULHU_DYNAMIC_CAST(const TpetraImportClass, importer, tImporter, "Cthulhu::TpetraImport::doImport only accept Cthulhu::TpetraImport as input arguments.");

      RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal,Node> > v = tSource.getTpetra_CrsMatrix();
      mtx_->doImport(*v, *tImporter.getTpetra_Import(), Cthulhu2Tpetra_CombineMode(CM));
    }

    void doExport(const DistObject<char,LocalOrdinal,GlobalOrdinal,Node> &dest,
                  const Import<LocalOrdinal,GlobalOrdinal,Node>& importer, CombineMode CM) {
      
      
      CTHULHU_DYNAMIC_CAST(const TpetraCrsMatrixClass, dest, tDest, "Cthulhu::TpetraCrsMatrix::doImport only accept Cthulhu::TpetraCrsMatrix as input arguments.");
      CTHULHU_DYNAMIC_CAST(const TpetraImportClass, importer, tImporter, "Cthulhu::TpetraImport::doImport only accept Cthulhu::TpetraImport as input arguments.");

      RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal,Node> > v = tDest.getTpetra_CrsMatrix();
      mtx_->doExport(*v, *tImporter.getTpetra_Import(), Cthulhu2Tpetra_CombineMode(CM)); 

    }

    void doImport(const DistObject<char,LocalOrdinal,GlobalOrdinal,Node> &source,
                  const Export<LocalOrdinal,GlobalOrdinal,Node>& exporter, CombineMode CM) {
      
      CTHULHU_DYNAMIC_CAST(const TpetraCrsMatrixClass, source, tSource, "Cthulhu::TpetraCrsMatrix::doImport only accept Cthulhu::TpetraCrsMatrix as input arguments.");
      CTHULHU_DYNAMIC_CAST(const TpetraExportClass, exporter, tExporter, "Cthulhu::TpetraImport::doImport only accept Cthulhu::TpetraImport as input arguments.");

      RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal,Node> > v = tSource.getTpetra_CrsMatrix();
      mtx_->doImport(*v, *tExporter.getTpetra_Export(), Cthulhu2Tpetra_CombineMode(CM));

    }


    void doExport(const DistObject<char,LocalOrdinal,GlobalOrdinal,Node> &dest,
                  const Export<LocalOrdinal,GlobalOrdinal,Node>& exporter, CombineMode CM) {
      
      
      CTHULHU_DYNAMIC_CAST(const TpetraCrsMatrixClass, dest, tDest, "Cthulhu::TpetraCrsMatrix::doImport only accept Cthulhu::TpetraCrsMatrix as input arguments.");
      CTHULHU_DYNAMIC_CAST(const TpetraExportClass, exporter, tExporter, "Cthulhu::TpetraImport::doImport only accept Cthulhu::TpetraImport as input arguments.");

      RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal,Node> > v = tDest.getTpetra_CrsMatrix();
      mtx_->doExport(*v, *tExporter.getTpetra_Export(), Cthulhu2Tpetra_CombineMode(CM)); 

    }

    // @}


  private:
    
    RCP< Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > mtx_;

  }; // class TpetraCrsMatrix

} // namespace Cthulhu

#define CTHULHU_TPETRACRSMATRIX_SHORT
#endif
