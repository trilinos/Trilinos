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

    TpetraCrsMatrix(const Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > &mtx) : mtx_(mtx) {  }

    // !Destructor.
    virtual ~TpetraCrsMatrix() { }

     //@}

     //! @name Insertion/Removal Methods
     //@{ 

     //! Insert matrix entries, using global IDs.
    inline void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &cols, const ArrayView<const Scalar> &vals) {  mtx_->insertGlobalValues(globalRow, cols, vals); }

     //! Insert matrix entries, using local IDs.
     inline void insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &cols, const ArrayView<const Scalar> &vals) {  mtx_->insertLocalValues(localRow, cols, vals); }

     //! \brief Replace matrix entries, using global IDs.
     inline void replaceGlobalValues(GlobalOrdinal globalRow, 
                                     const ArrayView<const GlobalOrdinal> &cols,
                                     const ArrayView<const Scalar>        &vals) {  mtx_->replaceGlobalValues(globalRow, cols, vals); }

     //! Replace matrix entries, using local IDs.
     inline void replaceLocalValues(LocalOrdinal localRow, 
                                    const ArrayView<const LocalOrdinal> &cols,
                                    const ArrayView<const Scalar>       &vals) {  mtx_->replaceLocalValues(localRow, cols, vals); }

     //! Sum into multiple entries, using global IDs.
     inline void sumIntoGlobalValues(GlobalOrdinal globalRow, 
                                     const ArrayView<const GlobalOrdinal> &cols,
                                     const ArrayView<const Scalar>        &vals) {  mtx_->sumIntoGlobalValues(globalRow, cols, vals); }


     //! Sum into multiple entries, using local IDs.
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

    //! Resume fill operations.
    inline void resumeFill() {  mtx_->resumeFill(); }

    //! Signal that data entry is complete, specifying domain and range maps.
    inline void fillComplete(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap, OptimizeOption os = DoOptimizeStorage) { 
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, domainMap, tDomainMap, "Cthulhu::TpetraCrsMatrix::fillComplete() only accept Cthulhu::TpetraMap as input arguments.");
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, rangeMap,  tRangeMap,  "Cthulhu::TpetraCrsMatrix:fillComplete() only accept Cthulhu::TpetraMap as input arguments.");
      mtx_->fillComplete(tDomainMap->getTpetra_Map(), tRangeMap->getTpetra_Map()); // TODO: os 
    }

    //! Signal that data entry is complete. 
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

     inline bool isStorageOptimized() const {  return mtx_->isStorageOptimized(); }

     //! Returns \c true if the matrix was allocated with static data structures.
    inline Cthulhu::ProfileType getProfileType() const {  return mtx_->getProfileType(); } // TODO Tpetra::ProfileType

     //! Indicates that the graph is static, so that new entries cannot be added to this matrix. */
     inline bool isStaticGraph() const {  return mtx_->isStaticGraph(); }

     //! Extract a list of entries in a specified global row of this matrix. Put into pre-allocated storage.
      inline void getGlobalRowCopy(GlobalOrdinal GlobalRow,
                                  const ArrayView<GlobalOrdinal> &Indices,
                                  const ArrayView<Scalar> &Values,
                                  size_t &NumEntries
                                  ) const {  mtx_->getGlobalRowCopy(GlobalRow, Indices, Values, NumEntries); }

     //! Extract a list of entries in a specified local row of the matrix. Put into storage allocated by calling routine.
      inline void getLocalRowCopy(LocalOrdinal LocalRow, 
                                 const ArrayView<LocalOrdinal> &Indices, 
                                 const ArrayView<Scalar> &Values,
                                 size_t &NumEntries
                                 ) const {  mtx_->getLocalRowCopy(LocalRow, Indices, Values, NumEntries); }

     //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
      inline void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &indices, ArrayView<const Scalar> &values) const {  mtx_->getGlobalRowView(GlobalRow, indices, values); }

     //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
      inline void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices, ArrayView<const Scalar> &values) const {  mtx_->getLocalRowView(LocalRow, indices, values); }

     //! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
     /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the 
       the zero and non-zero diagonals owned by this node. */
    inline void getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const { 
       
      
      CTHULHU_DYNAMIC_CAST(TpetraVectorClass, diag, tDiag, "Cthulhu::TpetraCrsMatrix.getLocalDiagCopy() only accept Cthulhu::TpetraVector as input arguments.");
      mtx_->getLocalDiagCopy(*tDiag.getTpetra_Vector()); 
    }
    
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
