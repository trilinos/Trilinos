#ifndef CTHULHU_EPETRACRSMATRIX_HPP
#define CTHULHU_EPETRACRSMATRIX_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifndef HAVE_CTHULHU_EPETRA
#error This file should be included only if HAVE_CTHULHU_EPETRA is defined.
#endif

#include <Teuchos_ArrayViewDecl.hpp>

#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>

#include "Cthulhu_CrsMatrix.hpp"
#include "Cthulhu_EpetraCrsGraph.hpp"

#include "Cthulhu_EpetraMap.hpp"
#include "Cthulhu_EpetraMultiVector.hpp"
#include "Cthulhu_EpetraVector.hpp"
#include "Cthulhu_Trans.hpp"

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

  class EpetraCrsMatrix: public CrsMatrix<double,int,int> {
  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructor specifying the number of non-zeros for all rows.
    EpetraCrsMatrix(const RCP<const Map<int,int> > &rowMap, size_t maxNumEntriesPerRow, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile) 
    { CTHULHU_DEBUG_ME;
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, rowMap, eRowMap, "Cthulhu::EpetraCrsMatrix constructors only accept Cthulhu::EpetraMap as input arguments.");

      //TODO: test constblk and blksize=1 for the Map ?
      mtx_ = rcp(new Epetra_CrsMatrix(Copy, eRowMap->getEpetra_Map(), maxNumEntriesPerRow, false)); // TODO Copy or View by default ? // TODO: bool StaticProfile
    }

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Constructor specifying the number of non-zeros for each row.
    EpetraCrsMatrix(const RCP<const Map<int,int> > &rowMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile)
    { CTHULHU_DEBUG_ME;

    }

    //! Constructor specifying a column map and the number of non-zeros for all rows.
    /** The column map will be used to filter any matrix entries inserted using insertLocalValues() or insertGlobalValues().
     */
    EpetraCrsMatrix(const RCP<const Map<int,int> > &rowMap, const RCP<const Map<int,int> > &colMap, size_t maxNumEntriesPerRow, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile)
    { CTHULHU_DEBUG_ME;

    }
    
    //! Constructor specifying a column map and the number of non-zeros for each row.
    /** The column map will be used to filter any matrix entries inserted using insertLocalValues() or insertGlobalValues().
     */
    EpetraCrsMatrix(const RCP<const Map<int,int> > &rowMap, const RCP<const Map<int,int> > &colMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile)
    { CTHULHU_DEBUG_ME;

    } 
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Constructor specifying a pre-constructed graph.
    // TODO: need a CrsGraph
    explicit EpetraCrsMatrix(const RCP<const CrsGraph<int,int,LocalMatOps> > &graph) { CTHULHU_DEBUG_ME; }
#endif // CTHULHU_NOT_IMPLEMENTED

    EpetraCrsMatrix(const Teuchos::RCP<Epetra_CrsMatrix> &mtx) : mtx_(mtx) { CTHULHU_DEBUG_ME; }

    // !Destructor.
    virtual ~EpetraCrsMatrix() { CTHULHU_DEBUG_ME; CTHULHU_DEBUG_ME_PRINT; }

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
    inline void insertGlobalValues(int globalRow, const ArrayView<const int> &cols, const ArrayView<const double> &vals) { 
      CTHULHU_DEBUG_ME; 
      CTHULHU_ERR_CHECK(mtx_->InsertGlobalValues(globalRow, vals.size(), vals.getRawPtr(), cols.getRawPtr())); 
    }

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
    inline void insertLocalValues(int localRow, const ArrayView<const int> &cols, const ArrayView<const double> &vals) { CTHULHU_DEBUG_ME; mtx_->insertLocalValues(localRow, cols, vals); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! \brief Replace matrix entries, using global IDs.
    /** All index values must be in the global space. 

    \pre \c globalRow is a global row belonging to the matrix on this node.

    \note If (globalRow,cols[i]) corresponds to an entry that is duplicated in this matrix row (likely because it was inserted more than once and fillComplete() has not been called in the interim), the behavior of this function is not defined. */
    inline void replaceGlobalValues(int globalRow, 
                                    const ArrayView<const int> &cols,
                                    const ArrayView<const double>        &vals) { CTHULHU_DEBUG_ME; mtx_->replaceGlobalValues(globalRow, cols, vals); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Replace matrix entries, using local IDs.
    /** All index values must be in the local space. 
     */
    inline void replaceLocalValues(int localRow, 
                                   const ArrayView<const int> &cols,
                                   const ArrayView<const double>       &vals) { CTHULHU_DEBUG_ME; mtx_->replaceLocalValues(localRow, cols, vals); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Sum into multiple entries, using global IDs.
    /** All index values must be in the global space. 

    \pre \c globalRow is a global row belonging to the matrix on this node.

    */
    inline void sumIntoGlobalValues(int globalRow, 
                                    const ArrayView<const int> &cols,
                                    const ArrayView<const double>        &vals) { CTHULHU_DEBUG_ME; mtx_->sumIntoGlobalValues(globalRow, cols, vals); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Sum into multiple entries, using local IDs.
    /** All index values must be in the local space. 

    \pre \c localRow is a local row belonging to the matrix on this node.

    */
    inline void sumIntoLocalValues(int globalRow, 
                                   const ArrayView<const int>  &cols,
                                   const ArrayView<const double>        &vals) { CTHULHU_DEBUG_ME; mtx_->sumIntoLocalValues(globalRow, cols, vals); } 
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Set all matrix entries equal to scalarThis.
    inline void setAllTodouble(const double &alpha) { CTHULHU_DEBUG_ME; mtx_->setAllTodouble(alpha); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //! Scale the current values of a matrix, this = alpha*this. 
    inline void scale(const double &alpha) { CTHULHU_DEBUG_ME; mtx_->Scale(alpha); }

    //@}

    //! @name Transformational Methods
    //@{ 

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! \brief Communicate non-local contributions to other nodes.
    inline void globalAssemble() { CTHULHU_DEBUG_ME; mtx_->globalAssemble(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    /*! Resume fill operations.
      After calling fillComplete(), resumeFill() must be called before initiating any changes to the matrix.

      resumeFill() may be called repeatedly. 

      \post  <tt>isFillActive() == true<tt>
      \post  <tt>isFillComplete() == false<tt>
    */
    inline void resumeFill() { CTHULHU_DEBUG_ME; mtx_->resumeFill(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

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
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, domainMap, tDomainMap, "Cthulhu::EpetraCrsMatrix::fillComplete() only accept Cthulhu::TpetraMap as input arguments.");
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, rangeMap,  tRangeMap,  "Cthulhu::EpetraCrsMatrix::fillComplete() only accept Cthulhu::TpetraMap as input arguments.");
      mtx_->FillComplete(tDomainMap->getEpetra_Map(), tRangeMap->getEpetra_Map()); // TODO: os 
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
    // TODO: Cthulhu::OptimizeOption
    inline void fillComplete(Cthulhu::OptimizeOption os = Cthulhu::DoOptimizeStorage) { CTHULHU_DEBUG_ME; if (os == Cthulhu::DoOptimizeStorage) mtx_->FillComplete(true); else mtx_->FillComplete(false); }

    //@}

    //! @name Methods implementing RowMatrix
    //@{ 

    //! Returns the communicator.
    inline const RCP<const Comm<int> > getComm() const {
      CTHULHU_DEBUG_ME; 
      
      RCP<const Epetra_Comm> rcpComm = rcpFromRef(mtx_->Comm());
      return Epetra2Teuchos_Comm(rcpComm);
    }

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns the underlying node.
    inline RCP<Node> getNode() const { CTHULHU_DEBUG_ME; return mtx_->getNode(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //! Returns the Map that describes the row distribution in this matrix.
    inline const RCP<const Cthulhu::Map<int,int> > getRowMap() const { 
      CTHULHU_DEBUG_ME; 

      // The Epetra map is owned by the underlying Epetra matrix and freed when the matrix is deleted.
      // I have to make a copy of the map to be sure that the RCP<Map> returned by this method will remain valid even if the matrix is deleted.
      // Note that it is not a problem because the copy constructor of Epetra_Map doesn't really copy the data (there is an reference count mecanism in Epetra).

      RCP<const Epetra_Map> map = rcp(new Epetra_Map(mtx_->RowMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }
     
    //! \brief Returns the Map that describes the column distribution in this matrix.
    inline const RCP<const Cthulhu::Map<int,int> > getColMap() const { 
      CTHULHU_DEBUG_ME; 

      // The Epetra map is owned by the underlying Epetra matrix and freed when the matrix is deleted.
      // I have to make a copy of the map to be sure that the RCP<Map> returned by this method will remain valid even if the matrix is deleted.
      // Note that it is not a problem because the copy constructor of Epetra_Map doesn't really copy the data (there is an reference count mecanism in Epetra).

      RCP<const Epetra_Map> map = rcp(new Epetra_Map(mtx_->ColMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Returns the RowGraph associated with this matrix. 
    inline RCP<const RowGraph<int,int> > getGraph() const { CTHULHU_DEBUG_ME; return mtx_->Graph(); }
#endif // CTHULHU_NOT_IMPLEMENTED

    //! Returns the CrsGraph associated with this matrix. 
    inline RCP<const CrsGraph<int,int> > getCrsGraph() const { CTHULHU_DEBUG_ME; 
      RCP<const Epetra_CrsGraph> const_graph = rcp(new Epetra_CrsGraph(mtx_->Graph()));

      RCP<Epetra_CrsGraph> graph = Teuchos::rcp_const_cast<Epetra_CrsGraph>(const_graph); //TODO: can I avoid the const_cast ?
      return rcp ( new Cthulhu::EpetraCrsGraph(graph) );
    }

    //! Returns the number of global rows in this matrix.
    /** Undefined if isFillActive().
     */
    inline global_size_t getGlobalNumRows() const { CTHULHU_DEBUG_ME; return mtx_->NumGlobalRows(); }

    //! \brief Returns the number of global columns in the matrix.
    /** Undefined if isFillActive().
     */
    inline global_size_t getGlobalNumCols() const { CTHULHU_DEBUG_ME; return mtx_->NumGlobalCols(); }

    //! Returns the number of matrix rows owned on the calling node.
    inline size_t getNodeNumRows() const { CTHULHU_DEBUG_ME; return mtx_->NumMyRows(); }

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns the number of columns connected to the locally owned rows of this matrix.
    /** Throws std::runtime_error if <tt>hasColMap() == false</tt>
     */
    inline size_t getNodeNumCols() const { CTHULHU_DEBUG_ME; return mtx_->getNodeNumCols(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns the index base for global indices for this matrix. 
    inline int getIndexBase() const { CTHULHU_DEBUG_ME; return mtx_->getIndexBase(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //! Returns the global number of entries in this matrix.
    inline global_size_t getGlobalNumEntries() const { CTHULHU_DEBUG_ME; return mtx_->NumGlobalNonzeros(); }

    //! Returns the local number of entries in this matrix.
    inline size_t getNodeNumEntries() const { CTHULHU_DEBUG_ME; return mtx_->NumMyNonzeros(); }

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! \brief Returns the current number of entries on this node in the specified global row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this matrix. */
    inline size_t getNumEntriesInGlobalRow(int globalRow) const { CTHULHU_DEBUG_ME; return mtx_->getNumEntriesInGlobalRow(globalRow); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //! Returns the current number of entries on this node in the specified local row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
    inline size_t getNumEntriesInLocalRow(int localRow) const { CTHULHU_DEBUG_ME; return mtx_->NumMyEntries(localRow); }

    //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
    /** Undefined if isFillActive().
     */
    inline global_size_t getGlobalNumDiags() const { CTHULHU_DEBUG_ME; return mtx_->NumGlobalDiagonals(); }

    //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
    /** Undefined if isFillActive().
     */
    inline size_t getNodeNumDiags() const { CTHULHU_DEBUG_ME; return mtx_->NumMyDiagonals(); }

    //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
    /** Undefined if isFillActive().
     */
    inline size_t getGlobalMaxNumRowEntries() const { CTHULHU_DEBUG_ME; return mtx_->GlobalMaxNumEntries(); }

    //! \brief Returns the maximum number of entries across all rows/columns on this node.
    /** Undefined if isFillActive().
     */
    inline size_t getNodeMaxNumRowEntries() const { CTHULHU_DEBUG_ME; return mtx_->MaxNumEntries(); }

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! \brief Indicates whether the matrix has a well-defined column map. 
    inline bool hasColMap() const { CTHULHU_DEBUG_ME; return mtx_->hasColMap(); } 
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! \brief Indicates whether the matrix is lower triangular.
    /** Undefined if isFillActive().
     */
    inline bool isLowerTriangular() const { CTHULHU_DEBUG_ME; return mtx_->isLowerTriangular(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! \brief Indicates whether the matrix is upper triangular.
    /** Undefined if isFillActive().
     */
    inline bool isUpperTriangular() const { CTHULHU_DEBUG_ME; return mtx_->isUpperTriangular(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
    inline bool isLocallyIndexed() const { CTHULHU_DEBUG_ME; return mtx_->IndicesAreLocal(); }

    //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
    inline bool isGloballyIndexed() const { CTHULHU_DEBUG_ME; return mtx_->IndicesAreGlobal(); }

    //! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
    inline bool isFillComplete() const { CTHULHU_DEBUG_ME; return mtx_->Filled(); }

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns \c true if resumeFill() has been called and the matrix is in edit mode.
    inline bool isFillActive() const { CTHULHU_DEBUG_ME; return mtx_->isFillActive(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! \brief Returns \c true if storage has been optimized.
    /**
       Optimized storage means that the allocation of each row is equal to the
       number of entries. The effect is that a pass through the matrix, i.e.,
       during a mat-vec, requires minimal memory traffic. One limitation of
       optimized storage is that no new indices can be added to the matrix.
    */
    inline bool isStorageOptimized() const { CTHULHU_DEBUG_ME; return mtx_->isStorageOptimized(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns \c true if the matrix was allocated with static data structures.
    inline Cthulhu::ProfileType getProfileType() const { CTHULHU_DEBUG_ME; return mtx_->getProfileType(); } // TODO Cthulhu::ProfileType
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Indicates that the graph is static, so that new entries cannot be added to this matrix. */
    inline bool isStaticGraph() const { CTHULHU_DEBUG_ME; return mtx_->isStaticGraph(); }
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
    inline void getGlobalRowCopy(int GlobalRow,
                                 const ArrayView<int> &Indices,
                                 const ArrayView<double> &Values,
                                 size_t &NumEntries
                                 ) const { CTHULHU_DEBUG_ME; mtx_->getGlobalRowCopy(GlobalRow, Indices, Values, NumEntries); }
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
    //TODO: throw same exception as Tpetra
    inline void getLocalRowCopy(int LocalRow, 
                                const ArrayView<int> &Indices, 
                                const ArrayView<double> &Values,
                                size_t &NumEntries
                                ) const { 
      CTHULHU_DEBUG_ME; 
      
      int numEntries=-1;
      CTHULHU_ERR_CHECK(mtx_->ExtractMyRowCopy(LocalRow, Indices.size(), numEntries, Values.getRawPtr(), Indices.getRawPtr()));
      NumEntries = numEntries;
    }

    //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
    /*!
      \param GlobalRow - (In) Global row number for which indices are desired.
      \param Indices   - (Out) Global column indices corresponding to values.
      \param Values    - (Out) Row values
      \pre <tt>isLocallyIndexed() == false</tt>
      \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

      Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
    */
    inline void getGlobalRowView(int GlobalRow, ArrayView<const int> &indices, ArrayView<const double> &values) const { 
      CTHULHU_DEBUG_ME; 
    
      int      numEntries;
      double * eValues;
      int    * eIndices;
      
      CTHULHU_ERR_CHECK(mtx_->ExtractGlobalRowView(GlobalRow, numEntries, eValues, eIndices));
      if (numEntries == 0) { eValues = NULL; eIndices = NULL; } // Cf. TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

      indices = ArrayView<const int>(eIndices, numEntries);
      values  = ArrayView<const double>(eValues, numEntries);
    }

    //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
    /*!
      \param LocalRow - (In) Local row number for which indices are desired.
      \param Indices  - (Out) Global column indices corresponding to values.
      \param Values   - (Out) Row values
      \pre <tt>isGloballyIndexed() == false</tt>
      \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

      Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
    */
    inline void getLocalRowView(int LocalRow, ArrayView<const int> &indices, ArrayView<const double> &values) const { 
      CTHULHU_DEBUG_ME; 

      int      numEntries;
      double * eValues;
      int    * eIndices;
      
      CTHULHU_ERR_CHECK(mtx_->ExtractMyRowView(LocalRow, numEntries, eValues, eIndices));
      if (numEntries == 0) { eValues = NULL; eIndices = NULL; } // Cf. TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

      indices = ArrayView<const int>(eIndices, numEntries);
      values  = ArrayView<const double>(eValues, numEntries);
    }

    //! \brief Get a copy of the diagonal entries owned by this node, with local row indices.
    /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the 
      the zero and non-zero diagonals owned by this node. */
    inline void getLocalDiagCopy(Vector<double,int,int> &diag) const { 
      CTHULHU_DEBUG_ME; 
      
      CTHULHU_DYNAMIC_CAST(EpetraVector, diag, eDiag, "Cthulhu::EpetraCrsMatrix.getLocalDiagCopy() only accept Cthulhu::EpetraVector as input arguments.");
      mtx_->ExtractDiagonalCopy(*eDiag.getEpetra_Vector()); 
    }

    //@}

    //! @name Advanced Matrix-vector multiplication and solve methods
    //@{

    //! Multiplies this matrix by a MultiVector.
    /*! \c X is required to be post-imported, i.e., described by the column map of the matrix. \c Y is required to be pre-exported, i.e., described by the row map of the matrix.

    Both are required to have constant stride, and they are not permitted to ocupy overlapping space. No runtime checking will be performed in a non-debug build.

    This method is templated on the scalar type of MultiVector objects, allowing this method to be applied to MultiVector objects of arbitrary type. However, it is recommended that multiply() not be called directly; instead, use the EpetraCrsMatrixMultiplyOp, as it will handle the import/exprt operations required to apply a matrix with non-trivial communication needs.
          
    If \c beta is equal to zero, the operation will enjoy overwrite semantics (\c Y will be overwritten with the result of the multiplication). Otherwise, the result of the multiplication
    will be accumulated into \c Y.
    */
#ifdef CTHULHU_NOT_IMPLEMENTED
    template <class Domaindouble, class Rangedouble>
    inline void multiply(const MultiVector<Domaindouble,int,int> & X, MultiVector<Rangedouble,int,int> &Y, Teuchos::ETransp trans, Rangedouble alpha, Rangedouble beta) const { CTHULHU_DEBUG_ME; mtx_->multiply(X, Y, trans, alpha, beta); }
#endif // CTHULHU_NOT_IMPLEMENTED

    // TODO Note: Do we need to use a Tpetra::CrsMatrixMultiplyOp ?? 
    //            (Epetra Doc of multiply: it is recommended that multiply() not be called directly; instead, use the CrsMatrixMultiplyOp, as it will handle the import/exprt operations required to apply a matrix with non-trivial communication needs.)

    // TODO : templated type
    inline void multiply(const MultiVector<double,int,int> & X, MultiVector<double,int,int> &Y, Teuchos::ETransp trans, double alpha, double beta) const {
      CTHULHU_DEBUG_ME; 

      TEST_FOR_EXCEPTION((alpha != 1) || (beta != 0), Cthulhu::Exceptions::NotImplemented, "Cthulhu::EpetraCrsMatrix.multiply() only accept alpha==1 and beta==0");
      
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, X, eX, "Cthulhu::EpetraCrsMatrix->multiply() only accept Cthulhu::EpetraMultiVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(      EpetraMultiVector, Y, eY, "Cthulhu::EpetraCrsMatrix->multiply() only accept Cthulhu::EpetraMultiVector as input arguments.");

      TEST_FOR_EXCEPTION((trans != Teuchos::NO_TRANS) && (trans == Teuchos::TRANS), Cthulhu::Exceptions::NotImplemented, "Cthulhu::EpetraCrsMatrix->multiply() only accept trans == NO_TRANS or trans == TRANS");
      bool eTrans = Teuchos2Epetra_Trans(trans);

      CTHULHU_ERR_CHECK(mtx_->Multiply(eTrans, *eX.getEpetra_MultiVector(), *eY.getEpetra_MultiVector()));
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Solves a linear system when the underlying matrix is triangular.
    /*! \c X is required to be post-imported, i.e., described by the column map of the matrix. \c Y is required to be pre-exported, i.e., described by the row map of the matrix.

    This method is templated on the scalar type of MultiVector objects, allowing this method to be applied to MultiVector objects of arbitrary type. However, it is recommended that solve() not be called directly; instead, use the EpetraCrsMatrixSolveOp, as it will handle the import/exprt operations required to apply a matrix with non-trivial communication needs.
          
    Both are required to have constant stride. However, unlike multiply(), it is permissible for <tt>&X == &Y</tt>. No runtime checking will be performed in a non-debug build.
    */
    template <class Domaindouble, class Rangedouble>
    inline void solve(const MultiVector<Rangedouble,int,int> & Y, MultiVector<Domaindouble,int,int> &X, Teuchos::ETransp trans) const { CTHULHU_DEBUG_ME; mtx_->solve(Y, X, trans); }
#endif // CTHULHU_NOT_IMPLEMENTED
          
    //@}

    //! @name Methods implementing Operator
    //@{ 

    //! \brief Computes the sparse matrix-multivector multiplication.
    /*! Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exceptions:
      - if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
    */
    inline void apply(const MultiVector<double,int,int> & X, MultiVector<double,int,int> &Y, 
                      Teuchos::ETransp mode = Teuchos::NO_TRANS,
                      double alpha = ScalarTraits<double>::one(),
                      double beta = ScalarTraits<double>::zero()) const { 
      CTHULHU_DEBUG_ME; 

      TEST_FOR_EXCEPTION((alpha != 1) || (beta != 0), Cthulhu::Exceptions::NotImplemented, "Cthulhu::EpetraCrsMatrix.multiply() only accept alpha==1 and beta==0");
      
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, X, eX, "Cthulhu::EpetraCrsMatrix->apply() only accept Cthulhu::EpetraMultiVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(      EpetraMultiVector, Y, eY, "Cthulhu::EpetraCrsMatrix->apply() only accept Cthulhu::EpetraMultiVector as input arguments.");

      TEST_FOR_EXCEPTION((mode != Teuchos::NO_TRANS) && (mode != Teuchos::TRANS), Cthulhu::Exceptions::NotImplemented, "Cthulhu::EpetraCrsMatrix->apply() only accept mode == NO_TRANS or mode == TRANS");
      bool eTrans = Teuchos2Epetra_Trans(mode);

      // /!\ UseTranspose value
      TEST_FOR_EXCEPTION(mtx_->UseTranspose(), Cthulhu::Exceptions::NotImplemented, "An exception is throw to let you know that Cthulhu::EpetraCrsMatrix->apply() do not take into account the UseTranspose() parameter of Epetra_CrsMatrix.");
      
      CTHULHU_ERR_CHECK(mtx_->Multiply(eTrans, *eX.getEpetra_MultiVector(), *eY.getEpetra_MultiVector()));
    }

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Indicates whether this operator supports applying the adjoint operator.
    inline bool hasTransposeApply() const { CTHULHU_DEBUG_ME; return mtx_->hasTransposeApply(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //! \brief Returns the Map associated with the domain of this operator.
    //! This will be <tt>null</tt> until fillComplete() is called.
    inline const RCP<const Map<int,int> > getDomainMap() const { 
      CTHULHU_DEBUG_ME; 

      // The Epetra map is owned by the underlying Epetra matrix and freed when the matrix is deleted.
      // I have to make a copy of the map to be sure that the RCP<Map> returned by this method will remain valid even if the matrix is deleted.
      // Note that it is not a problem because the copy constructor of Epetra_Map doesn't really copy the data (there is an reference count mecanism in Epetra).

      RCP<const Epetra_Map> map = rcp(new Epetra_Map(mtx_->DomainMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }
    
    //! Returns the Map associated with the domain of this operator.
    //! This will be <tt>null</tt> until fillComplete() is called.
    inline const RCP<const Map<int,int> > getRangeMap() const { 
      CTHULHU_DEBUG_ME; 

      // The Epetra map is owned by the underlying Epetra matrix and freed when the matrix is deleted.
      // I have to make a copy of the map to be sure that the RCP<Map> returned by this method will remain valid even if the matrix is deleted.
      // Note that it is not a problem because the copy constructor of Epetra_Map doesn't really copy the data (there is an reference count mecanism in Epetra).

      RCP<const Epetra_Map> map = rcp(new Epetra_Map(mtx_->RangeMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }

    //@}

    //! @name Overridden from Teuchos::Describable 
    //@{

    /** \brief Return a simple one-line description of this object. */
    inline std::string description() const { 
      CTHULHU_DEBUG_ME; 

      // This implementation come from Tpetra_CrsMatrix_def.hpp (without modification)
      std::ostringstream oss;
      //TODO: oss << DistObject<char, LocalOrdinal,GlobalOrdinal>::description();
      if (isFillComplete()) {
        oss << "{status = fill complete"
            << ", global rows = " << getGlobalNumRows()
            << ", global cols = " << getGlobalNumCols()
            << ", global num entries = " << getGlobalNumEntries()
            << "}";
      }
      else {
        oss << "{status = fill not complete"
            << ", global rows = " << getGlobalNumRows()
            << "}";
      }
      return oss.str();
      
    } 
    
    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    inline void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const { 
      CTHULHU_DEBUG_ME; 

      typedef int LocalOrdinal;
      typedef int GlobalOrdinal;
      typedef double Scalar;

      // This implementation come from Tpetra_CrsMatrix_def.hpp (without modification)
      using std::endl;
      using std::setw;
      using Teuchos::VERB_DEFAULT;
      using Teuchos::VERB_NONE;
      using Teuchos::VERB_LOW;
      using Teuchos::VERB_MEDIUM;
      using Teuchos::VERB_HIGH;
      using Teuchos::VERB_EXTREME;
      Teuchos::EVerbosityLevel vl = verbLevel;
      if (vl == VERB_DEFAULT) vl = VERB_LOW;
      RCP<const Comm<int> > comm = this->getComm();
      const int myImageID = comm->getRank(),
        numImages = comm->getSize();
      size_t width = 1;
      for (size_t dec=10; dec<getGlobalNumRows(); dec *= 10) {
        ++width;
      }
      width = std::max<size_t>(width,11) + 2;
      Teuchos::OSTab tab(out);
      //    none: print nothing
      //     low: print O(1) info from node 0
      //  medium: print O(P) info, num entries per node
      //    high: print O(N) info, num entries per row
      // extreme: print O(NNZ) info: print indices and values
      // 
      // for medium and higher, print constituent objects at specified verbLevel
      if (vl != VERB_NONE) {
        if (myImageID == 0) out << this->description() << std::endl; 
        // O(1) globals, minus what was already printed by description()
        if (isFillComplete() && myImageID == 0) {
          out << "Global number of diagonals = " << getGlobalNumDiags() << std::endl;
          out << "Global max number of entries = " << getGlobalMaxNumRowEntries() << std::endl;
        }
        // constituent objects
        if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
          if (myImageID == 0) out << "\nRow map: " << std::endl;
          getRowMap()->describe(out,vl);
          //
          if (getColMap() != null) {
            if (getColMap() == getRowMap()) {
              if (myImageID == 0) out << "\nColumn map is row map.";
            }
            else {
              if (myImageID == 0) out << "\nColumn map: " << std::endl;
              getColMap()->describe(out,vl);
            }
          }
          if (getDomainMap() != null) {
            if (getDomainMap() == getRowMap()) {
              if (myImageID == 0) out << "\nDomain map is row map.";
            }
            else if (getDomainMap() == getColMap()) {
              if (myImageID == 0) out << "\nDomain map is row map.";
            }
            else {
              if (myImageID == 0) out << "\nDomain map: " << std::endl;
              getDomainMap()->describe(out,vl);
            }
          }
          if (getRangeMap() != null) {
            if (getRangeMap() == getDomainMap()) {
              if (myImageID == 0) out << "\nRange map is domain map." << std::endl;
            }
            else if (getRangeMap() == getRowMap()) {
              if (myImageID == 0) out << "\nRange map is row map." << std::endl;
            }
            else {
              if (myImageID == 0) out << "\nRange map: " << std::endl;
              getRangeMap()->describe(out,vl);
            }
          }
          if (myImageID == 0) out << std::endl;
        }
        // O(P) data
        if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
          for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
            if (myImageID == imageCtr) {
              out << "Node ID = " << imageCtr << std::endl;
// TODO: need a graph
//               if (staticGraph_->indicesAreAllocated() == false) {
//                 out << "Node not allocated" << std::endl;
//               }
//               else {
//                 out << "Node number of allocated entries = " << staticGraph_->getNodeAllocationSize() << std::endl;
//               }

// TMP:
//            const Epetra_CrsGraph & staticGraph_ = mtx_->Graph();
// End of TMP

              out << "Node number of entries = " << getNodeNumEntries() << std::endl;
              if (isFillComplete()) {
                out << "Node number of diagonals = " << getNodeNumDiags() << std::endl;
              }
              out << "Node max number of entries = " << getNodeMaxNumRowEntries() << std::endl;
            }
            comm->barrier();
            comm->barrier();
            comm->barrier();
          }
        }
        // O(N) and O(NNZ) data
        if (vl == VERB_HIGH || vl == VERB_EXTREME) {
          for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
            if (myImageID == imageCtr) {
              out << std::setw(width) << "Node ID"
                  << std::setw(width) << "Global Row" 
                  << std::setw(width) << "Num Entries";
              if (vl == VERB_EXTREME) {
                out << std::setw(width) << "(Index,Value)";
              }
              out << std::endl;
              for (size_t r=0; r < getNodeNumRows(); ++r) {
                const size_t nE = getNumEntriesInLocalRow(r);
                GlobalOrdinal gid = getRowMap()->getGlobalElement(r);
                out << std::setw(width) << myImageID 
                    << std::setw(width) << gid
                    << std::setw(width) << nE;
                if (vl == VERB_EXTREME) {
                  if (isGloballyIndexed()) {
                    ArrayView<const GlobalOrdinal> rowinds;
                    ArrayView<const Scalar> rowvals;
                    getGlobalRowView(gid,rowinds,rowvals);
                    for (size_t j=0; j < nE; ++j) {
                      out << " (" << rowinds[j]
                          << ", " << rowvals[j]
                          << ") ";
                    }
                  }
                  else if (isLocallyIndexed()) {
                    ArrayView<const LocalOrdinal> rowinds;
                    ArrayView<const Scalar> rowvals;
                    getLocalRowView(r,rowinds,rowvals);
                    for (size_t j=0; j < nE; ++j) {
                      out << " (" << getColMap()->getGlobalElement(rowinds[j]) 
                          << ", " << rowvals[j]
                          << ") ";
                    }
                  }
                }
                out << std::endl;
              }
            }
            comm->barrier();
            comm->barrier();
            comm->barrier();
          }
        }
      }
    
    }
    
    //@}

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! @name Methods implementing Cthulhu::DistObject
    //@{
    inline bool checkSizes(const DistObject<char, int,int>& source) { CTHULHU_DEBUG_ME; return mtx_->checkSizes(); }

    inline void copyAndPermute(const DistObject<char, int,int>& source,
                               size_t numSameIDs,
                               const ArrayView<const int> &permuteToLIDs,
                               const ArrayView<const int> &permuteFromLIDs) { CTHULHU_DEBUG_ME; mtx_->copyAndPermute(); }

    inline void packAndPrepare(const DistObject<char, int,int>& source,
                               const ArrayView<const int> &exportLIDs,
                               Array<char> &exports,
                               const ArrayView<size_t> & numPacketsPerLID,
                               size_t& constantNumPackets,
                               Distributor &distor) { CTHULHU_DEBUG_ME; mtx_->packAndPrepare(); }

    inline void unpackAndCombine(const ArrayView<const int> &importLIDs,
                                 const ArrayView<const char> &imports,
                                 const ArrayView<size_t> &numPacketsPerLID,
                                 size_t constantNumPackets,
                                 Distributor &distor,
                                 CombineMode CM) { CTHULHU_DEBUG_ME; mtx_->unpackAndCombine(); }
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
    inline CTHULHU_DEPRECATED void optimizeStorage() { CTHULHU_DEBUG_ME; }
    
    //! Deprecated. Get a persisting const view of the entries in a specified global row of this matrix.
    inline CTHULHU_DEPRECATED void getGlobalRowView(int GlobalRow, ArrayRCP<const int> &indices, ArrayRCP<const double> &values) const { CTHULHU_DEBUG_ME; } 
    
    //! Deprecated. Get a persisting const view of the entries in a specified local row of this matrix.
    inline CTHULHU_DEPRECATED void getLocalRowView(int LocalRow, ArrayRCP<const int> &indices, ArrayRCP<const double> &values) const { CTHULHU_DEBUG_ME; } 
    
    //@}
#endif // CTHULHU_NOT_IMPLEMENTED

    RCP<const Epetra_CrsMatrix> getEpetra_CrsMatrix() const { CTHULHU_DEBUG_ME; return mtx_; }

    RCP<Epetra_CrsMatrix> getEpetra_CrsMatrixNonConst() const { CTHULHU_DEBUG_ME; return mtx_; }

    /** TODO : interface of Teuchos_LabeledObject.hpp **/
    void setObjectLabel (const std::string &objectLabel) { CTHULHU_DEBUG_ME; //mtx_->setObjectLabel(objectLabel); TODO
    }

  private:
    
    RCP<Epetra_CrsMatrix> mtx_;

  }; // class EpetraCrsMatrix

} // namespace Cthulhu

#endif
