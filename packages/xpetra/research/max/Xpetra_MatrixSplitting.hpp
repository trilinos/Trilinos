#ifndef XPETRA_MATRIXSPLITTING_HPP
#define XPETRA_MATRIXSPLITTING_HPP

#include "Xpetra_Map.hpp"
#include "Xpetra_Matrix.hpp"

/** \file Xpetra_MatrixSplitting.hpp

Declarations for the class Xpetra::MatrixSplitting.
*/
namespace Xpetra {

    /*!
     @class Xpetra::MatrixSplitting class.
     @brief Xpetra-specific matrix class.

     This class is specific to Xpetra and has no analogue in Epetra or Tpetra.  The main motivation for this class is to be able to access matrix data in a manner different than how it is stored.
     For example, it might be more convenient to treat ("view") a matrix stored in compressed row storage as if it were a block matrix.  The Xpetra::MatrixSplitting class is intended to manage these "views".

     <B>How to create a Matrix from an existing CrsMatrix</B>

    */

  typedef std::string viewLabel_t;

  template <class Scalar        = Operator<>::scalar_type,
            class LocalOrdinal  = Operator<>::local_ordinal_type,
            class GlobalOrdinal = typename Operator<LocalOrdinal>::global_ordinal_type,
            class Node          = typename Operator<LocalOrdinal, GlobalOrdinal>::node_type>
  class MatrixSplitting : public Matrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > {

    typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
		typedef Xpetra::Matrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > Matrix;
    typedef Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> CrsGraph;
    typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrix;
    typedef Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrixFactory;
    typedef Xpetra::MatrixView<Scalar, LocalOrdinal, GlobalOrdinal, Node> MatrixView;

    //! @name Constructor/Destructor Methods
    //@{

		//! Constructor specifying fixed number of entries for each row.
		MatrixSplitting (const RCP<const Map>& rowMap,
		               size_t maxNumEntriesPerRow,
		               Xpetra::ProfileType pftype = Xpetra::DynamicProfile)
		  : finalDefaultView_ (false)
		{
		  matrixData_ = CrsMatrixFactory::Build (rowMap, maxNumEntriesPerRow, pftype);
		  CreateDefaultView ();
		}

		//! Constructor specifying (possibly different) number of entries in each row.
		MatrixSplitting (const RCP<const Map>& rowMap,
		               const ArrayRCP<const size_t>& NumEntriesPerRowToAlloc,
		               ProfileType pftype = Xpetra::DynamicProfile)
		  : finalDefaultView_ (false)
		{
		  matrixData_ = CrsMatrixFactory::Build(rowMap, NumEntriesPerRowToAlloc, pftype);
		  CreateDefaultView ();
		}

		//! Constructor specifying fixed number of entries for each row and column map
		MatrixSplitting(const RCP<const Map> &rowMap, const RCP<const Map>& colMap, size_t maxNumEntriesPerRow, Xpetra::ProfileType pftype = Xpetra::DynamicProfile)
		  : finalDefaultView_(false)
		{
		  // Set matrix data
		  matrixData_ = CrsMatrixFactory::Build(rowMap, colMap, maxNumEntriesPerRow, pftype);

		  // Default view
		  CreateDefaultView();
		}

		//! Constructor specifying fixed number of entries for each row and column map
		MatrixSplitting(const RCP<const Map> &rowMap, const RCP<const Map>& colMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, Xpetra::ProfileType pftype = Xpetra::DynamicProfile)
		  : finalDefaultView_(false)
		{
		  // Set matrix data
		  matrixData_ = CrsMatrixFactory::Build(rowMap, colMap, NumEntriesPerRowToAlloc, pftype);

		  // Default view
		  CreateDefaultView();
		}

	#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
	#ifdef HAVE_XPETRA_TPETRA
		//! Constructor specifying fixed number of entries for each row and column map
		MatrixSplitting(const RCP<const Map> &rowMap, const RCP<const Map>& colMap, const local_matrix_type& lclMatrix, const Teuchos::RCP<Teuchos::ParameterList>& params = null)
		  : finalDefaultView_(false)
		{
		  // Set matrix data
		  matrixData_ = CrsMatrixFactory::Build(rowMap, colMap, lclMatrix, params);

		  // Default view
		  CreateDefaultView();
		}
	#else
	#ifdef __GNUC__
	#warning "Xpetra Kokkos interface for CrsMatrix is enabled (HAVE_XPETRA_KOKKOS_REFACTOR) but Tpetra is disabled. The Kokkos interface needs Tpetra to be enabled, too."
	#endif
	#endif
	#endif

		MatrixSplitting(RCP<CrsMatrix> matrix)
		  : finalDefaultView_(matrix->isFillComplete())
		{
		  // Set matrix data
		  matrixData_ = matrix;

		  // Default view
		  CreateDefaultView();
		}

		MatrixSplitting(const RCP<const CrsGraph>& graph, const RCP<ParameterList>& paramList = Teuchos::null)
		  : finalDefaultView_(false)
		{
		  // Set matrix data
		  matrixData_ = CrsMatrixFactory::Build(graph, paramList);

		  // Default view
		  CreateDefaultView();
		}


    //! Destructor
    virtual ~MatrixSplitting() { }

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
		void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &cols, const ArrayView<const Scalar> &vals) {
		  matrixData_->insertGlobalValues(globalRow, cols, vals);
		}

		//! Insert matrix entries, using local IDs.
		/** All index values must be in the local space.
		    \pre \c localRow exists as an ID in the global row map
		    \pre <tt>isGloballyIndexed() == false</tt>
		    \pre <tt>isStorageOptimized() == false</tt>

		    \post <tt>isLocallyIndexed() == true</tt>
		*/
		void insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &cols, const ArrayView<const Scalar> &vals) {
		  matrixData_->insertLocalValues(localRow, cols, vals);
		}

		//! \brief Replace matrix entries, using global IDs.
		/** All index values must be in the global space.

		\pre \c globalRow is a global row belonging to the matrix on this node.

		\note If (globalRow,cols[i]) corresponds to an entry that is duplicated in this matrix row (likely because it was inserted more than once and fillComplete() has not been called in the interim), the behavior of this function is not defined. */
		void replaceGlobalValues(GlobalOrdinal globalRow,
		                         const ArrayView<const GlobalOrdinal> &cols,
		                         const ArrayView<const Scalar>        &vals) { matrixData_->replaceGlobalValues(globalRow, cols, vals); }

		//! Replace matrix entries, using local IDs.
		/** All index values must be in the local space.
		    Note that if a value is not already present for the specified location in the matrix, the input value will be ignored silently.
		*/
		void replaceLocalValues(LocalOrdinal localRow,
		                        const ArrayView<const LocalOrdinal> &cols,
		                        const ArrayView<const Scalar>       &vals) { matrixData_->replaceLocalValues(localRow, cols, vals); }

		//! Set all matrix entries equal to scalar
		virtual void setAllToScalar(const Scalar &alpha) { matrixData_->setAllToScalar(alpha); }

		//! Scale the current values of a matrix, this = alpha*this.
		void scale(const Scalar &alpha) {
		  matrixData_->scale(alpha);
		}

		//@}

		//! @name Transformational Methods
		//@{

		/*! Resume fill operations.
		  After calling fillComplete(), resumeFill() must be called before initiating any changes to the matrix.

		  resumeFill() may be called repeatedly.

		  \post  <tt>isFillActive() == true<tt>
		  \post  <tt>isFillComplete() == false<tt>
		*/
		void resumeFill(const RCP< ParameterList > &params=null) {
		  matrixData_->resumeFill(params);
		}

		/*! \brief Signal that data entry is complete, specifying domain and range maps.

		Off-node indices are distributed (via globalAssemble()), indices are sorted, redundant indices are eliminated, and global indices are transformed to local indices.

		\pre  <tt>isFillActive() == true<tt>
		\pre <tt>isFillComplete()() == false<tt>

		\post <tt>isFillActive() == false<tt>
		\post <tt>isFillComplete() == true<tt>
		\post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
		*/
		void fillComplete(const RCP<const Map> &domainMap, const RCP<const Map> &rangeMap, const RCP<Teuchos::ParameterList> &params = null) {
		  matrixData_->fillComplete(domainMap, rangeMap, params);

		  // Update default view with the colMap because colMap can be <tt>null</tt> until fillComplete() is called.
		  updateDefaultView();
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
		//TODO : Get ride of "Tpetra"::OptimizeOption
		void fillComplete(const RCP<ParameterList> &params = null) {
		  matrixData_->fillComplete(params);

		  // Update default view with the colMap because colMap can be <tt>null</tt> until fillComplete() is called.
		  updateDefaultView();
		}

		//@}

		//! Returns the number of global rows in this matrix.
		/** Undefined if isFillActive().
		 */
		global_size_t getGlobalNumRows() const {
		  return matrixData_->getGlobalNumRows();
		}

		//! \brief Returns the number of global columns in the matrix.
		/** Undefined if isFillActive().
		 */
		global_size_t getGlobalNumCols() const {
		  return matrixData_->getGlobalNumCols();
		}

		//! Returns the number of matrix rows owned on the calling node.
		size_t getNodeNumRows() const {
		  return matrixData_->getNodeNumRows();
		}

		//! Returns the global number of entries in this matrix.
		global_size_t getGlobalNumEntries() const {
		  return matrixData_->getGlobalNumEntries();
		}

		//! Returns the local number of entries in this matrix.
		size_t getNodeNumEntries() const {
		  return matrixData_->getNodeNumEntries();
		}

		//! Returns the current number of entries on this node in the specified local row.
		/*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
		size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const {
		  return matrixData_->getNumEntriesInLocalRow(localRow);
		}

		//! \brief Returns the number of global diagonal entries, based on global row/column index comparisons.
		/** Undefined if isFillActive().
		 */
		global_size_t getGlobalNumDiags() const {
		  return matrixData_->getGlobalNumDiags();
		}

		//! \brief Returns the number of local diagonal entries, based on global row/column index comparisons.
		/** Undefined if isFillActive().
		 */
		size_t getNodeNumDiags() const {
		  return matrixData_->getNodeNumDiags();
		}

		//! \brief Returns the maximum number of entries across all rows/columns on all nodes.
		/** Undefined if isFillActive().
		 */
		size_t getGlobalMaxNumRowEntries() const {
		  return matrixData_->getGlobalMaxNumRowEntries();
		}

		//! \brief Returns the maximum number of entries across all rows/columns on this node.
		/** Undefined if isFillActive().
		 */
		size_t getNodeMaxNumRowEntries() const {
		  return matrixData_->getNodeMaxNumRowEntries();
		}

		//! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
		bool isLocallyIndexed() const {
		  return matrixData_->isLocallyIndexed();
		}

		//! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
		bool isGloballyIndexed() const {
		  return matrixData_->isGloballyIndexed();
		}

		//! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
		bool isFillComplete() const {
		  return matrixData_->isFillComplete();
		}

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
		void getLocalRowCopy(LocalOrdinal LocalRow,
		                     const ArrayView<LocalOrdinal> &Indices,
		                     const ArrayView<Scalar> &Values,
		                     size_t &NumEntries
		                     ) const {
		  matrixData_->getLocalRowCopy(LocalRow, Indices, Values, NumEntries);
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
		void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &indices, ArrayView<const Scalar> &values) const {
		   matrixData_->getGlobalRowView(GlobalRow, indices, values);
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
		void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices, ArrayView<const Scalar> &values) const {
		   matrixData_->getLocalRowView(LocalRow, indices, values);
		}

		//! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
		/*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the
		  the zero and non-zero diagonals owned by this node. */
		void getLocalDiagCopy(Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const {
		  matrixData_->getLocalDiagCopy(diag);
		}

		//! Get offsets of the diagonal entries in the matrix.
		void getLocalDiagOffsets(Teuchos::ArrayRCP<size_t> &offsets) const {
		  matrixData_->getLocalDiagOffsets(offsets);
		}

		//! Get a copy of the diagonal entries owned by this node, with local row indices, using row offsets.
		void getLocalDiagCopy(Xpetra::Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &diag, const Teuchos::ArrayView<const size_t> &offsets) const {
		  matrixData_->getLocalDiagCopy(diag,offsets);
		}

		//! Get Frobenius norm of the matrix
		typename ScalarTraits<Scalar>::magnitudeType getFrobeniusNorm() const {
		  return matrixData_->getFrobeniusNorm();
		}

		//! Left scale matrix using the given vector entries
		void leftScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
		  matrixData_->leftScale(x);
		}

		//! Right scale matrix using the given vector entries
		void rightScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
		  matrixData_->rightScale(x);
		}

		//! Returns true if globalConstants have been computed; false otherwise
		bool haveGlobalConstants() const {
		  return matrixData_->haveGlobalConstants();
		}

		//@}

	private:

		// Default view is created after fillComplete()
		// Because ColMap might not be available before fillComplete().
		void CreateDefaultView() {

		  // Create default view
		  this->defaultViewLabel_ = "point";
		  this->CreateView(this->GetDefaultViewLabel(), matrixData_->getRowMap(), matrixData_->getColMap());

		  // Set current view
		  this->currentViewLabel_ = this->GetDefaultViewLabel();
		}

		// The colMap can be <tt>null</tt> until fillComplete() is called. The default view of the Matrix have to be updated when fillComplete() is called.
		// If CrsMatrix::fillComplete() have been used instead of CrsMatrixWrap::fillComplete(), the default view is updated when getColMap() is called.
		void updateDefaultView() const {
		  if ((finalDefaultView_ == false) &&  matrixData_->isFillComplete() ) {
		    // Update default view with the colMap
		    Matrix::operatorViewTable_.get(Matrix::GetDefaultViewLabel())->SetColMap(matrixData_->getColMap());
		    finalDefaultView_ = true;
		  }
		}

  // The boolean finalDefaultView_ keep track of the status of the default view (= already updated or not)
  // See also MatrixSplitting::updateDefaultView()
  mutable bool finalDefaultView_;

	RCP<CrsMatrix> matrixData_;

  }; //class MatrixSplitting

} //namespace Xpetra

#define XPETRA_MATRIX_SHORT
#endif //XPETRA_MATRIX_DECL_HPP
