#ifndef XPETRA_MATRIXSPLITTING_HPP
#define XPETRA_MATRIXSPLITTING_HPP

//Xpetra
#include "Xpetra_Map.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_Matrix.hpp"
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include "Xpetra_IO.hpp"
#include "Xpetra_SplittingDriver_def.hpp"

//Epetra
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"


#include "Ifpack2_OverlappingPartitioner_def.hpp"

/** \file Xpetra_MatrixSplitting.hpp

Declarations for the class Xpetra::MatrixSplitting.
*/
namespace Xpetra {


//Definition of the predicate for the node_ structure.
//Given a tuple made of node index and a specific region it belongs to,
//this predicate returns true if the node has global index which coincides with the index specified in input.
template<class GlobalOrdinal>
class checkerRegionToAll { 
 
	public:  

		//Constructor
		checkerRegionToAll( GlobalOrdinal node_index){node_index_ = node_index;};

		//Unary Operator
  		bool operator()(const std::tuple<GlobalOrdinal, GlobalOrdinal > &node)  
  		{ return (std::get<0>(node) == node_index_); }  

	private:

		GlobalOrdinal node_index_;

}; 


//Definition of the predicate for the node_ structure.
//Given a tuple made of node index and a specific region it belongs to,
//this predicate returns true if the node has global index which coincides with the index specified in input.
//It does the same thing as checkerRegionToAll but it works on a different data structure 
template<class GlobalOrdinal>
class checkerAllToRegion { 
 
	public:  

		//Constructor
		checkerAllToRegion( GlobalOrdinal node_index){node_index_ = node_index;};

		//Unary Operator
  		bool operator()(const std::tuple<GlobalOrdinal, GlobalOrdinal > &node)  
  		{ return (std::get<1>(node) == node_index_); }  

	private:

		GlobalOrdinal node_index_;

}; 

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
	    class Node          = typename Operator<LocalOrdinal, GlobalOrdinal>::node_type,
	    UnderlyingLib lib		= Xpetra::UseEpetra>
	class MatrixSplitting : public Matrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > {

	typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
	typedef Xpetra::Matrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > Matrix;
	typedef Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> CrsGraph;
	typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrix;
	typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrixWrap;
	typedef Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrixFactory;
	typedef Xpetra::MatrixView<Scalar, LocalOrdinal, GlobalOrdinal, Node> MatrixView;

	public:

		//! @name Constructor/Destructor Methods
		//@{

		//! Constructor specifying fixed number of entries for each row.
		MatrixSplitting(RCP<Matrix> matrix, RCP< Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > > nodes)
		{
			std::cout<<"This version of MatrixSplitting constructor is NOT currently supported \n";			
		}
		//
		//
		MatrixSplitting(const char* matrix_file_name, const char* elements_file_name, RCP<const Teuchos::Comm<int> > comm)
		{
			comm_ = comm;
			driver_ = rcp( new Xpetra::SplittingDriver<Scalar, LocalOrdinal, GlobalOrdinal, Node> (elements_file_name, comm_) );
			Array<GlobalOrdinal> elementlist = driver_->GetGlobalRowMap();
			num_total_elements_ = driver_->GetNumGlobalElements();
			num_total_regions_ = driver_->GetNumTotalRegions();

			//Create Xpetra map for global stiffness matrix
			RCP<const Xpetra::Map<int,GlobalOrdinal,Node> > xpetraMap;
			xpetraMap = Xpetra::MapFactory<int,GlobalOrdinal,Node>::Build(lib, num_total_elements_, elementlist, 0, comm); 

			//Import matrix from an .mtx file into an Xpetra wrapper for an Epetra matrix
			globalMatrixData_ = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Read(matrix_file_name, xpetraMap);

			CreateLocalMatrices( driver_->GetLocalRowMaps() );
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
		  globalMatrixData_->insertGlobalValues(globalRow, cols, vals);
		}

		//! Insert matrix entries, using local IDs.
		/** All index values must be in the local space.
		    \pre \c localRow exists as an ID in the global row map
		    \pre <tt>isGloballyIndexed() == false</tt>
		    \pre <tt>isStorageOptimized() == false</tt>

		    \post <tt>isLocallyIndexed() == true</tt>
		*/
		void insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &cols, const ArrayView<const Scalar> &vals) {
		  globalMatrixData_->insertLocalValues(localRow, cols, vals);
		}

		//! \brief Replace matrix entries, using global IDs.
		/** All index values must be in the global space.

		\pre \c globalRow is a global row belonging to the matrix on this node.

		\note If (globalRow,cols[i]) corresponds to an entry that is duplicated in this matrix row (likely because it was inserted more than once and fillComplete() has not been called in the interim), the behavior of this function is not defined. */
		void replaceGlobalValues(GlobalOrdinal globalRow,
		                         const ArrayView<const GlobalOrdinal> &cols,
		                         const ArrayView<const Scalar>        &vals) { globalMatrixData_->replaceGlobalValues(globalRow, cols, vals); }

		//! Replace matrix entries, using local IDs.
		/** All index values must be in the local space.
		    Note that if a value is not already present for the specified location in the matrix, the input value will be ignored silently.
		*/
		void replaceLocalValues(LocalOrdinal localRow,
		                        const ArrayView<const LocalOrdinal> &cols,
		                        const ArrayView<const Scalar>       &vals) { globalMatrixData_->replaceLocalValues(localRow, cols, vals); }

		//! Set all matrix entries equal to scalar
		virtual void setAllToScalar(const Scalar &alpha) { globalMatrixData_->setAllToScalar(alpha); }

		//! Scale the current values of a matrix, this = alpha*this.
		void scale(const Scalar &alpha) {
		  globalMatrixData_->scale(alpha);
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
		  globalMatrixData_->resumeFill(params);
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
		  globalMatrixData_->fillComplete(domainMap, rangeMap, params);

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
		  globalMatrixData_->fillComplete(params);

		  // Update default view with the colMap because colMap can be <tt>null</tt> until fillComplete() is called.
		  updateDefaultView();
		}

		//@}

		//! Returns the number of global rows in this matrix.
		/** Undefined if isFillActive().
		 */
		global_size_t getGlobalNumRows() const {
		  return globalMatrixData_->getGlobalNumRows();
		}

		//! \brief Returns the number of global columns in the matrix.
		/** Undefined if isFillActive().
		 */
		global_size_t getGlobalNumCols() const {
		  return globalMatrixData_->getGlobalNumCols();
		}

		//! Returns the number of matrix rows owned on the calling node.
		size_t getNodeNumRows() const {
		  return globalMatrixData_->getNodeNumRows();
		}

		//! Returns the global number of entries in this matrix.
		global_size_t getGlobalNumEntries() const {
		  return globalMatrixData_->getGlobalNumEntries();
		}

		//! Returns the local number of entries in this matrix.
		size_t getNodeNumEntries() const {
		  return globalMatrixData_->getNodeNumEntries();
		}

		//! Returns the current number of entries on this node in the specified local row.
		/*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
		size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const {
		  return globalMatrixData_->getNumEntriesInLocalRow(localRow);
		}

		//! \brief Returns the number of global diagonal entries, based on global row/column index comparisons.
		/** Undefined if isFillActive().
		 */
		global_size_t getGlobalNumDiags() const {
		  return globalMatrixData_->getGlobalNumDiags();
		}

		//! \brief Returns the number of local diagonal entries, based on global row/column index comparisons.
		/** Undefined if isFillActive().
		 */
		size_t getNodeNumDiags() const {
		  return globalMatrixData_->getNodeNumDiags();
		}

		//! \brief Returns the maximum number of entries across all rows/columns on all nodes.
		/** Undefined if isFillActive().
		 */
		size_t getGlobalMaxNumRowEntries() const {
		  return globalMatrixData_->getGlobalMaxNumRowEntries();
		}

		//! \brief Returns the maximum number of entries across all rows/columns on this node.
		/** Undefined if isFillActive().
		 */
		size_t getNodeMaxNumRowEntries() const {
		  return globalMatrixData_->getNodeMaxNumRowEntries();
		}

		//! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
		bool isLocallyIndexed() const {
		  return globalMatrixData_->isLocallyIndexed();
		}

		//! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
		bool isGloballyIndexed() const {
		  return globalMatrixData_->isGloballyIndexed();
		}

		//! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
		bool isFillComplete() const {
		  return globalMatrixData_->isFillComplete();
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
		  globalMatrixData_->getLocalRowCopy(LocalRow, Indices, Values, NumEntries);
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
		   globalMatrixData_->getGlobalRowView(GlobalRow, indices, values);
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
		   globalMatrixData_->getLocalRowView(LocalRow, indices, values);
		}

		//! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
		/*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the
		  the zero and non-zero diagonals owned by this node. */
		void getLocalDiagCopy(Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const {
		  globalMatrixData_->getLocalDiagCopy(diag);
		}

		//! Get offsets of the diagonal entries in the matrix.
		void getLocalDiagOffsets(ArrayRCP<size_t> &offsets) const {
		  globalMatrixData_->getLocalDiagOffsets(offsets);
		}

		//! Get a copy of the diagonal entries owned by this node, with local row indices, using row offsets.
		void getLocalDiagCopy(Xpetra::Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &diag, const ArrayView<const size_t> &offsets) const {
		  globalMatrixData_->getLocalDiagCopy(diag,offsets);
		}

		//! Get Frobenius norm of the matrix
		typename ScalarTraits<Scalar>::magnitudeType getFrobeniusNorm() const {
		  return globalMatrixData_->getFrobeniusNorm();
		}

		//! Left scale matrix using the given vector entries
		void leftScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
		  globalMatrixData_->leftScale(x);
		}

		//! Right scale matrix using the given vector entries
		void rightScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
		  globalMatrixData_->rightScale(x);
		}

		//! Returns true if globalConstants have been computed; false otherwise
		bool haveGlobalConstants() const {
		  return globalMatrixData_->haveGlobalConstants();
		}

		//@}

		//! @name Methods implementing Matrix
		//@{

		//! \brief Computes the sparse matrix-multivector multiplication.
		/*! Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exceptions:
		- if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
		*/
		virtual void apply(const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
			   Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
			   Teuchos::ETransp mode = Teuchos::NO_TRANS,
			   Scalar alpha = ScalarTraits<Scalar>::one(),
			   Scalar beta = ScalarTraits<Scalar>::zero()) const {

			globalMatrixData_->apply(X,Y,mode,alpha,beta);
		}

		//! \brief Returns the Map associated with the domain of this operator.
		//! This will be <tt>null</tt> until fillComplete() is called.
		RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const {
		return globalMatrixData_->getDomainMap();
		}

		//! Returns the Map associated with the domain of this operator.
		//! This will be <tt>null</tt> until fillComplete() is called.
		RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const {
		return globalMatrixData_->getRangeMap();
		}

		//! \brief Returns the Map that describes the column distribution in this matrix.
		//! This might be <tt>null</tt> until fillComplete() is called.
		const RCP<const Map> & getColMap() const { return getColMap(Matrix::GetCurrentViewLabel()); }

		//! \brief Returns the Map that describes the column distribution in this matrix.
		const RCP<const Map> & getColMap(viewLabel_t viewLabel) const {
		TEUCHOS_TEST_FOR_EXCEPTION(Matrix::operatorViewTable_.containsKey(viewLabel) == false, Xpetra::Exceptions::RuntimeError, "Xpetra::Matrix.GetColMap(): view '" + viewLabel + "' does not exist.");
		updateDefaultView(); // If CrsMatrix::fillComplete() have been used instead of MatrixSplitting::fillComplete(), the default view is updated.
		return Matrix::operatorViewTable_.get(viewLabel)->GetColMap();
		}

		void removeEmptyProcessesInPlace(const RCP<const Map>& newMap) {
		globalMatrixData_->removeEmptyProcessesInPlace(newMap);
		this->operatorViewTable_.get(this->GetCurrentViewLabel())->SetRowMap(globalMatrixData_->getRowMap());
		this->operatorViewTable_.get(this->GetCurrentViewLabel())->SetColMap(globalMatrixData_->getColMap());
		}

		//@}

		//! Implements DistObject interface
		//{@

		//! Access function for the Tpetra::Map this DistObject was constructed with.
		const RCP< const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > getMap() const {
		return globalMatrixData_->getMap();
		}

		//! Import.
		void doImport(const Matrix &source,
			const Xpetra::Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM) {
			std::cout<<"Import not implemented"<<std::endl;
		//const MatrixSplitting & sourceWrp = dynamic_cast<const MatrixSplitting &>(source);
		//globalMatrixData_->doImport(*sourceWrp.getCrsMatrix(), importer, CM);
		}

		//! Export.
		void doExport(const Matrix &dest,
			const Xpetra::Import< LocalOrdinal, GlobalOrdinal, Node >& importer, CombineMode CM) {
			std::cout<<"Export not implemented"<<std::endl;
		//const MatrixSplitting & destWrp = dynamic_cast<const MatrixSplitting &>(dest);
		//globalMatrixData_->doExport(*destWrp.getCrsMatrix(), importer, CM);
		}

		//! Import (using an Exporter).
		void doImport(const Matrix &source,
			const Xpetra::Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM) {
			std::cout<<"Import not implemented"<<std::endl;
		//const MatrixSplitting & sourceWrp = dynamic_cast<const MatrixSplitting &>(source);
		//globalMatrixData_->doImport(*sourceWrp.getCrsMatrix(), exporter, CM);
		}

		//! Export (using an Importer).
		void doExport(const Matrix &dest,
			const Xpetra::Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM) {
			std::cout<<"Export not implemented"<<std::endl;
		//const MatrixSplitting & destWrp = dynamic_cast<const MatrixSplitting &>(dest);
		//globalMatrixData_->doExport(*destWrp.getCrsMatrix(), exporter, CM);
		}

		// @}


		//! @name Overridden from Teuchos::Describable
		//@{

		/** \brief Return a simple one-line description of this object. */
		std::string description() const {
		return "Xpetra::MatrixSplitting";
		}

		/** \brief Print the object with some verbosity level to an FancyOStream object. */
		void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const {
		//     Teuchos::EVerbosityLevel vl = verbLevel;
		//     if (vl == VERB_DEFAULT) vl = VERB_LOW;
		//     RCP<const Comm<int> > comm = this->getComm();
		//     const int myImageID = comm->getRank(),
		//       numImages = comm->getSize();

		//     if (myImageID == 0) out << this->description() << std::endl;

		globalMatrixData_->describe(out,verbLevel);
		}

		//! Returns the CrsGraph associated with this matrix.
		RCP<const CrsGraph> getCrsGraph() const { return globalMatrixData_->getCrsGraph(); }

		RCP<CrsMatrix> getCrsMatrix() const {  return globalMatrixData_; }


		//! Implements DistObject interface
		//{@
		void write(const char* output_file_name)
		{
                	Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write(output_file_name, *globalMatrixData_);
		}
		// @}
		
	private:

		// Default view is created after fillComplete()
		// Because ColMap might not be available before fillComplete().
		void CreateDefaultView() {

		  // Create default view
		  this->defaultViewLabel_ = "point";
		  this->CreateView(this->GetDefaultViewLabel(), globalMatrixData_->getRowMap(), globalMatrixData_->getColMap());

		  // Set current view
		  this->currentViewLabel_ = this->GetDefaultViewLabel();
		}

		// The colMap can be <tt>null</tt> until fillComplete() is called. The default view of the Matrix have to be updated when fillComplete() is called.
		// If CrsMatrix::fillComplete() have been used instead of MatrixSplitting::fillComplete(), the default view is updated when getColMap() is called.
		void updateDefaultView() const {
		  if ((finalDefaultView_ == false) &&  globalMatrixData_->isFillComplete() ) {
		    // Update default view with the colMap
		    Matrix::operatorViewTable_.get(Matrix::GetDefaultViewLabel())->SetColMap(globalMatrixData_->getColMap());
		    finalDefaultView_ = true;
		  }
		}

		//! Support Functions
		//{@
		void RegionalMatrix(GlobalOrdinal region_idx)
		{	
			TEUCHOS_TEST_FOR_EXCEPTION( num_total_regions_!=localMatrixData_.size(), Exceptions::RuntimeError, "Number of regions does not match with the size of localMatrixData_ structure \n");
			RCP<Matrix> regional_matrix = localMatrixData_[region_idx];
			Array< std::tuple<GlobalOrdinal,GlobalOrdinal> >   regionToAll = driver_->GetRegionToAll(region_idx);

			//THIS IS THE CORE OF THE PROBLEM WHERE ONE NEEDS TO POPULATE THE REGIONAL MATRICES BY ACCESSING ENTRIES OF THE GLOBAL MATRIX
			//
			/*regional_matrix -> resumeFill();
      			ArrayView<const GlobalOrdinal> MyGlobalElements = globalMatrixData_->getRowMap()->getNodeElementList();
 			for( typename ArrayView<const GlobalOrdinal>::iterator iter = MyGlobalElements.begin(); iter!=MyGlobalElements.end(); ++iter )	
			{
				checkerAllToRegion<GlobalOrdinal> unaryPredicate(*iter);
				typename Array< std::tuple<GlobalOrdinal, GlobalOrdinal > >::iterator region_iterator;
				region_iterator = std::find_if<typename Array< std::tuple< GlobalOrdinal,GlobalOrdinal > >::iterator, checkerAllToRegion<GlobalOrdinal> >(regionToAll.begin(), regionToAll.end(), unaryPredicate);
				if( region_iterator!=regionToAll.end() ) 
				{
					GlobalOrdinal node_regional_idx = std::get<0>( *region_iterator );
					ArrayView<const GlobalOrdinal> inds;
					ArrayView<const Scalar> vals;
					globalMatrixData_ -> getGlobalRowView( *iter, inds, vals );
				}
			}
			regional_matrix -> fillComplete();*/
		};
		// @}


		//! @name Creation of Local matrices
		//@{
		void CreateLocalMatrices( Array<Array<GlobalOrdinal> > local_maps){ 

			TEUCHOS_TEST_FOR_EXCEPTION( num_total_regions_!=local_maps.size(), Exceptions::RuntimeError, "Number of regions does not match with the size of local_maps structure \n");

			localMatrixData_.clear( );		

			for( int i = 0; i<num_total_regions_; ++i )
			{
				//Create Xpetra map for local stiffness matrix
				RCP<const Xpetra::Map<int,GlobalOrdinal,Node> > xpetraMap;
				xpetraMap = Xpetra::MapFactory<int,GlobalOrdinal,Node>::Build(lib, driver_->num_regional_nodes_[i], local_maps[i], 0, comm_); 			
				int num_elements = xpetraMap->getGlobalNumElements();
				
				RCP<CrsMatrix> crs_matrix;
				if( Xpetra::UseEpetra==lib )
					crs_matrix = rcp( new EpetraCrsMatrix(xpetraMap, num_elements) );
				else if( Xpetra::UseTpetra==lib )
					crs_matrix = rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>( xpetraMap, num_elements ) );	
				else
					std::cerr<<" The library to build matrices must be either Epetra or Tpetra \n";

				RCP<Matrix> matrixPointer = rcp(new CrsMatrixWrap(crs_matrix));
				localMatrixData_.push_back( matrixPointer );
			}

			for( GlobalOrdinal i = 0; i<num_total_regions_; ++i )
				RegionalMatrix(i);

		};
		//@}

		// The boolean finalDefaultView_ keep track of the status of the default view (= already updated or not)
		// See also MatrixSplitting::updateDefaultView()
		mutable bool finalDefaultView_;

		RCP<const Teuchos::Comm<int> > comm_;

		RCP<SplittingDriver<Scalar, LocalOrdinal, GlobalOrdinal, Node> > driver_;
		RCP<Matrix> globalMatrixData_;
		Array<RCP<Matrix> > localMatrixData_;

		GlobalOrdinal num_total_elements_;
		GlobalOrdinal num_total_regions_;

	}; //class MatrixSplitting

} //namespace Xpetra

#define XPETRA_MATRIX_SHORT
#endif //XPETRA_MATRIXSPLITTING_HPP
