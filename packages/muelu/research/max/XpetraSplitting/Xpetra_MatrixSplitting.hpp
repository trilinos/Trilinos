// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * Xpetra_MatrixSplitting.hpp
 *
 * Created on: August 17, 2017
 * 	Author: Massimiliano Lupo Pasini (massimiliano.lupo.pasini@gmail.com)
 */
#ifndef XPETRA_MATRIXSPLITTING_HPP
#define XPETRA_MATRIXSPLITTING_HPP

// Xpetra
#include "Xpetra_Map.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_Matrix.hpp"
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include "Xpetra_IO.hpp"
#include "Xpetra_RegionHandler_def.hpp"

// Ifpack2
#include "Ifpack2_OverlappingRowMatrix_def.hpp"

// MueLu
#include <MueLu_Utilities.hpp>

/** \file Xpetra_MatrixSplitting.hpp

Declarations for the class Xpetra::MatrixSplitting.
 */
namespace Xpetra {

// Definition of the predicate for the regionToAll structure.
// Given a tuple made of node index and a specific region it belongs to,
// this predicate returns true if the node has composite index which coincides with the index specified in input.
template <class GlobalOrdinal>
class checkerRegionToAll {
 public:
  // Constructor
  checkerRegionToAll(GlobalOrdinal node_index) { node_index_ = node_index; };

  // Unary Operator
  bool operator()(const std::tuple<GlobalOrdinal, GlobalOrdinal> &node) { return (std::get<0>(node) == node_index_); }

 private:
  GlobalOrdinal node_index_;
};

// Definition of the predicate for the node_ structure.
// Given a tuple made of node index and a specific region it belongs to,
// this predicate returns true if the node has composite index which coincides with the index specified in input.
// It does the same thing as checkerRegionToAll but it works on a different data structure
template <class GlobalOrdinal>
class checkerAllToRegion {
 public:
  // Constructor
  checkerAllToRegion(GlobalOrdinal node_index) { node_index_ = node_index; };

  // Unary Operator
  bool operator()(const std::tuple<GlobalOrdinal, GlobalOrdinal> &node) { return (std::get<1>(node) == node_index_); }

 private:
  GlobalOrdinal node_index_;
};

// Definition of the predicate for the node_ structure.
// Given a tuple made of node index and a specific region it belongs to,
// this predicate returns true if the node has composite index which coincides with the index specified in input.
// This checker is specifically used only for nodes lying on the interface
template <class GlobalOrdinal>
class checkerInterfaceNodes {
 public:
  // Constructor
  checkerInterfaceNodes(GlobalOrdinal node_index) { node_index_ = node_index; };

  // Unary Operator
  bool operator()(const std::tuple<int, Array<GlobalOrdinal> > &node) { return (std::get<0>(node) == node_index_); }

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
          UnderlyingLib lib   = Xpetra::UseEpetra,
          bool collapse       = false>
class MatrixSplitting : public Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;
  typedef Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> CrsGraph;
  typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrix;
  typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrixWrap;
  typedef Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrixFactory;
  typedef Xpetra::MatrixView<Scalar, LocalOrdinal, GlobalOrdinal, Node> MatrixView;

  // Xpetra structures must be converted into Tpetra specialized ones to construct an Ifpack2::OverlappingRowMatrix object
  // Once the Ifpack2::OverlappingRowMatrix class is transferred into the Xpetra directory, the following 6 lines can be changed/removed
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> tpetra_crs_matrix;
  typedef Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> tpetra_row_matrix;

 public:
  //! @name Constructor/Destructor Methods
  //@{

  //! Constructor specifying fixed number of entries for each row.
  MatrixSplitting(RCP<Matrix> matrix, RCP<Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > > nodes) {
    std::cout << "This version of MatrixSplitting constructor is NOT currently supported \n";
  }
  //
  //
  MatrixSplitting(const char *matrix_file_name,
                  Teuchos::RCP<Xpetra::RegionHandler<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionHandler,
                  RCP<const Teuchos::Comm<int> > comm) {
    comm_                            = comm;
    regionHandler_                   = regionHandler;
    Array<GlobalOrdinal> elementlist = regionHandler_->GetGlobalRowMap();
    num_total_elements_              = regionHandler_->GetNumGlobalElements();
    num_total_regions_               = regionHandler_->GetNumTotalRegions();

    if (comm_->getRank() == 0)
      std::cout << "MatrixSplitting constructor initialized" << std::endl;
    region_matrix_initialized_.clear();
    for (int i = 0; i < num_total_regions_; ++i)
      region_matrix_initialized_.push_back(false);

    // Create Xpetra map for composite stiffness matrix
    if (comm_->getRank() == 0)
      std::cout << "Starting construction of Composite Map" << std::endl;
    RCP<const Xpetra::Map<int, GlobalOrdinal, Node> > xpetraMap;
    xpetraMap = Xpetra::MapFactory<int, GlobalOrdinal, Node>::Build(lib, num_total_elements_, elementlist, 0, comm);
    if (comm_->getRank() == 0)
      std::cout << "Finished construction of Composite Map" << std::endl;

    if (comm_->getRank() == 0)
      std::cout << "Started reading composite matrix" << std::endl;
    // Import matrix from an .mm file into an Xpetra wrapper for an Epetra matrix
    compositeMatrixData_ = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read(matrix_file_name, xpetraMap);
    if (comm_->getRank() == 0)
      std::cout << "Finished reading composite matrix" << std::endl;

    CreateRegionMatrices(regionHandler_->GetRegionRowMaps());
  }

  //! Destructor
  virtual ~MatrixSplitting() {}
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
    compositeMatrixData_->insertGlobalValues(globalRow, cols, vals);
  }

  //! Insert matrix entries, using local IDs.
  /** All index values must be in the local space.
        \pre \c localRow exists as an ID in the globalrow map
        \pre <tt>isGloballyIndexed() == false</tt>
        \pre <tt>isStorageOptimized() == false</tt>

        \post <tt>isLocallyIndexed() == true</tt>
   */
  void insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &cols, const ArrayView<const Scalar> &vals) {
    compositeMatrixData_->insertLocalValues(localRow, cols, vals);
  }

  //! \brief Replace matrix entries, using global IDs.
  /** All index values must be in the global space.

    \pre \c globalRow is a global row belonging to the matrix on this node.

    \note If (globalRow,cols[i]) corresponds to an entry that is duplicated in this matrix row (likely because it was inserted more than once and fillComplete() has not been called in the interim), the behavior of this function is not defined. */
  void replaceGlobalValues(GlobalOrdinal globalRow,
                           const ArrayView<const GlobalOrdinal> &cols,
                           const ArrayView<const Scalar> &vals) { compositeMatrixData_->replaceGlobalValues(globalRow, cols, vals); }

  //! Replace matrix entries, using local IDs.
  /** All index values must be in the local space.
        Note that if a value is not already present for the specified location in the matrix, the input value will be ignored silently.
   */
  void replaceLocalValues(LocalOrdinal localRow,
                          const ArrayView<const LocalOrdinal> &cols,
                          const ArrayView<const Scalar> &vals) { compositeMatrixData_->replaceLocalValues(localRow, cols, vals); }

  //! Set all matrix entries equal to scalar
  virtual void setAllToScalar(const Scalar &alpha) { compositeMatrixData_->setAllToScalar(alpha); }

  //! Scale the current values of a matrix, this = alpha*this.
  void scale(const Scalar &alpha) {
    compositeMatrixData_->scale(alpha);
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
  void resumeFill(const RCP<ParameterList> &params = null) {
    compositeMatrixData_->resumeFill(params);
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
    compositeMatrixData_->fillComplete(domainMap, rangeMap, params);

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
  // TODO : Get ride of "Tpetra"::OptimizeOption
  void fillComplete(const RCP<ParameterList> &params = null) {
    compositeMatrixData_->fillComplete(params);

    // Update default view with the colMap because colMap can be <tt>null</tt> until fillComplete() is called.
    updateDefaultView();
  }

  //@}

  //! Returns the number of global rows in this matrix.
  /** Undefined if isFillActive().
   */
  global_size_t getGlobalNumRows() const {
    return compositeMatrixData_->getGlobalNumRows();
  }

  //! \brief Returns the number of global columns in the matrix.
  /** Undefined if isFillActive().
   */
  global_size_t getGlobalNumCols() const {
    return compositeMatrixData_->getGlobalNumCols();
  }

  //! Returns the number of matrix rows owned on the calling node.
  size_t getLocalNumRows() const {
    return compositeMatrixData_->getLocalNumRows();
  }

  //! Returns the global number of entries in this matrix.
  global_size_t getGlobalNumEntries() const {
    return compositeMatrixData_->getGlobalNumEntries();
  }

  //! Returns the local number of entries in this matrix.
  size_t getLocalNumEntries() const {
    return compositeMatrixData_->getLocalNumEntries();
  }

  //! Returns the current number of entries on this node in the specified local row.
  /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
  size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const {
    return compositeMatrixData_->getNumEntriesInLocalRow(localRow);
  }

  //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
  /** Undefined if isFillActive().
   */
  size_t getGlobalMaxNumRowEntries() const {
    return compositeMatrixData_->getGlobalMaxNumRowEntries();
  }

  //! \brief Returns the maximum number of entries across all rows/columns on this node.
  /** Undefined if isFillActive().
   */
  size_t getLocalMaxNumRowEntries() const {
    return compositeMatrixData_->getLocalMaxNumRowEntries();
  }

  //! \brief Returns the number of regions in the composite domain.
  global_size_t getNumRegions() const {
    TEUCHOS_TEST_FOR_EXCEPTION(region_matrix_initialized_.size() == 0, Exceptions::RuntimeError, "Regions have not been initialized yet \n");
    for (int i = 0; i < region_matrix_initialized_.size(); ++i)
      TEUCHOS_TEST_FOR_EXCEPTION(!region_matrix_initialized_[i], Exceptions::RuntimeError, "Region matrix for region " << i + 1 << "has not been initialized yet \n");

    return num_total_regions_;
  }

  //! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
  bool isLocallyIndexed() const {
    return compositeMatrixData_->isLocallyIndexed();
  }

  //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
  bool isGloballyIndexed() const {
    return compositeMatrixData_->isGloballyIndexed();
  }

  //! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
  bool isFillComplete() const {
    return compositeMatrixData_->isFillComplete();
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
                       size_t &NumEntries) const {
    compositeMatrixData_->getLocalRowCopy(LocalRow, Indices, Values, NumEntries);
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
    compositeMatrixData_->getGlobalRowView(GlobalRow, indices, values);
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
    compositeMatrixData_->getLocalRowView(LocalRow, indices, values);
  }

  //! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
  /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the
                  the zero and non-zero diagonals owned by this node. */
  void getLocalDiagCopy(Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &diag) const {
    compositeMatrixData_->getLocalDiagCopy(diag);
  }

  //! Get offsets of the diagonal entries in the matrix.
  void getLocalDiagOffsets(ArrayRCP<size_t> &offsets) const {
    compositeMatrixData_->getLocalDiagOffsets(offsets);
  }

  //! Get a copy of the diagonal entries owned by this node, with local row indices, using row offsets.
  void getLocalDiagCopy(Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &diag, const ArrayView<const size_t> &offsets) const {
    compositeMatrixData_->getLocalDiagCopy(diag, offsets);
  }

  //! Get Frobenius norm of the matrix
  typename ScalarTraits<Scalar>::magnitudeType getFrobeniusNorm() const {
    return compositeMatrixData_->getFrobeniusNorm();
  }

  //! Left scale matrix using the given vector entries
  void leftScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &x) {
    compositeMatrixData_->leftScale(x);
  }

  //! Neighbor2 scale matrix using the given vector entries
  void rightScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &x) {
    compositeMatrixData_->rightScale(x);
  }

  //! Returns true if compositeConstants have been computed; false otherwise
  bool haveGlobalConstants() const {
    return compositeMatrixData_->haveGlobalConstants();
  }

  //@}

  //! @name Methods implementing Matrix
  //@{

  //! \brief Computes the sparse matrix-multivector multiplication.
  /*! Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exceptions:
                - if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
   */
  virtual void apply(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &X,
                     Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &Y,
                     Teuchos::ETransp mode = Teuchos::NO_TRANS,
                     Scalar alpha          = ScalarTraits<Scalar>::one(),
                     Scalar beta           = ScalarTraits<Scalar>::zero()) const {
    compositeMatrixData_->apply(X, Y, mode, alpha, beta);
  }

  //! \brief Returns the Map associated with the domain of this operator.
  //! This will be <tt>null</tt> until fillComplete() is called.
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getDomainMap() const {
    return compositeMatrixData_->getDomainMap();
  }

  //! Returns the Map associated with the domain of this operator.
  //! This will be <tt>null</tt> until fillComplete() is called.
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getRangeMap() const {
    return compositeMatrixData_->getRangeMap();
  }

  //! \brief Returns the Map that describes the column distribution in this matrix.
  //! This might be <tt>null</tt> until fillComplete() is called.
  const RCP<const Map> &getColMap() const { return getColMap(Matrix::GetCurrentViewLabel()); }

  //! \brief Returns the Map that describes the column distribution in this matrix.
  const RCP<const Map> &getColMap(viewLabel_t viewLabel) const {
    TEUCHOS_TEST_FOR_EXCEPTION(Matrix::operatorViewTable_.containsKey(viewLabel) == false, Xpetra::Exceptions::RuntimeError, "Xpetra::Matrix.GetColMap(): view '" + viewLabel + "' does not exist.");
    updateDefaultView();  // If CrsMatrix::fillComplete() have been used instead of MatrixSplitting::fillComplete(), the default view is updated.
    return Matrix::operatorViewTable_.get(viewLabel)->GetColMap();
  }

  void removeEmptyProcessesInPlace(const RCP<const Map> &newMap) {
    compositeMatrixData_->removeEmptyProcessesInPlace(newMap);
    this->operatorViewTable_.get(this->GetCurrentViewLabel())->SetRowMap(compositeMatrixData_->getRowMap());
    this->operatorViewTable_.get(this->GetCurrentViewLabel())->SetColMap(compositeMatrixData_->getColMap());
  }

  //@}

  //! Implements DistObject interface
  //{@

  //! Access function for the Tpetra::Map this DistObject was constructed with.
  const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getMap() const {
    return compositeMatrixData_->getMap();
  }

  //! Import.
  void doImport(const Matrix &source,
                const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> &importer, CombineMode CM) {
    std::cout << "Import not implemented" << std::endl;
    // const MatrixSplitting & sourceWrp = dynamic_cast<const MatrixSplitting &>(source);
    // compositeMatrixData_->doImport(*sourceWrp.getCrsMatrix(), importer, CM);
  }

  //! Export.
  void doExport(const Matrix &dest,
                const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> &importer, CombineMode CM) {
    std::cout << "Export not implemented" << std::endl;
    // const MatrixSplitting & destWrp = dynamic_cast<const MatrixSplitting &>(dest);
    // compositeMatrixData_->doExport(*destWrp.getCrsMatrix(), importer, CM);
  }

  //! Import (using an Exporter).
  void doImport(const Matrix &source,
                const Xpetra::Export<LocalOrdinal, GlobalOrdinal, Node> &exporter, CombineMode CM) {
    std::cout << "Import not implemented" << std::endl;
    // const MatrixSplitting & sourceWrp = dynamic_cast<const MatrixSplitting &>(source);
    // compositeMatrixData_->doImport(*sourceWrp.getCrsMatrix(), exporter, CM);
  }

  //! Export (using an Importer).
  void doExport(const Matrix &dest,
                const Xpetra::Export<LocalOrdinal, GlobalOrdinal, Node> &exporter, CombineMode CM) {
    std::cout << "Export not implemented" << std::endl;
    // const MatrixSplitting & destWrp = dynamic_cast<const MatrixSplitting &>(dest);
    // compositeMatrixData_->doExport(*destWrp.getCrsMatrix(), exporter, CM);
  }

  // @}

  //! @name Overridden from Teuchos::Describable
  //@{

  /** \brief Return a simple one-line description of this object. */
  std::string description() const {
    return "Xpetra::MatrixSplitting";
  }

  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const {
    compositeMatrixData_->describe(out, verbLevel);
  }
  //@}

  //! Get methods
  //@{

  //! Returns the CrsGraph associated with this matrix.
  RCP<const CrsGraph> getCrsGraph() const { return compositeMatrixData_->getCrsGraph(); }

  //! Returns an Xpetra::CrsMatrix pointer to the the composite matrix
  RCP<CrsMatrix> getCrsMatrix() const { return compositeMatrixData_; }

  //! Returns an Xpetra::Matrix pointer to the composite matrix
  RCP<Matrix> getMatrix() const { return compositeMatrixData_; }

  //! Returns an Xpetra::Matrix pointer to the region matrix associated with a specific region index
  RCP<Matrix> getRegionMatrix(GlobalOrdinal region_idx) const {
    // The region index is assumed to start from 0
    TEUCHOS_TEST_FOR_EXCEPTION(num_total_regions_ <= 0, Exceptions::RuntimeError, "Regions not initialized yet ( total number of regions is <=0 ) \n");
    TEUCHOS_TEST_FOR_EXCEPTION(region_idx >= num_total_regions_, Exceptions::RuntimeError, "Region index not valid \n");

    return regionMatrixData_[region_idx];
  }

  //! Return a ppointer to the underlying regionHandler object used for the matrix splitting
  RCP<RegionHandler<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getRegionHandler() const {
    return regionHandler_;
  }
  //@}

  //! Write methods
  //{@
  void writeGlobalMatrix() {
    std::string file_name;
    file_name += "./output/A_composite.mm";
    Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write(file_name, *compositeMatrixData_);
  }

  void writeRegionMatrices() {
    for (int i = 0; i < num_total_regions_; ++i) {
      std::string region_str = std::to_string(i + 1);
      std::string file_name;
      file_name += "./output/A_region_";
      file_name += region_str;
      file_name += ".mm";
      Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write(file_name.c_str(), *regionMatrixData_[i]);
    }
  }
  // @}

 private:
  // Default view is created after fillComplete()
  // Because ColMap might not be available before fillComplete().
  void CreateDefaultView() {
    // Create default view
    this->defaultViewLabel_ = "point";
    this->CreateView(this->GetDefaultViewLabel(), compositeMatrixData_->getRowMap(), compositeMatrixData_->getColMap());

    // Set current view
    this->currentViewLabel_ = this->GetDefaultViewLabel();
  }

  // The colMap can be <tt>null</tt> until fillComplete() is called. The default view of the Matrix have to be updated when fillComplete() is called.
  // If CrsMatrix::fillComplete() have been used instead of MatrixSplitting::fillComplete(), the default view is updated when getColMap() is called.
  void updateDefaultView() const {
    if ((finalDefaultView_ == false) && compositeMatrixData_->isFillComplete()) {
      // Update default view with the colMap
      Matrix::operatorViewTable_.get(Matrix::GetDefaultViewLabel())->SetColMap(compositeMatrixData_->getColMap());
      finalDefaultView_ = true;
    }
  }

  //! Private Methods
  //{@

  //! Creation of region matrices
  //@{
  void RegionMatrix(GlobalOrdinal region_idx) {
    TEUCHOS_TEST_FOR_EXCEPTION(num_total_regions_ != regionMatrixData_.size(), Exceptions::RuntimeError, "Number of regions does not match with the size of regionMatrixData_ structure \n");
    RCP<Matrix> region_matrix = regionMatrixData_[region_idx];

    RCP<tpetra_crs_matrix> tpetraGlobalMatrix = MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraCrs(compositeMatrixData_);
    Ifpack2::OverlappingRowMatrix<tpetra_row_matrix> enlargedMatrix(tpetraGlobalMatrix, 2);

    region_matrix->resumeFill();

    // Region matrices are initially built to be a chopped version of the composite matrix
    InitializeRegionMatrices(region_idx, region_matrix, enlargedMatrix);

    // If the template paramater is set to collapse by the user, then interface entries of the region matrix are modified to collapse
    // information coming from adjacent regions. If the collapsing is not done, then the splitting is calculated
    if (collapse)
      RegionCollapse(region_idx, region_matrix, enlargedMatrix);
    else
      RegionSplitting(region_idx, region_matrix, enlargedMatrix);

    region_matrix->fillComplete();
  };
  // @}

  //! @name Initialization of Region matrices
  //@{
  void InitializeRegionMatrices(GlobalOrdinal region_idx, RCP<Matrix> &region_matrix, Ifpack2::OverlappingRowMatrix<tpetra_row_matrix> &enlargedMatrix) {
    // Region matrices are initially built to be a chopped version of the composite matrix
    TEUCHOS_TEST_FOR_EXCEPTION(region_matrix_initialized_[region_idx], Exceptions::RuntimeError, "Surrogate region stiffness matrices are already initialized by chopping the composite stiffness matrix \n");

    Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > regionToAll = regionHandler_->GetRegionToAll(region_idx);

    // THIS IS THE CORE OF THE PROBLEM WHERE ONE NEEDS TO POPULATE THE REGIONAL MATRICES BY ACCESSING ENTRIES OF THE GLOBAL MATRIX
    //
    ArrayView<const GlobalOrdinal> MyRegionElements = region_matrix->getRowMap()->getLocalElementList();
    for (typename ArrayView<const GlobalOrdinal>::iterator iter = MyRegionElements.begin(); iter != MyRegionElements.end(); ++iter) {
      // Nodes are saved in data structures with 1 as base index
      checkerRegionToAll<GlobalOrdinal> unaryPredicate(*iter + 1);
      typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator composite_iterator;
      composite_iterator = std::find_if<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerRegionToAll<GlobalOrdinal> >(regionToAll.begin(), regionToAll.end(), unaryPredicate);
      TEUCHOS_TEST_FOR_EXCEPTION(composite_iterator == regionToAll.end(), Exceptions::RuntimeError, "Process ID: " << comm_->getRank() << " - Region: " << region_idx << " - "
                                                                                                                   << " node with region index: " << *iter + 1 << " is not in regionToAll[" << region_idx << "]"
                                                                                                                   << "\n");
      GlobalOrdinal node_idx      = std::get<1>(*composite_iterator);
      LocalOrdinal node_local_idx = enlargedMatrix.getRowMap()->getLocalElement(node_idx - 1);

      ArrayView<const GlobalOrdinal> inds;
      ArrayView<const Scalar> vals;
      enlargedMatrix.getLocalRowView(node_local_idx, inds, vals);

      std::vector<GlobalOrdinal> region_inds_vector(0);
      std::vector<Scalar> region_vals_vector(0);

      for (LocalOrdinal i = 0; i < inds.size(); ++i) {
        // Nodes are saved in data structures with 1 as base index
        GlobalOrdinal composite_col_ind = enlargedMatrix.getColMap()->getGlobalElement(inds[i]) + 1;
        checkerAllToRegion<GlobalOrdinal> unaryPredicate2(composite_col_ind);
        typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator region_iterator;
        region_iterator = std::find_if<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerAllToRegion<GlobalOrdinal> >(regionToAll.begin(), regionToAll.end(), unaryPredicate2);
        if (region_iterator != regionToAll.end()) {
          region_inds_vector.push_back(std::get<0>(*region_iterator) - 1);
          region_vals_vector.push_back(vals[i]);
        }
      }
      ArrayView<GlobalOrdinal> region_inds(region_inds_vector);
      ArrayView<Scalar> region_vals(region_vals_vector);
      region_matrix->insertGlobalValues(*iter, region_inds, region_vals);
    }
    region_matrix_initialized_[region_idx] = true;
  }
  //@}

  //! @name Collapse of External neighbouring nodes information on Region matrices
  //@{
  void RegionCollapse(GlobalOrdinal region_idx, RCP<Matrix> &region_matrix, Ifpack2::OverlappingRowMatrix<tpetra_row_matrix> &enlargedMatrix) {
    TEUCHOS_TEST_FOR_EXCEPTION(!region_matrix_initialized_[region_idx], Exceptions::RuntimeError, "The composite stiffness matrix must be chopped into surrogate region matrices before collapsing \n");
    TEUCHOS_TEST_FOR_EXCEPTION(regionHandler_->GetNumRegionNodes(region_idx) != regionMatrixData_[region_idx]->getGlobalNumRows(), Exceptions::RuntimeError, "Process ID: " << comm_->getRank() << " - Number of region nodes in region " << region_idx + 1 << " does not coincide with the value returned by regionMatrixData_[" << region_idx + 1 << "]->getGlobalNumRows() \n");

    Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > regionToAll = regionHandler_->GetRegionToAll(region_idx);

    // This portion of the code assumes that the number of region nodes is the same on each direction of the domain
    // Foir a 2D problem we have then nx = ny = sqrt( num_region_ndoes_ )
    GlobalOrdinal n;
    GlobalOrdinal nx;
    GlobalOrdinal ny;

    n  = regionHandler_->GetNumRegionNodes(region_idx);
    nx = std::sqrt(n);
    ny = nx;
    TEUCHOS_TEST_FOR_EXCEPTION(static_cast<double>(nx - std::floor(static_cast<double>(std::sqrt(static_cast<double>(n))))) != 0.0, Exceptions::RuntimeError, "The code assumes that the regions are 2D and that the number of region nodes is the same on each direction of the domain \n");

    // interfaceNodes contains nodes on an interface between any regions
    Array<std::tuple<int, Array<GlobalOrdinal> > > interfaceNodes = regionHandler_->GetInterfaceNodes();

    ArrayView<const GlobalOrdinal> MyRegionElements = region_matrix->getRowMap()->getLocalElementList();
    for (typename ArrayView<const GlobalOrdinal>::iterator iter = MyRegionElements.begin(); iter != MyRegionElements.end(); ++iter) {
      // Nodes are saved in data structures with 1 as base index
      GlobalOrdinal region_node_idx = *iter + 1;
      checkerRegionToAll<GlobalOrdinal> unaryPredicate(region_node_idx);
      typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator composite_iterator;
      composite_iterator = std::find_if<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerRegionToAll<GlobalOrdinal> >(regionToAll.begin(), regionToAll.end(), unaryPredicate);
      TEUCHOS_TEST_FOR_EXCEPTION(composite_iterator == regionToAll.end(), Exceptions::RuntimeError, "Process ID: " << comm_->getRank() << " - Region: " << region_idx << " - "
                                                                                                                   << " node with region index: " << region_node_idx << " is not in regionToAll[" << region_idx << "]"
                                                                                                                   << "\n");

      GlobalOrdinal composite_node_idx = std::get<1>(*composite_iterator);
      checkerInterfaceNodes<GlobalOrdinal> unaryPredicate2(composite_node_idx);
      typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator interface_iterator;
      interface_iterator = std::find_if<typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator, checkerInterfaceNodes<GlobalOrdinal> >(interfaceNodes.begin(), interfaceNodes.end(), unaryPredicate2);

      // Here we assuming that a specific labeling choice is adopted region wise and we use it to distinguish coarse node from fine nodes
      bool coarse_point = false;
      if (region_node_idx % 3 == 1)
        coarse_point = true;

      GlobalOrdinal region_node_idx_neighbor1 = 0;
      GlobalOrdinal region_node_idx_neighbor2 = 0;

      // Horizontal-Vertical Collapse
      if (interface_iterator != interfaceNodes.end() && region_node_idx > ny && region_node_idx <= (nx - 1) * ny && !coarse_point) {
        region_node_idx_neighbor1 = region_node_idx - ny;
        region_node_idx_neighbor2 = region_node_idx + ny;
      } else if (interface_iterator != interfaceNodes.end() && region_node_idx % ny > 1 && !coarse_point) {
        region_node_idx_neighbor1 = region_node_idx - 1;
        region_node_idx_neighbor2 = region_node_idx + 1;
      }

      if (region_node_idx_neighbor1 != 0 && region_node_idx_neighbor2 != 0) {
        // Computation of composite index for neighbor1 node
        checkerRegionToAll<GlobalOrdinal> unaryPredicateLeft(region_node_idx_neighbor1);
        typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator composite_iterator_neighbor1;
        composite_iterator_neighbor1 = std::find_if<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerRegionToAll<GlobalOrdinal> >(regionToAll.begin(), regionToAll.end(), unaryPredicateLeft);

        // Computation of composite index for neighbor2 node
        checkerRegionToAll<GlobalOrdinal> unaryPredicateRight(region_node_idx_neighbor2);
        typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator composite_iterator_neighbor2;
        composite_iterator_neighbor2 = std::find_if<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerRegionToAll<GlobalOrdinal> >(regionToAll.begin(), regionToAll.end(), unaryPredicateRight);

        TEUCHOS_TEST_FOR_EXCEPTION(composite_iterator_neighbor1 == regionToAll.end() || composite_iterator_neighbor2 == regionToAll.end(), Exceptions::RuntimeError, "Process ID: " << comm_->getRank() << " - Region: " << region_idx << " - "
                                                                                                                                                                                    << " node with region index: " << region_node_idx << " lies on the interface between regions: " << std::get<1>(*interface_iterator) << " BUT has compositely mislabeled neighbouring nodes missing from regionToAll \n");

        // Check to see if neighbor1 node lies on a coarse line
        GlobalOrdinal composite_node_idx_neighbor1 = std::get<1>(*composite_iterator_neighbor1);
        checkerInterfaceNodes<GlobalOrdinal> unaryPredicate2neighbor1(composite_node_idx_neighbor1);
        typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator interface_iterator_neighbor1;
        interface_iterator_neighbor1 = std::find_if<typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator, checkerInterfaceNodes<GlobalOrdinal> >(interfaceNodes.begin(), interfaceNodes.end(), unaryPredicate2neighbor1);

        // Check to see if neighbor2 node lies on a coarse line
        GlobalOrdinal composite_node_idx_neighbor2 = std::get<1>(*composite_iterator_neighbor2);
        checkerInterfaceNodes<GlobalOrdinal> unaryPredicate2neighbor2(composite_node_idx_neighbor2);
        typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator interface_iterator_neighbor2;
        interface_iterator_neighbor2 = std::find_if<typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator, checkerInterfaceNodes<GlobalOrdinal> >(interfaceNodes.begin(), interfaceNodes.end(), unaryPredicate2neighbor2);

        // I apply the collapse only if the current node is a fine node which lies on a coarse line
        // This means that the neighbor1 node and neighbor2 node must both lie on the coarse line as well
        if (interface_iterator_neighbor1 != interfaceNodes.end() && interface_iterator_neighbor2 != interfaceNodes.end()) {
          // For each fine node on a horixontal coarse line on the interface, I extract the rows from the composite matrix
          LocalOrdinal node_idx           = enlargedMatrix.getRowMap()->getLocalElement(composite_node_idx - 1);
          LocalOrdinal node_idx_neighbor1 = enlargedMatrix.getRowMap()->getLocalElement(composite_node_idx_neighbor1 - 1);
          LocalOrdinal node_idx_neighbor2 = enlargedMatrix.getRowMap()->getLocalElement(composite_node_idx_neighbor2 - 1);
          ArrayView<const LocalOrdinal> inds;
          ArrayView<const Scalar> vals;
          enlargedMatrix.getLocalRowView(node_idx, inds, vals);
          ArrayView<const LocalOrdinal> inds_neighbor1;
          ArrayView<const Scalar> vals_neighbor1;
          enlargedMatrix.getLocalRowView(node_idx_neighbor1, inds_neighbor1, vals_neighbor1);
          ArrayView<const LocalOrdinal> inds_neighbor2;
          ArrayView<const Scalar> vals_neighbor2;
          enlargedMatrix.getLocalRowView(node_idx_neighbor2, inds_neighbor2, vals_neighbor2);

          std::vector<LocalOrdinal> inds_vector           = createVector(inds);
          std::vector<LocalOrdinal> inds_neighbor1_vector = createVector(inds_neighbor1);
          std::vector<LocalOrdinal> inds_neighbor2_vector = createVector(inds_neighbor2);

          std::vector<GlobalOrdinal> composite_inds_vector(0);
          std::vector<GlobalOrdinal> composite_inds_neighbor1_vector(0);
          std::vector<GlobalOrdinal> composite_inds_neighbor2_vector(0);

          for (typename std::vector<LocalOrdinal>::iterator iter_node = inds_vector.begin(); iter_node != inds_vector.end(); ++iter_node)
            composite_inds_vector.push_back(enlargedMatrix.getRowMap()->getGlobalElement(*iter_node));
          std::sort(composite_inds_vector.begin(), composite_inds_vector.end());

          for (typename std::vector<LocalOrdinal>::iterator iter_node = inds_neighbor1_vector.begin(); iter_node != inds_neighbor1_vector.end(); ++iter_node)
            composite_inds_neighbor1_vector.push_back(enlargedMatrix.getRowMap()->getGlobalElement(*iter_node));

          std::sort(composite_inds_neighbor1_vector.begin(), composite_inds_neighbor1_vector.end());

          for (typename std::vector<LocalOrdinal>::iterator iter_node = inds_neighbor2_vector.begin(); iter_node != inds_neighbor2_vector.end(); ++iter_node)
            composite_inds_neighbor2_vector.push_back(enlargedMatrix.getRowMap()->getGlobalElement(*iter_node));

          std::sort(composite_inds_neighbor2_vector.begin(), composite_inds_neighbor2_vector.end());

          // IDENTIFICATION OF EXTERNAL NODES THROUGH COMPOSITE INDICES STARTS HERE
          std::vector<GlobalOrdinal> composite_node_idx_neighbor1_extra;
          std::vector<GlobalOrdinal> composite_node_idx_neighbor2_extra;
          std::vector<GlobalOrdinal> composite_node_idx_extra(composite_inds_vector);

          // The follolwing triple of vector is expected to EVENTUALLY contain only one entry:
          // the label of the external node with information to collapse close to neighbor1, neighbor2 and central node
          std::vector<GlobalOrdinal> diff_neighbor1;
          std::vector<GlobalOrdinal> diff_neighbor2;
          std::vector<GlobalOrdinal> diff_center;

          // Identification of external node from the side of neighbor1
          {
            // Compute the intersection between neoghbourhood of neighbor1 node and neighbourhood of central node
            std::set_intersection(composite_inds_vector.begin(), composite_inds_vector.end(), composite_inds_neighbor1_vector.begin(), composite_inds_neighbor1_vector.end(), std::back_inserter(composite_node_idx_neighbor1_extra));
            for (typename std::vector<GlobalOrdinal>::iterator iter_node = composite_node_idx_neighbor1_extra.begin(); iter_node != composite_node_idx_neighbor1_extra.end(); ++iter_node) {
              checkerAllToRegion<GlobalOrdinal> unaryPredicateExtra(*iter_node + 1);
              typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator region_iterator_extra;
              region_iterator_extra = std::find_if<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerAllToRegion<GlobalOrdinal> >(regionToAll.begin(), regionToAll.end(), unaryPredicateExtra);

              // Invalidation of node indices for nodes belonging to the current region
              if (region_iterator_extra != regionToAll.end())
                *iter_node = -1;
            }

            // Removal of invalidated indices associated with nodes belonging to current region: (external nodes do not belong to this region)
            composite_node_idx_neighbor1_extra.erase(std::remove(composite_node_idx_neighbor1_extra.begin(), composite_node_idx_neighbor1_extra.end(), -1), composite_node_idx_neighbor1_extra.end());

            // External node from neighbor1 side does not belong to the neighborhood of neighbor2
            std::set_difference(composite_node_idx_neighbor1_extra.begin(), composite_node_idx_neighbor1_extra.end(), composite_inds_neighbor2_vector.begin(), composite_inds_neighbor2_vector.end(), std::inserter(diff_neighbor1, diff_neighbor1.begin()));
            TEUCHOS_TEST_FOR_EXCEPTION(diff_neighbor1.size() != 1, Exceptions::RuntimeError, "Process ID: " << comm_->getRank() << " - Region: " << region_idx << " - "
                                                                                                            << "Mislabeling of nodes obstructed the identification of the extra node: region node " << region_node_idx << " leads to diff_neighbor1.size()= " << diff_neighbor1.size() << " \n");
          }

          // Identification of external node from the side of neighbor2
          {
            // Compute the intersection between neighbourhood of neighbor2 node and neighbourhood of central node
            std::set_intersection(composite_inds_vector.begin(), composite_inds_vector.end(), composite_inds_neighbor2_vector.begin(), composite_inds_neighbor2_vector.end(), std::back_inserter(composite_node_idx_neighbor2_extra));

            for (typename std::vector<GlobalOrdinal>::iterator iter_node = composite_node_idx_neighbor2_extra.begin(); iter_node != composite_node_idx_neighbor2_extra.end(); ++iter_node) {
              checkerAllToRegion<GlobalOrdinal> unaryPredicateExtra(*iter_node + 1);
              typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator region_iterator_extra;
              region_iterator_extra = std::find_if<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerAllToRegion<GlobalOrdinal> >(regionToAll.begin(), regionToAll.end(), unaryPredicateExtra);

              // Invalidation of node indices for nodes belonging to the current region
              if (region_iterator_extra != regionToAll.end())
                *iter_node = -1;
            }

            // Removal of invalidated indices associated with nodes belonging to current region: (external nodes do not belong to this region)
            composite_node_idx_neighbor2_extra.erase(std::remove(composite_node_idx_neighbor2_extra.begin(), composite_node_idx_neighbor2_extra.end(), -1), composite_node_idx_neighbor2_extra.end());

            // External node from neighbor2 side does not belong to the neighborhood of neighbor1
            std::set_difference(composite_node_idx_neighbor2_extra.begin(), composite_node_idx_neighbor2_extra.end(), composite_inds_neighbor1_vector.begin(), composite_inds_neighbor1_vector.end(), std::inserter(diff_neighbor2, diff_neighbor2.begin()));
            TEUCHOS_TEST_FOR_EXCEPTION(diff_neighbor2.size() != 1, Exceptions::RuntimeError, "Process ID: " << comm_->getRank() << " - Region: " << region_idx << " - "
                                                                                                            << "Mislabeling of nodes obstructed the identification of the extra node: region node " << region_node_idx << " leads to diff_neighbor2.size()= " << diff_neighbor2.size() << " \n");
          }

          // Identification of external node from the side of central node
          {
            for (typename std::vector<GlobalOrdinal>::iterator iter_node = composite_node_idx_extra.begin(); iter_node != composite_node_idx_extra.end(); ++iter_node) {
              checkerAllToRegion<GlobalOrdinal> unaryPredicateExtra(*iter_node + 1);
              typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator region_iterator_extra;
              region_iterator_extra = std::find_if<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerAllToRegion<GlobalOrdinal> >(regionToAll.begin(), regionToAll.end(), unaryPredicateExtra);

              // Invalidation of node indices for nodes belonging to the current region
              if (region_iterator_extra != regionToAll.end())
                *iter_node = -1;
            }

            // Removal of invalidated indices associated with nodes belonging to current region: (external nodes do not belong to this region)
            composite_node_idx_extra.erase(std::remove(composite_node_idx_extra.begin(), composite_node_idx_extra.end(), -1), composite_node_idx_extra.end());
            std::vector<GlobalOrdinal> diff_center_temp;

            // At thie point composite_node_idx_extra contains indices of all the three external nodes: two of these must be removed since they are already tracked
            // External nodes from neighbors1's and neighbor2's side must be removed
            std::set_difference(composite_node_idx_extra.begin(), composite_node_idx_extra.end(), diff_neighbor1.begin(), diff_neighbor1.end(), std::inserter(diff_center_temp, diff_center_temp.begin()));
            std::set_difference(diff_center_temp.begin(), diff_center_temp.end(), diff_neighbor2.begin(), diff_neighbor2.end(), std::inserter(diff_center, diff_center.begin()));
            TEUCHOS_TEST_FOR_EXCEPTION(diff_center.size() != 1, Exceptions::RuntimeError, "Process ID: " << comm_->getRank() << " - Region: " << region_idx << " - "
                                                                                                         << "Mislabeling of nodes obstructed the identification of the extra node: region node " << region_node_idx << " leads to diff_center.size()= " << diff_center.size() << " \n");
          }

          // Computation of local indices for central node and its neighbors
          LocalOrdinal local_region_node_idx           = region_matrix->getRowMap()->getLocalElement(region_node_idx);
          LocalOrdinal local_region_node_idx_neighbor1 = region_matrix->getRowMap()->getLocalElement(region_node_idx_neighbor1);
          LocalOrdinal local_region_node_idx_neighbor2 = region_matrix->getRowMap()->getLocalElement(region_node_idx_neighbor2);

          // Computation of local indices for external nodes
          LocalOrdinal local_extra_central   = enlargedMatrix.getRowMap()->getLocalElement(diff_center[0]);
          LocalOrdinal local_extra_neighbor1 = enlargedMatrix.getRowMap()->getLocalElement(diff_neighbor1[0]);
          LocalOrdinal local_extra_neighbor2 = enlargedMatrix.getRowMap()->getLocalElement(diff_neighbor2[0]);

          ArrayView<const GlobalOrdinal> region_row;
          ArrayView<const GlobalOrdinal> region_col;
          ArrayView<const Scalar> region_val;

          // Extract Row view of the region matrix before collapsing
          if (region_matrix->isLocallyIndexed())
            region_matrix->getLocalRowView(local_region_node_idx, region_col, region_val);
          else
            region_matrix->getGlobalRowView(region_node_idx, region_col, region_val);

          // Extract Row of overlapped composite matrix to detect node with information to collapse
          ArrayView<const GlobalOrdinal> external_row;
          ArrayView<const GlobalOrdinal> external_col;
          ArrayView<const Scalar> external_val;
          enlargedMatrix.getLocalRowView(node_idx, external_col, external_val);

          // neighbor1 collapse
          {
            Scalar initial_value = 0;
            for (typename ArrayView<const GlobalOrdinal>::iterator iter_view = region_col.begin(); iter_view != region_col.end(); ++iter_view) {
              if (region_matrix->isLocallyIndexed())
                if (*iter_view == local_region_node_idx_neighbor1)
                  initial_value = region_val[iter_view - region_col.begin()];
              if (region_matrix->isGloballyIndexed())
                if (*iter_view == region_node_idx_neighbor1)
                  initial_value = region_val[iter_view - region_col.begin()];

              if (initial_value != 0)
                break;
            }

            Scalar external_value = 0;
            for (typename ArrayView<const GlobalOrdinal>::iterator iter_view = external_col.begin(); iter_view != external_col.end(); ++iter_view) {
              if (*iter_view == local_extra_neighbor1)
                external_value = external_val[iter_view - external_col.begin()];

              if (external_value != 0)
                break;
            }

            Scalar new_value = external_value;  // new matrix entry generated with the collapsing
            std::vector<GlobalOrdinal> new_entry_ind;
            std::vector<Scalar> new_entry_val;
            new_entry_ind.push_back(region_node_idx_neighbor1 - 1);
            new_entry_val.push_back(new_value);

            // If a nonzero value is already stored in the specified position, the new values is SUMMED to the already existing one
            // See description of insertGlobalValues(...)
            region_matrix->insertGlobalValues(region_node_idx - 1, new_entry_ind, new_entry_val);
          }
          // neighbor2 collapse
          {
            Scalar initial_value = 0;
            for (typename ArrayView<const GlobalOrdinal>::iterator iter_view = region_col.begin(); iter_view != region_col.end(); ++iter_view) {
              if (region_matrix->isLocallyIndexed())
                if (*iter_view == local_region_node_idx_neighbor2)
                  initial_value = region_val[iter_view - region_col.begin()];
              if (region_matrix->isGloballyIndexed())
                if (*iter_view == region_node_idx_neighbor2)
                  initial_value = region_val[iter_view - region_col.begin()];

              if (initial_value != 0)
                break;
            }

            Scalar external_value = 0;
            for (typename ArrayView<const GlobalOrdinal>::iterator iter_view = external_col.begin(); iter_view != external_col.end(); ++iter_view) {
              if (*iter_view == local_extra_neighbor2)
                external_value = external_val[iter_view - external_col.begin()];

              if (external_value != 0)
                break;
            }

            Scalar new_value = external_value;  // new matrix entry generated with the collapsing
            std::vector<GlobalOrdinal> new_entry_ind;
            std::vector<Scalar> new_entry_val;
            new_entry_ind.push_back(region_node_idx_neighbor2 - 1);
            new_entry_val.push_back(new_value);

            // If a nonzero value is already stored in the specified position, the new values is SUMMED to the already existing one
            // See description of insertGlobalValues(...)
            region_matrix->insertGlobalValues(region_node_idx - 1, new_entry_ind, new_entry_val);
          }
          // central node collapse
          {
            Scalar initial_value = 0;
            for (typename ArrayView<const GlobalOrdinal>::iterator iter_view = region_col.begin(); iter_view != region_col.end(); ++iter_view) {
              if (region_matrix->isLocallyIndexed())
                if (*iter_view == local_region_node_idx)
                  initial_value = region_val[iter_view - region_col.begin()];
              if (region_matrix->isGloballyIndexed())
                if (*iter_view == region_node_idx)
                  initial_value = region_val[iter_view - region_col.begin()];

              if (initial_value != 0)
                break;
            }

            Scalar external_value = 0;
            for (typename ArrayView<const GlobalOrdinal>::iterator iter_view = external_col.begin(); iter_view != external_col.end(); ++iter_view) {
              if (*iter_view == local_extra_central)
                external_value = external_val[iter_view - external_col.begin()];

              if (external_value != 0)
                break;
            }

            Scalar new_value = external_value;  // new matrix entry generated with the collapsing
            std::vector<GlobalOrdinal> new_entry_ind;
            std::vector<Scalar> new_entry_val;
            new_entry_ind.push_back(region_node_idx - 1);
            new_entry_val.push_back(new_value);

            // If a nonzero value is already stored in the specified position, the new values is SUMMED to the already existing one
            // See description of insertGlobalValues(...)
            region_matrix->insertGlobalValues(region_node_idx - 1, new_entry_ind, new_entry_val);
          }
        }
      }
    }
  }
  //@}

  //! @name Creation of Region Splitting
  //@{
  void RegionSplitting(GlobalOrdinal region_idx, RCP<Matrix> &region_matrix, Ifpack2::OverlappingRowMatrix<tpetra_row_matrix> &enlargedMatrix) {
    TEUCHOS_TEST_FOR_EXCEPTION(!region_matrix_initialized_[region_idx], Exceptions::RuntimeError, "The composite stiffness matrix must be chopped into surrogate region matrices before collapsing \n");
    TEUCHOS_TEST_FOR_EXCEPTION(regionHandler_->GetNumRegionNodes(region_idx) != regionMatrixData_[region_idx]->getGlobalNumRows(), Exceptions::RuntimeError, "Process ID: " << comm_->getRank() << " - Number of region nodes in region " << region_idx + 1 << " does not coincide with the value returned by regionMatrixData_[" << region_idx + 1 << "]->getGlobalNumRows() \n");

    Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > regionToAll = regionHandler_->GetRegionToAll(region_idx);

    // This portion of the code assumes that the number of region nodes is the same on each direction of the domain
    // Foir a 2D problem we have then nx = ny = sqrt( num_region_ndoes_ )
    GlobalOrdinal n;
    GlobalOrdinal nx;
    GlobalOrdinal ny;

    n  = regionHandler_->GetNumRegionNodes(region_idx);
    nx = std::sqrt(n);
    ny = nx;
    TEUCHOS_TEST_FOR_EXCEPTION(static_cast<double>(nx - std::floor(static_cast<double>(std::sqrt(static_cast<double>(n))))) != 0.0, Exceptions::RuntimeError, "The code assumes that the regions are 2D and that the number of region nodes is the same on each direction of the domain \n");

    // interfaceNodes contains nodes on an interface between any regions
    Array<std::tuple<int, Array<GlobalOrdinal> > > interfaceNodes = regionHandler_->GetInterfaceNodes();

    ArrayView<const GlobalOrdinal> MyRegionElements = region_matrix->getRowMap()->getLocalElementList();
    for (typename ArrayView<const GlobalOrdinal>::iterator iter = MyRegionElements.begin(); iter != MyRegionElements.end(); ++iter) {
      // Nodes are saved in data structures with 1 as base index
      GlobalOrdinal region_node_idx = *iter + 1;
      checkerRegionToAll<GlobalOrdinal> unaryPredicate(region_node_idx);
      typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator composite_iterator;
      composite_iterator = std::find_if<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerRegionToAll<GlobalOrdinal> >(regionToAll.begin(), regionToAll.end(), unaryPredicate);
      TEUCHOS_TEST_FOR_EXCEPTION(composite_iterator == regionToAll.end(), Exceptions::RuntimeError, "Process ID: " << comm_->getRank() << " - Region: " << region_idx << " - "
                                                                                                                   << " node with region index: " << region_node_idx << " is not in regionToAll[" << region_idx << "]"
                                                                                                                   << "\n");

      GlobalOrdinal composite_node_idx = std::get<1>(*composite_iterator);
      checkerInterfaceNodes<GlobalOrdinal> unaryPredicate2(composite_node_idx);
      typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator interface_iterator;
      interface_iterator = std::find_if<typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator, checkerInterfaceNodes<GlobalOrdinal> >(interfaceNodes.begin(), interfaceNodes.end(), unaryPredicate2);

      GlobalOrdinal region_node_idx_neighbor_w = 0;
      GlobalOrdinal region_node_idx_neighbor_e = 0;
      GlobalOrdinal region_node_idx_neighbor_n = 0;
      GlobalOrdinal region_node_idx_neighbor_s = 0;

      int count_neighbours = 0;

      if (interface_iterator != interfaceNodes.end() && region_node_idx > ny) {
        region_node_idx_neighbor_w = region_node_idx - ny;
        count_neighbours++;
      }

      if (interface_iterator != interfaceNodes.end() && region_node_idx <= (nx - 1) * ny) {
        region_node_idx_neighbor_e = region_node_idx + ny;
        count_neighbours++;
      }

      if (interface_iterator != interfaceNodes.end() && region_node_idx % ny != 1) {
        region_node_idx_neighbor_s = region_node_idx - 1;
        count_neighbours++;
      }

      if (interface_iterator != interfaceNodes.end() && region_node_idx % ny != 0) {
        region_node_idx_neighbor_n = region_node_idx + 1;
        count_neighbours++;
      }

      bool interface_line   = false;
      bool interface_corner = false;

      if (3 == count_neighbours)
        interface_line = true;
      else if (2 == count_neighbours)
        interface_corner = true;

      typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator composite_iterator_neighbor_e;
      typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator interface_iterator_neighbor_e;
      typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator composite_iterator_neighbor_w;
      typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator interface_iterator_neighbor_w;
      typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator composite_iterator_neighbor_s;
      typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator interface_iterator_neighbor_s;
      typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator composite_iterator_neighbor_n;
      typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator interface_iterator_neighbor_n;

      if (interface_line || interface_corner) {
        // Computation of composite index for East node
        if (region_node_idx_neighbor_e != 0) {
          checkerRegionToAll<GlobalOrdinal> unaryPredicateEast(region_node_idx_neighbor_e);
          composite_iterator_neighbor_e = std::find_if<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerRegionToAll<GlobalOrdinal> >(regionToAll.begin(), regionToAll.end(), unaryPredicateEast);

          // Check to see if neighbor_e node lies on a coarse line
          GlobalOrdinal composite_node_idx_neighbor_e = std::get<1>(*composite_iterator_neighbor_e);
          checkerInterfaceNodes<GlobalOrdinal> unaryPredicate2neighborEast(composite_node_idx_neighbor_e);
          interface_iterator_neighbor_e = std::find_if<typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator, checkerInterfaceNodes<GlobalOrdinal> >(interfaceNodes.begin(), interfaceNodes.end(), unaryPredicate2neighborEast);
        } else
          interface_iterator_neighbor_e = interfaceNodes.end();

        // Computation of composite index for West node
        if (region_node_idx_neighbor_w != 0) {
          checkerRegionToAll<GlobalOrdinal> unaryPredicateWest(region_node_idx_neighbor_w);
          composite_iterator_neighbor_w = std::find_if<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerRegionToAll<GlobalOrdinal> >(regionToAll.begin(), regionToAll.end(), unaryPredicateWest);

          // Check to see if neighbor_w node lies on a coarse line
          GlobalOrdinal composite_node_idx_neighbor_w = std::get<1>(*composite_iterator_neighbor_w);
          checkerInterfaceNodes<GlobalOrdinal> unaryPredicate2neighborWest(composite_node_idx_neighbor_w);
          interface_iterator_neighbor_w = std::find_if<typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator, checkerInterfaceNodes<GlobalOrdinal> >(interfaceNodes.begin(), interfaceNodes.end(), unaryPredicate2neighborWest);
        } else
          interface_iterator_neighbor_w = interfaceNodes.end();

        // Computation of composite index for South node
        if (region_node_idx_neighbor_s != 0) {
          checkerRegionToAll<GlobalOrdinal> unaryPredicateSouth(region_node_idx_neighbor_s);
          composite_iterator_neighbor_s = std::find_if<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerRegionToAll<GlobalOrdinal> >(regionToAll.begin(), regionToAll.end(), unaryPredicateSouth);

          // Check to see if neighbor_s node lies on a coarse line
          GlobalOrdinal composite_node_idx_neighbor_s = std::get<1>(*composite_iterator_neighbor_s);
          checkerInterfaceNodes<GlobalOrdinal> unaryPredicate2neighborSouth(composite_node_idx_neighbor_s);
          interface_iterator_neighbor_s = std::find_if<typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator, checkerInterfaceNodes<GlobalOrdinal> >(interfaceNodes.begin(), interfaceNodes.end(), unaryPredicate2neighborSouth);
        } else
          interface_iterator_neighbor_s = interfaceNodes.end();

        // Computation of composite index for North node
        if (region_node_idx_neighbor_n != 0) {
          checkerRegionToAll<GlobalOrdinal> unaryPredicateNorth(region_node_idx_neighbor_n);
          composite_iterator_neighbor_n = std::find_if<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerRegionToAll<GlobalOrdinal> >(regionToAll.begin(), regionToAll.end(), unaryPredicateNorth);

          // Check to see if neighbor_n node lies on a coarse line
          GlobalOrdinal composite_node_idx_neighbor_n = std::get<1>(*composite_iterator_neighbor_n);
          checkerInterfaceNodes<GlobalOrdinal> unaryPredicate2neighborNorth(composite_node_idx_neighbor_n);
          interface_iterator_neighbor_n = std::find_if<typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator, checkerInterfaceNodes<GlobalOrdinal> >(interfaceNodes.begin(), interfaceNodes.end(), unaryPredicate2neighborNorth);
        } else
          interface_iterator_neighbor_n = interfaceNodes.end();

        int count_neighbours_interface = 0;
        if (interface_iterator_neighbor_e != interfaceNodes.end())
          count_neighbours_interface++;
        if (interface_iterator_neighbor_w != interfaceNodes.end())
          count_neighbours_interface++;
        if (interface_iterator_neighbor_s != interfaceNodes.end())
          count_neighbours_interface++;
        if (interface_iterator_neighbor_n != interfaceNodes.end())
          count_neighbours_interface++;

        TEUCHOS_TEST_FOR_EXCEPTION(count_neighbours_interface > count_neighbours, Exceptions::RuntimeError, "Process ID: " << comm_->getRank() << " - Region: " << region_idx << " - "
                                                                                                                           << " node with region index: " << region_node_idx << " has inconsistent information on the number of neighbours: count_neighbours = " << count_neighbours << "but count_neighbours_interface =" << count_neighbours_interface << " \n");

        // First the splitting is applied on extradiagonal entries

        // Computation of local indices for central node and its neighbors
        // Node index base start from 1 in the structures used, but Trilinos maps start from 0, so
        // indices must be shifted by 1
        LocalOrdinal local_region_node_idx            = region_matrix->getRowMap()->getLocalElement(region_node_idx - 1);
        LocalOrdinal local_region_node_idx_neighbor_e = region_matrix->getRowMap()->getLocalElement(region_node_idx_neighbor_e - 1);
        LocalOrdinal local_region_node_idx_neighbor_w = region_matrix->getRowMap()->getLocalElement(region_node_idx_neighbor_w - 1);
        LocalOrdinal local_region_node_idx_neighbor_s = region_matrix->getRowMap()->getLocalElement(region_node_idx_neighbor_s - 1);
        LocalOrdinal local_region_node_idx_neighbor_n = region_matrix->getRowMap()->getLocalElement(region_node_idx_neighbor_n - 1);

        ArrayView<const GlobalOrdinal> region_col;
        ArrayView<const Scalar> region_val;

        // Extract Row view of the region matrix
        if (region_matrix->isLocallyIndexed())
          region_matrix->getLocalRowView(local_region_node_idx, region_col, region_val);
        else
          region_matrix->getGlobalRowView(*iter, region_col, region_val);

        std::vector<LocalOrdinal> region_col_vector = createVector(region_col);
        std::vector<LocalOrdinal> ind_vector(0);
        std::vector<Scalar> val_vector(0);

        // Extraction of the info about East neighbour to halve the associated entry in the matrix
        if (interface_iterator_neighbor_e != interfaceNodes.end()) {
          typename std::vector<GlobalOrdinal>::iterator iter_east_vector;
          GlobalOrdinal east_ind;
          if (region_matrix->isLocallyIndexed()) {
            iter_east_vector = std::find(region_col_vector.begin(), region_col_vector.end(), local_region_node_idx_neighbor_e);
            east_ind         = region_matrix->getRowMap()->getGlobalElement(*iter_east_vector);
          } else {
            iter_east_vector = std::find(region_col_vector.begin(), region_col_vector.end(), region_node_idx_neighbor_e - 1);
            east_ind         = *iter_east_vector;
          }
          Scalar east_val = -0.5 * region_val[iter_east_vector - region_col_vector.begin()];
          ind_vector.push_back(east_ind);
          val_vector.push_back(east_val);
        }

        // Extraction of the info about West neighbour to halve the associated entry in the matrix
        if (interface_iterator_neighbor_w != interfaceNodes.end()) {
          typename std::vector<GlobalOrdinal>::iterator iter_west_vector;
          GlobalOrdinal west_ind;
          if (region_matrix->isLocallyIndexed()) {
            iter_west_vector = std::find(region_col_vector.begin(), region_col_vector.end(), local_region_node_idx_neighbor_w);
            west_ind         = region_matrix->getRowMap()->getGlobalElement(*iter_west_vector);
          } else {
            iter_west_vector = std::find(region_col_vector.begin(), region_col_vector.end(), region_node_idx_neighbor_w - 1);
            west_ind         = *iter_west_vector;
          }
          Scalar west_val = -0.5 * region_val[iter_west_vector - region_col_vector.begin()];
          ind_vector.push_back(west_ind);
          val_vector.push_back(west_val);
        }

        // Extraction of the info about South neighbour to halve the associated entry in the matrix
        if (interface_iterator_neighbor_s != interfaceNodes.end()) {
          typename std::vector<GlobalOrdinal>::iterator iter_south_vector;
          GlobalOrdinal south_ind;
          if (region_matrix->isLocallyIndexed()) {
            iter_south_vector = std::find(region_col_vector.begin(), region_col_vector.end(), local_region_node_idx_neighbor_s);
            south_ind         = region_matrix->getRowMap()->getGlobalElement(*iter_south_vector);
          } else {
            iter_south_vector = std::find(region_col_vector.begin(), region_col_vector.end(), region_node_idx_neighbor_s - 1);
            south_ind         = *iter_south_vector;
          }
          Scalar south_val = -0.5 * region_val[iter_south_vector - region_col_vector.begin()];
          ind_vector.push_back(south_ind);
          val_vector.push_back(south_val);
        }

        // Extraction of the info about North neighbour to halve the associated entry in the matrix
        if (interface_iterator_neighbor_n != interfaceNodes.end()) {
          typename std::vector<GlobalOrdinal>::iterator iter_north_vector;
          GlobalOrdinal north_ind;
          if (region_matrix->isLocallyIndexed()) {
            iter_north_vector = std::find(region_col_vector.begin(), region_col_vector.end(), local_region_node_idx_neighbor_n);
            north_ind         = region_matrix->getRowMap()->getGlobalElement(*iter_north_vector);
          } else {
            iter_north_vector = std::find(region_col_vector.begin(), region_col_vector.end(), region_node_idx_neighbor_n - 1);
            north_ind         = *iter_north_vector;
          }
          Scalar north_val = -0.5 * region_val[iter_north_vector - region_col_vector.begin()];
          ind_vector.push_back(north_ind);
          val_vector.push_back(north_val);
        }

        // Extraction of the info about my Node ID to split the associated entry in the matrix
        // The ratio used for the splitting depends on the num ber of regions this current node
        // belongs to
        typename std::vector<LocalOrdinal>::iterator iter_center_vector;
        GlobalOrdinal center_ind;
        if (region_matrix->isLocallyIndexed()) {
          iter_center_vector = std::find(region_col_vector.begin(), region_col_vector.end(), local_region_node_idx);
          center_ind         = region_matrix->getRowMap()->getGlobalElement(*iter_center_vector);
        } else {
          iter_center_vector = std::find(region_col_vector.begin(), region_col_vector.end(), region_node_idx - 1);
          center_ind         = *iter_center_vector;
        }

        // Count of the nubmer of regions the current node belogns to
        GlobalOrdinal region_belonging = std::get<1>(*interface_iterator).size();
        TEUCHOS_TEST_FOR_EXCEPTION(region_belonging < 2, Exceptions::RuntimeError, "Process ID: " << comm_->getRank() << " - Region: " << region_idx << " - "
                                                                                                  << " node with composite index: " << std::get<0>(*interface_iterator) << " should lie on an interface between regions but the nubmer of regions it belongs to is only " << region_belonging << "\n");

        Scalar center_val;
        // If a node is on a corner between four itnerfaces, then each the entry A(node_idx,node_idx) must be split in four parts
        // otherwise the entry must be divided by two, similarly to what done for the neighbours

        center_val = -(1 - static_cast<Scalar>(1 / static_cast<Scalar>(region_belonging))) * region_val[iter_center_vector - region_col_vector.begin()];
        ind_vector.push_back(center_ind);
        val_vector.push_back(center_val);

        region_matrix->insertGlobalValues(region_node_idx - 1, ind_vector, val_vector);
      }
    }
  }
  //@}

  //! @name Creation of Region matrices
  //@{
  void CreateRegionMatrices(Array<Array<GlobalOrdinal> > region_maps) {
    TEUCHOS_TEST_FOR_EXCEPTION(num_total_regions_ != region_maps.size(), Exceptions::RuntimeError, "Number of regions does not match with the size of region_maps structure \n");

    regionMatrixData_.clear();

    for (int i = 0; i < num_total_regions_; ++i) {
      // Create Xpetra map for region stiffness matrix
      RCP<const Xpetra::Map<int, GlobalOrdinal, Node> > xpetraMap;
      xpetraMap        = Xpetra::MapFactory<int, GlobalOrdinal, Node>::Build(lib, regionHandler_->GetNumRegionNodes(i), region_maps[i], 0, comm_);
      int num_elements = xpetraMap->getGlobalNumElements();

      RCP<CrsMatrix> crs_matrix;
      if (Xpetra::UseEpetra == lib)
        crs_matrix = rcp(new EpetraCrsMatrix(xpetraMap, num_elements));
      else if (Xpetra::UseTpetra == lib)
        crs_matrix = rcp(new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xpetraMap, num_elements));
      else
        std::cerr << " The library to build matrices must be either Epetra or Tpetra \n";

      RCP<Matrix> matrixPointer = rcp(new CrsMatrixWrap(crs_matrix));
      regionMatrixData_.push_back(matrixPointer);
    }

    for (GlobalOrdinal i = 0; i < num_total_regions_; ++i)
      RegionMatrix(i);
  };
  //@}

  //@}

  //! Private variables
  //@{

  // The boolean finalDefaultView_ keep track of the status of the default view (= already updated or not)
  // See also MatrixSplitting::updateDefaultView()
  mutable bool finalDefaultView_;

  Array<bool> region_matrix_initialized_;

  RCP<const Teuchos::Comm<int> > comm_;

  //! Handling node-to-region mappings etc.
  RCP<RegionHandler<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionHandler_;

  //! The original (non-splitted) matrix
  RCP<Matrix> compositeMatrixData_;

  //! The matrix after splitting according to region layout
  Array<RCP<Matrix> > regionMatrixData_;

  GlobalOrdinal num_total_elements_ = 0;
  GlobalOrdinal num_total_regions_  = 0;

  //@}

};  // class MatrixSplitting

}  // namespace Xpetra

#endif  // XPETRA_MATRIXSPLITTING_HPP
