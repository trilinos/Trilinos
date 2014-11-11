// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_BLOCKEDCRSMATRIX_HPP
#define XPETRA_BLOCKEDCRSMATRIX_HPP

#include <Kokkos_DefaultNode.hpp>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Hashtable.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_CrsMatrixFactory.hpp"

#include "Xpetra_MapExtractor.hpp"

#include "Xpetra_Matrix.hpp"

/** \file Xpetra_BlockedCrsMatrix.hpp

  Declarations for the class Xpetra::BlockedCrsMatrix.
*/
namespace Xpetra {

  typedef std::string viewLabel_t;

  template <class Scalar = Matrix<>::scalar_type,
            class LocalOrdinal =
              typename Matrix<Scalar>::local_ordinal_type,
            class GlobalOrdinal =
              typename Matrix<Scalar, LocalOrdinal>::global_ordinal_type,
            class Node =
              typename Matrix<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
  class BlockedCrsMatrix :
    public Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  public:
    typedef Scalar scalar_type;
    typedef LocalOrdinal local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node node_type;

  private:
#undef XPETRA_BLOCKEDCRSMATRIX_SHORT
#include "Xpetra_UseShortNames.hpp"

  public:

    //! @name Constructor/Destructor Methods
    //@{

    //! Constructor
    /*!
     * \param rangeMaps range maps for all blocks
     * \param domainMaps domain maps for all blocks
     * \param npr extimated number of entries per row in each block(!)
     * \param pftype Xpetra profile type
     */
    BlockedCrsMatrix(Teuchos::RCP<const MapExtractor>& rangeMaps,
                     Teuchos::RCP<const MapExtractor>& domainMaps,
                     size_t npr, Xpetra::ProfileType pftype = Xpetra::DynamicProfile)
    : domainmaps_(domainMaps), rangemaps_(rangeMaps)
    {
      blocks_.reserve(Rows()*Cols());

      // add CrsMatrix objects in row,column order
      for (size_t r = 0; r < Rows(); ++r)
        for (size_t c = 0; c < Cols(); ++c)
          blocks_.push_back(CrsMatrixFactory::Build(getRangeMap(r), npr, pftype));

      // Default view
      CreateDefaultView();
    }

    //! Destructor
    virtual ~BlockedCrsMatrix() {}

    //@}


    //! @name Insertion/Removal Methods
    //@{

    //! Insert matrix entries, using global IDs.
    /** All index values must be in the global space.
      \pre \c globalRow exists as an ID in the global row map
      \pre <tt>isLocallyIndexed() == false</tt>
      \pre <tt>isStorageOptimized() == false</tt>

      \post <tt>isGloballyIndexed() == true</tt>

      \note If \c globalRow does not belong to the matrix on this node, then it
      will be communicated to the appropriate node when globalAssemble() is
      called (which will, at the latest, occur during the next call to
      fillComplete().) Otherwise, the entries will be inserted in the local
      matrix.

      \note If the matrix row already contains values at the indices
      corresponding to values in \c cols, then the new values will be summed
      with the old values; this may happen at insertion or during the next call
      to fillComplete().

      \note If <tt>hasColMap() == true</tt>, only (cols[i],vals[i]) where
      cols[i] belongs to the column map on this node will be inserted into the
      matrix.
      */
    void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal>& cols, const ArrayView<const Scalar>& vals) {
      throw Xpetra::Exceptions::RuntimeError("insertGlobalValues not supported by BlockedCrsMatrix");
    }

    //! Insert matrix entries, using local IDs.
    /** All index values must be in the local space.
      \pre \c localRow exists as an ID in the local row map
      \pre <tt>isGloballyIndexed() == false</tt>
      \pre <tt>isStorageOptimized() == false</tt>

      \post <tt>isLocallyIndexed() == true</tt>
      */
    void insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal>& cols, const ArrayView<const Scalar>& vals) {
      throw Xpetra::Exceptions::RuntimeError("insertLocalValues not supported by BlockedCrsMatrix");
    }

    void removeEmptyProcessesInPlace(const Teuchos::RCP<const Map>& newMap) {
      throw Xpetra::Exceptions::RuntimeError("removeEmptyProcesses not supported by BlockedCrsMatrix");
    }

    //! \brief Replace matrix entries, using global IDs.
    /** All index values must be in the global space.

      \pre \c globalRow is a global row belonging to the matrix on this node.

      \note If (globalRow,cols[i]) corresponds to an entry that is duplicated in
      this matrix row (likely because it was inserted more than once and
      fillComplete() has not been called in the interim), the behavior of this
      function is not defined. */
    void replaceGlobalValues(GlobalOrdinal globalRow,
                             const ArrayView<const GlobalOrdinal> &cols,
                             const ArrayView<const Scalar>        &vals) {
      throw Xpetra::Exceptions::RuntimeError("replaceGlobalValues not supported by BlockedCrsMatrix");
    }

    //! Replace matrix entries, using local IDs.
    /** All index values must be in the local space.
      Note that if a value is not already present for the specified location in
      the matrix, the input value will be ignored silently.
      */
    void replaceLocalValues(LocalOrdinal localRow,
                            const ArrayView<const LocalOrdinal> &cols,
                            const ArrayView<const Scalar>       &vals) {
      throw Xpetra::Exceptions::RuntimeError("replaceLocalValues not supported by BlockedCrsMatrix");
    }

    //! Set all matrix entries equal to scalar
    virtual void setAllToScalar(const Scalar& alpha) {
      throw Xpetra::Exceptions::RuntimeError("setAllToScalar not supported by BlockedCrsMatrix");
    }

    //! Scale the current values of a matrix, this = alpha*this.
    void scale(const Scalar& alpha) {
      throw Xpetra::Exceptions::RuntimeError("scale not supported by BlockedCrsMatrix");
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
    void resumeFill(const RCP< ParameterList >& params = null) {
      throw Xpetra::Exceptions::RuntimeError("resumeFill not supported for block matrices");
    }

    /*! \brief Signal that data entry is complete, specifying domain and range maps.

      Off-node indices are distributed (via globalAssemble()), indices are sorted,
      redundant indices are eliminated, and global indices are transformed to local
      indices.

      \pre  <tt>isFillActive() == true<tt>
      \pre <tt>isFillComplete()() == false<tt>

      \post <tt>isFillActive() == false<tt>
      \post <tt>isFillComplete() == true<tt>
      \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
      */
    void fillComplete(const RCP<const Map>& domainMap, const RCP<const Map>& rangeMap, const RCP<ParameterList>& params = null) {
      throw Xpetra::Exceptions::RuntimeError("fillComplete with arguments not supported for block matrices");
    }

    /*! \brief Signal that data entry is complete.

      Off-node entries are distributed (via globalAssemble()), repeated entries are
      summed, and global indices are transformed to local indices.

      \note This method calls fillComplete( getRowMap(), getRowMap(), os ).

      \pre  <tt>isFillActive() == true<tt>
      \pre <tt>isFillComplete()() == false<tt>

      \post <tt>isFillActive() == false<tt>
      \post <tt>isFillComplete() == true<tt>
      \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
      */
    void fillComplete(const RCP<ParameterList>& params = null) {
      for (size_t r = 0; r < Rows(); ++r)
        for (size_t c = 0; c < Cols(); ++c) {
          Teuchos::RCP<CrsMatrix> Ablock = getMatrix(r,c);

          if (Ablock != Teuchos::null && !Ablock->isFillComplete())
            Ablock->fillComplete(getDomainMap(c), getRangeMap(r), params);
        }

      // get full row map
      RCP<const Map> rangeMap = rangemaps_->getFullMap();
      fullrowmap_ = MapFactory::Build(rangeMap()->lib(), rangeMap()->getGlobalNumElements(), rangeMap()->getNodeElementList(), rangeMap()->getIndexBase(), rangeMap()->getComm());

      // build full col map
      fullcolmap_ = Teuchos::null; // delete old full column map

      std::vector<GO> colmapentries;
      for (size_t c = 0; c < Cols(); ++c) {
        // copy all local column lids of all block rows to colset
        std::set<GO> colset;
        for (size_t r = 0; r < Rows(); ++r) {
          Teuchos::RCP<CrsMatrix> Ablock = getMatrix(r,c);

          if (Ablock != Teuchos::null) {
            Teuchos::ArrayView<const GO> colElements = Ablock->getColMap()->getNodeElementList();
            Teuchos::RCP<const Map>      colmap      = Ablock->getColMap();
            copy(colElements.getRawPtr(), colElements.getRawPtr() + colElements.size(), inserter(colset, colset.begin()));
          }
        }

        // remove duplicates (entries which are in column maps of more than one block row)
        colmapentries.reserve(colmapentries.size() + colset.size());
        copy(colset.begin(), colset.end(), back_inserter(colmapentries));
        sort(colmapentries.begin(), colmapentries.end());
        typename std::vector<GO>::iterator gendLocation;
        gendLocation = std::unique(colmapentries.begin(), colmapentries.end());
        colmapentries.erase(gendLocation,colmapentries.end());
      }

      // sum up number of local elements
      size_t numGlobalElements = 0;
      Teuchos::reduceAll(*(rangeMap->getComm()), Teuchos::REDUCE_SUM, colmapentries.size(), Teuchos::outArg(numGlobalElements));

      // store global full column map
      const Teuchos::ArrayView<const GO> aView = Teuchos::ArrayView<const GO>(colmapentries);
      fullcolmap_ = MapFactory::Build(rangeMap->lib(), numGlobalElements, aView, 0, rangeMap->getComm());

    }

    //@}

    //! Returns the number of global rows.
    /** Undefined if isFillActive().
    */
    global_size_t getGlobalNumRows() const {
      global_size_t globalNumRows = 0;

      for (size_t row = 0; row < Rows(); row++)
        for (size_t col = 0; col < Cols(); col++)
          if (!getMatrix(row,col).is_null()) {
            globalNumRows += getMatrix(row,col)->getGlobalNumRows();
            break; // we need only one non-null matrix in a row
          }

      return globalNumRows;
    }

    //! \brief Returns the number of global columns in the matrix.
    /** Undefined if isFillActive().
    */
    global_size_t getGlobalNumCols() const {
      global_size_t globalNumCols = 0;

      for (size_t col = 0; col < Cols(); col++)
        for (size_t row = 0; row < Rows(); row++)
          if (!getMatrix(row,col).is_null()) {
            globalNumCols += getMatrix(row,col)->getGlobalNumCols();
            break; // we need only one non-null matrix in a col
          }

      return globalNumCols;
    }

    //! Returns the number of matrix rows owned on the calling node.
    size_t getNodeNumRows() const {
      global_size_t nodeNumRows = 0;

      for (size_t row = 0; row < Rows(); ++row)
        for (size_t col = 0; col < Cols(); col++)
          if (!getMatrix(row,col).is_null()) {
            nodeNumRows += getMatrix(row,col)->getNodeNumRows();
            break; // we need only one non-null matrix in a row
          }

      return nodeNumRows;
    }

    //! Returns the global number of entries in this matrix.
    global_size_t getGlobalNumEntries() const {
      global_size_t globalNumEntries = 0;

      for (size_t row = 0; row < Rows(); ++row)
        for (size_t col = 0; col < Cols(); ++col)
          if (!getMatrix(row,col).is_null())
            globalNumEntries += getMatrix(row,col)->getGlobalNumEntries();

      return globalNumEntries;
    }

    //! Returns the local number of entries in this matrix.
    size_t getNodeNumEntries() const {
      global_size_t nodeNumEntries = 0;

      for (size_t row = 0; row < Rows(); ++row)
        for (size_t col = 0; col < Cols(); ++col)
          if (!getMatrix(row,col).is_null())
            nodeNumEntries += getMatrix(row,col)->getNodeNumEntries();

      return nodeNumEntries;
    }

    //! Returns the current number of entries on this node in the specified local row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
    size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const {
      throw Xpetra::Exceptions::RuntimeError("getNumEntriesInLocalRow not supported by BlockedCrsMatrix");
    }

    //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons.
    /** Undefined if isFillActive().
    */
    global_size_t getGlobalNumDiags() const {
      throw Xpetra::Exceptions::RuntimeError("getGlobalNumDiags() not supported by BlockedCrsMatrix");
    }

    //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons.
    /** Undefined if isFillActive().
    */
    size_t getNodeNumDiags() const {
      throw Xpetra::Exceptions::RuntimeError("getNodeNumDiags() not supported by BlockedCrsMatrix");
    }

    //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
    /** Undefined if isFillActive().
    */
    size_t getGlobalMaxNumRowEntries() const {
      throw Xpetra::Exceptions::RuntimeError("getGlobalMaxNumRowEntries() not supported by BlockedCrsMatrix");
    }

    //! \brief Returns the maximum number of entries across all rows/columns on this node.
    /** Undefined if isFillActive().
    */
    size_t getNodeMaxNumRowEntries() const {
      throw Xpetra::Exceptions::RuntimeError("getNodeMaxNumRowEntries() not supported by BlockedCrsMatrix");
    }

    //! \brief If matrix indices of all matrix blocks are in the local range, this function returns true. Otherwise, this function returns false.
    /** if false, then this does not automatically mean that all blocks are globally indexed. The user has to make sure, that all matrix blocks
     * are indexed in the same way (locally or globally). Otherwise the block matrix is not valid...
     */
    bool isLocallyIndexed() const {
      for (size_t i = 0; i < blocks_.size(); ++i)
        if (blocks_[i] != Teuchos::null && !blocks_[i]->isLocallyIndexed())
          return false;
      return true;
    }

    //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false.
    /** if false, then this does not automatically mean that all blocks are locally indexed. The user has to make sure, that all matrix blocks
     * are indexed in the same way (locally or globally). Otherwise the block matrix is not valid...
     */
    bool isGloballyIndexed() const {
      for (size_t i = 0; i < blocks_.size(); i++)
        if (blocks_[i] != Teuchos::null && !blocks_[i]->isGloballyIndexed())
          return false;
      return true;
    }

    //! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
    bool isFillComplete() const {
      for (size_t i = 0; i < blocks_.size(); i++)
        if (blocks_[i] != Teuchos::null && !blocks_[i]->isFillComplete())
          return false;
      return true;
    }

    //! Extract a list of entries in a specified local row of the matrix. Put into storage allocated by calling routine.
    /*!
      \param LocalRow - (In) Local row number for which indices are desired.
      \param Indices - (Out) Local column indices corresponding to values.
      \param Values - (Out) Matrix values.
      \param NumIndices - (Out) Number of indices.

      Note: A std::runtime_error exception is thrown if either \c Indices or \c
      Values is not large enough to hold the data associated with row \c
      LocalRow. If \c LocalRow is not valid for this node, then \c Indices and
      \c Values are unchanged and \c NumIndices is returned as
      OrdinalTraits<size_t>::invalid().

      \pre <tt>isLocallyIndexed()==true</tt> or <tt>hasColMap() == true</tt>
      */
    virtual void getLocalRowCopy(LocalOrdinal LocalRow,
                                 const ArrayView<LocalOrdinal>& Indices,
                                 const ArrayView<Scalar>& Values,
                                 size_t &NumEntries) const {
      throw Xpetra::Exceptions::RuntimeError("getLocalRowCopy not supported by BlockedCrsMatrix" );
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
    void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal>& indices, ArrayView<const Scalar>& values) const {
      throw Xpetra::Exceptions::RuntimeError("getGlobalRowView not supported by BlockedCrsMatrix");
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
    void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal>& indices, ArrayView<const Scalar>& values) const {
      throw Xpetra::Exceptions::RuntimeError("getLocalRowView not supported by BlockedCrsMatrix");
    }

    //! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
    /*! Returns a distributed Vector object partitioned according to this
      matrix's row map, containing the
      the zero and non-zero diagonals owned by this node. */
    void getLocalDiagCopy(Vector& diag) const {
      throw Xpetra::Exceptions::RuntimeError("getLocalDiagCopy not supported by BlockedCrsMatrix" );
    }

    //! Get Frobenius norm of the matrix
    virtual typename ScalarTraits<Scalar>::magnitudeType getFrobeniusNorm() const {
      throw Xpetra::Exceptions::RuntimeError("getFrobeniusNorm() not supported by BlockedCrsMatrix, yet");
    }

    //@}

    //! @name Advanced Matrix-vector multiplication and solve methods
    //@{

    //! Multiplies this matrix by a MultiVector.
    /*! \c X is required to be post-imported, i.e., described by the column map
     * of the matrix. \c Y is required to be pre-exported, i.e., described by the
     * row map of the matrix.

     Both are required to have constant stride, and they are not permitted to
     ocupy overlapping space. No runtime checking will be performed in a
     non-debug build.

     This method is templated on the scalar type of MultiVector objects, allowing
     this method to be applied to MultiVector objects of arbitrary type. However,
     it is recommended that multiply() not be called directly; instead, use the
     CrsMatrixMultiplyOp, as it will handle the import/exprt operations required
     to apply a matrix with non-trivial communication needs.

     If \c beta is equal to zero, the operation will enjoy overwrite semantics
     (\c Y will be overwritten with the result of the multiplication). Otherwise,
     the result of the multiplication will be accumulated into \c Y.
     */

    //@}

    //! @name Methods implementing Matrix
    //@{

    //! \brief Computes the sparse matrix-multivector multiplication.
    /*! Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exceptions:
      - if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
      */
    virtual void apply(const MultiVector& X, MultiVector& Y,
                       Teuchos::ETransp mode = Teuchos::NO_TRANS,
                       Scalar alpha = ScalarTraits<Scalar>::one(),
                       Scalar beta  = ScalarTraits<Scalar>::zero()) const
    {
      using Teuchos::RCP;

      TEUCHOS_TEST_FOR_EXCEPTION(mode != Teuchos::NO_TRANS && mode != Teuchos::TRANS, Xpetra::Exceptions::RuntimeError,
                                 "apply() only supports the following modes: NO_TRANS and TRANS." );

      RCP<const MultiVector> refX = rcpFromRef(X);
      RCP<MultiVector>       tmpY = MultiVectorFactory::Build(Y.getMap(), Y.getNumVectors());

      SC one = ScalarTraits<SC>::one();

      if (mode == Teuchos::NO_TRANS) {
        for (size_t row = 0; row < Rows(); row++) {
          RCP<MultiVector>    Yblock = rangemaps_->getVector(row, Y.getNumVectors());
          RCP<MultiVector> tmpYblock = rangemaps_->getVector(row, Y.getNumVectors());

          for (size_t col = 0; col < Cols(); col++) {
            RCP<const MultiVector> Xblock = domainmaps_->ExtractVector(refX, col);
            RCP<CrsMatrix>         Ablock = getMatrix(row, col);

            if (Ablock.is_null())
              continue;

            Ablock->apply(*Xblock, *tmpYblock);

            Yblock->update(one, *tmpYblock, one);
          }
          rangemaps_->InsertVector(Yblock, row, tmpY);
        }

      } else if (mode == Teuchos::TRANS) {
        // TODO: test me!
        for (size_t col = 0; col < Cols(); col++) {
          RCP<MultiVector>    Yblock = domainmaps_->getVector(col, Y.getNumVectors());
          RCP<MultiVector> tmpYblock = domainmaps_->getVector(col, Y.getNumVectors());

          for (size_t row = 0; row<Rows(); row++) {
            RCP<const MultiVector> Xblock = rangemaps_->ExtractVector(refX, row);
            RCP<CrsMatrix>         Ablock = getMatrix(row, col);

            if (Ablock.is_null())
              continue;

            Ablock->apply(*Xblock, *tmpYblock, Teuchos::TRANS);

            Yblock->update(one, *tmpYblock, one);
          }
          domainmaps_->InsertVector(Yblock, col, tmpY);
        }
      }

      Y.update(alpha, *tmpY, beta);
    }

    //! \brief Returns the Map associated with the full domain of this operator.
    //! This will be <tt>null</tt> until fillComplete() is called.
    RCP<const Map > getDomainMap() const            { return domainmaps_->getFullMap(); }

    //! \brief Returns the Map associated with the i'th block domain of this operator.
    //! This will be <tt>null</tt> until fillComplete() is called.
    RCP<const Map > getDomainMap(size_t i) const    { return domainmaps_->getMap(i); }

    //! Returns the Map associated with the full range of this operator.
    //! This will be <tt>null</tt> until fillComplete() is called.
    RCP<const Map > getRangeMap() const             { return rangemaps_->getFullMap(); }

    //! Returns the Map associated with the i'th block range of this operator.
    //! This will be <tt>null</tt> until fillComplete() is called.
    RCP<const Map > getRangeMap(size_t i) const     { return rangemaps_->getMap(i); }

    //! Returns map extractor class for range map
    RCP<const MapExtractor> getRangeMapExtractor()  { return rangemaps_; }

    //! Returns map extractor for domain map
    RCP<const MapExtractor> getDomainMapExtractor() { return domainmaps_; }

    //@}

    //! Implements DistObject interface
    //{@

    //! Access function for the Tpetra::Map this DistObject was constructed with.
    const Teuchos::RCP< const Map > getMap() const {
      throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::getMap(): operation not supported.");
    }

    //! Import.
    void doImport(const Matrix &source, const Import& importer, CombineMode CM) {
      throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::doImport(): operation not supported.");
    }

    //! Export.
    void doExport(const Matrix& dest, const Import& importer, CombineMode CM) {
      throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::doExport(): operation not supported.");
    }

    //! Import (using an Exporter).
    void doImport(const Matrix& source, const Export& exporter, CombineMode CM) {
      throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::doImport(): operation not supported.");
    }

    //! Export (using an Importer).
    void doExport(const Matrix& dest, const Export& exporter, CombineMode CM) {
      throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::doExport(): operation not supported.");
    }

    // @}

    //! @name Overridden from Teuchos::Describable
    //@{

    /** \brief Return a simple one-line description of this object. */
    std::string description() const                       { return "Xpetra_BlockedCrsMatrix.description()"; }

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const {
      out << "Xpetra::BlockedCrsMatrix: " << Rows() << " x " << Cols() << std::endl;

      if (isFillComplete()) {
        out << "BlockMatrix is fillComplete" << std::endl;

        out << "fullRowMap" << std::endl;
        fullrowmap_->describe(out,verbLevel);

        out << "fullColMap" << std::endl;
        fullcolmap_->describe(out,verbLevel);

      } else {
        out << "BlockMatrix is NOT fillComplete" << std::endl;
      }

      for (size_t r = 0; r < Rows(); ++r)
        for (size_t c = 0; c < Cols(); ++c) {
          out << "Block(" << r << "," << c << ")" << std::endl;
          getMatrix(r,c)->describe(out,verbLevel);
        }
    }

    //! Returns the CrsGraph associated with this matrix.
    RCP<const CrsGraph> getCrsGraph() const {
      throw Xpetra::Exceptions::RuntimeError("getCrsGraph() not supported by BlockedCrsMatrix");
    }

    //@}

    //! @name Block matrix access
    //@{

    /// number of row blocks
    virtual size_t Rows() const                                       { return rangemaps_->NumMaps(); }

    /// number of column blocks
    virtual size_t Cols() const                                       { return domainmaps_->NumMaps(); }

    /// return block (r,c)
    Teuchos::RCP<CrsMatrix> getMatrix(size_t r, size_t c) const       { return blocks_[r*Cols()+c]; }

    /// set matrix block
    void setMatrix(size_t r, size_t c, Teuchos::RCP<CrsMatrix> mat) {
      // TODO: if filled -> return error

      TEUCHOS_TEST_FOR_EXCEPTION(r > Rows(), std::out_of_range, "Error, r = " << Rows() << " is too big");
      TEUCHOS_TEST_FOR_EXCEPTION(c > Cols(), std::out_of_range, "Error, c = " << Cols() << " is too big");

      // set matrix
      blocks_[r*Cols() + c] = mat;
    }

    /// merge BlockedCrsMatrix blocks in a CrsMatrix
    // NOTE: This is a rather expensive operation, since all blocks are copied
    // into a new big CrsMatrix
    Teuchos::RCP<CrsMatrix> Merge() const {
      using Teuchos::RCP;
      using Teuchos::rcp_dynamic_cast;
      Scalar one = ScalarTraits<SC>::one();

      RCP<CrsMatrix> sparse = CrsMatrixFactory::Build(fullrowmap_, 33);

      for (size_t i = 0; i < blocks_.size(); i++)
        if (blocks_[i] != Teuchos::null)
          this->Add(*blocks_[i], one, *sparse, one);

      sparse->fillComplete(getDomainMap(), getRangeMap());

      return sparse;
    }
    //@}


  private:

    /** \name helper functions */
    //@{

    /// Add a Xpetra::CrsMatrix to another: B = B*scalarB + A*scalarA
    /**
     * Note, that this routine works only correctly if A only has entries which are empty (zero) in B.
     * We use the insertGlobalValues routine for inserting the new values from A in B. The sumIntoGlobalValues
     * routine is not implemented in Xpetra (and would not extend the graph of B for new entries).
     * Here we need something to catch the exceptions of a future implementation of sumIntoGlobalValues that
     * then adds the remaining new entries with insertGlobal Values.
     *
     * This routine is private and used only by Merge. Since the blocks in BlockedCrsMatrix are seperated,
     * this routine works for merging a BlockedCrsMatrix.
     */
    void Add(const CrsMatrix& A, const Scalar scalarA, CrsMatrix& B, const Scalar scalarB) const {
      TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Xpetra::Exceptions::RuntimeError,
                                 "Matrix A is not completed");
      using Teuchos::Array;
      using Teuchos::ArrayView;

      B.scale(scalarB);

      Scalar one  = ScalarTraits<SC>::one();
      Scalar zero = ScalarTraits<SC>::zero();

      if (scalarA == zero)
        return;

      size_t maxNumEntries = A.getNodeMaxNumRowEntries();

      size_t    numEntries;
      Array<GO> inds(maxNumEntries);
      Array<SC> vals(maxNumEntries);

      RCP<const Map> rowMap = A.getRowMap();
      RCP<const Map> colMap = A.getColMap();

      ArrayView<const GO> rowGIDs = A.getRowMap()->getNodeElementList();
      for (size_t i = 0; i < A.getNodeNumRows(); i++) {
        GO row = rowGIDs[i];
        A.getGlobalRowCopy(row, inds(), vals(), numEntries);

        if (scalarA != one)
          for (size_t j = 0; j < numEntries; ++j)
            vals[j] *= scalarA;

        B.insertGlobalValues(row, inds(0, numEntries), vals(0, numEntries)); // insert should be ok, since blocks in BlockedCrsOpeartor do not overlap!
      }
    }

    //@}

    // Default view is created after fillComplete()
    // Because ColMap might not be available before fillComplete().
    void CreateDefaultView() {
      // Create default view
      this->defaultViewLabel_ = "point";
      this->CreateView(this->GetDefaultViewLabel(), getRangeMap(), getDomainMap());

      // Set current view
      this->currentViewLabel_ = this->GetDefaultViewLabel();
    }

  private:
    Teuchos::RCP<const MapExtractor>      domainmaps_;        // full        domain map together with all partial domain maps
    Teuchos::RCP<const MapExtractor>      rangemaps_;         // full         range map together with all partial domain maps
    Teuchos::RCP<Map>                     fullrowmap_;        // full matrix    row map
    Teuchos::RCP<Map>                     fullcolmap_;        // full matrix column map

    std::vector<Teuchos::RCP<CrsMatrix> > blocks_;            // row major matrix block storage
};

} //namespace Xpetra

#define XPETRA_BLOCKEDCRSMATRIX_SHORT
#endif /* XPETRA_BLOCKEDCRSMATRIX_HPP */
