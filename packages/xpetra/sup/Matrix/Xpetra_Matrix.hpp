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

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_MATRIX_HPP
#define XPETRA_MATRIX_HPP

#include <Kokkos_DefaultNode.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_MultiVector.hpp"
#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_CrsMatrixFactory.hpp"
#include "Xpetra_MatrixView.hpp"
#include "Xpetra_Operator.hpp"
#include "Xpetra_StridedMap.hpp"
#include "Xpetra_StridedMapFactory.hpp"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Hashtable.hpp>

/** \file Xpetra_Matrix.hpp

Declarations for the class Xpetra::Matrix.
*/
namespace Xpetra {

    /*!
     @class Xpetra::Matrix class.
     @brief Xpetra-specific matrix class.

     This class is specific to Xpetra and has no analogue in Epetra or Tpetra.  The main motivation for this class is to be able to access matrix data in a manner different than how it is stored.
     For example, it might be more convenient to treat ("view") a matrix stored in compressed row storage as if it were a block matrix.  The Xpetra::Matrix class is intended to manage these "views".

     <B>How to create a Matrix from an existing CrsMatrix</B>

     @code
     RCP<Xpetra::CrsMatrix> crsA;
     RCP<Xpetra::Matrix>    A  = rcp(new CrsMatrixWrap(crsA));
     @endcode

    */

  typedef std::string viewLabel_t;

  template <class Scalar        = Operator<>::scalar_type,
            class LocalOrdinal  = Operator<>::local_ordinal_type,
            class GlobalOrdinal = typename Operator<LocalOrdinal>::global_ordinal_type,
            class Node          = typename Operator<LocalOrdinal, GlobalOrdinal>::node_type>
  class Matrix : public Xpetra::Operator< Scalar, LocalOrdinal, GlobalOrdinal, Node > {
    typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
    typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrix;
    typedef Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> CrsGraph;
    typedef Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrixFactory;
    typedef Xpetra::MatrixView<Scalar, LocalOrdinal, GlobalOrdinal, Node> MatrixView;

  public:
    typedef Scalar          scalar_type;
    typedef LocalOrdinal    local_ordinal_type;
    typedef GlobalOrdinal   global_ordinal_type;
    typedef Node            node_type;

#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
#ifdef HAVE_XPETRA_TPETRA
    typedef typename CrsMatrix::local_matrix_type local_matrix_type;
#endif
#endif

    //! @name Constructor/Destructor Methods
    //@{

    Matrix() { }

    //! Destructor
    virtual ~Matrix() { }

    //@}

    //! @name View management methods
    //@{
    void CreateView(viewLabel_t viewLabel, const RCP<const Map> & rowMap, const RCP<const Map> & colMap) {
      TEUCHOS_TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == true, Xpetra::Exceptions::RuntimeError, "Xpetra::Matrix.CreateView(): a view labeled '" + viewLabel + "' already exist.");
      RCP<MatrixView> view = rcp(new MatrixView(rowMap, colMap));
      operatorViewTable_.put(viewLabel, view);
    }

    // JG TODO: why this is a member function??
    void CreateView(const viewLabel_t viewLabel, const RCP<const Matrix>& A, bool transposeA = false, const RCP<const Matrix>& B = Teuchos::null, bool transposeB = false) {
      RCP<const Map> domainMap = Teuchos::null;
      RCP<const Map> rangeMap  = Teuchos::null;

      const size_t        blkSize = 1;
      std::vector<size_t> stridingInfo(1, blkSize);
      LocalOrdinal        stridedBlockId = -1;


      if (A->IsView(viewLabel)) {
        rangeMap  = transposeA ? A->getColMap(viewLabel) : A->getRowMap(viewLabel);
        domainMap = transposeA ? A->getRowMap(viewLabel) : A->getColMap(viewLabel); // will be overwritten if B != Teuchos::null

      } else {
        rangeMap  = transposeA ? A->getDomainMap()       : A->getRangeMap();
        domainMap = transposeA ? A->getRangeMap()        : A->getDomainMap();

        if (viewLabel == "stridedMaps") {
          rangeMap  = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(rangeMap,  stridingInfo, stridedBlockId);
          domainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(domainMap, stridingInfo, stridedBlockId);
        }
      }

      if (B != Teuchos::null ) {
        // B has strided Maps

        if (B->IsView(viewLabel)) {
          domainMap = transposeB ? B->getRowMap(viewLabel) : B->getColMap(viewLabel);

        } else {
          domainMap = transposeB ? B->getRangeMap()        : B->getDomainMap();

          if (viewLabel == "stridedMaps")
            domainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(domainMap, stridingInfo, stridedBlockId);
        }
      }


      if (IsView(viewLabel))
        RemoveView(viewLabel);

      CreateView(viewLabel, rangeMap, domainMap);
    }

    //! Print all of the views associated with the Matrix.
    void PrintViews(Teuchos::FancyOStream &out) const {
      int last = out.getOutputToRootOnly();
      Teuchos::OSTab tab(out);
      out.setOutputToRootOnly(0);
      Teuchos::Array<viewLabel_t> viewLabels;
      Teuchos::Array<RCP<MatrixView> > viewList;
      operatorViewTable_.arrayify(viewLabels,viewList);
      out << "views associated with this operator" << std::endl;
      for (int i=0; i<viewLabels.size(); ++i)
        out << viewLabels[i] << std::endl;
      out.setOutputToRootOnly(last);
    }


    void RemoveView(const viewLabel_t viewLabel) {
      TEUCHOS_TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == false, Xpetra::Exceptions::RuntimeError, "Xpetra::Matrix.RemoveView(): view '" + viewLabel + "' does not exist.");
      TEUCHOS_TEST_FOR_EXCEPTION(viewLabel == GetDefaultViewLabel(), Xpetra::Exceptions::RuntimeError, "Xpetra::Matrix.RemoveView(): view '" + viewLabel + "' is the default view and cannot be removed.");
      operatorViewTable_.remove(viewLabel);
    }

    const viewLabel_t SwitchToView(const viewLabel_t viewLabel) {
      TEUCHOS_TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == false, Xpetra::Exceptions::RuntimeError, "Xpetra::Matrix.SwitchToView(): view '" + viewLabel + "' does not exist.");
      viewLabel_t oldViewLabel = GetCurrentViewLabel();
      currentViewLabel_ = viewLabel;
      return oldViewLabel;
    }

    bool IsView(const viewLabel_t viewLabel) const {
      return operatorViewTable_.containsKey(viewLabel);
    }

    const viewLabel_t SwitchToDefaultView() { return SwitchToView(GetDefaultViewLabel()); }

    const viewLabel_t & GetDefaultViewLabel() const { return defaultViewLabel_; }

    const viewLabel_t & GetCurrentViewLabel() const { return currentViewLabel_; }

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
    virtual void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &cols, const ArrayView<const Scalar> &vals) = 0;

    //! Insert matrix entries, using local IDs.
    /** All index values must be in the local space.
        \pre \c localRow exists as an ID in the local row map
        \pre <tt>isGloballyIndexed() == false</tt>
        \pre <tt>isStorageOptimized() == false</tt>

        \post <tt>isLocallyIndexed() == true</tt>
    */
    virtual void insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &cols, const ArrayView<const Scalar> &vals) = 0;

    //! \brief Replace matrix entries, using global IDs.
    /** All index values must be in the global space.

        \pre \c globalRow is a global row belonging to the matrix on this node.

        \note If (globalRow,cols[i]) corresponds to an entry that is duplicated in this matrix row (likely because it was inserted more than once and fillComplete() has not been called in the interim), the behavior of this function is not defined. */
    virtual void replaceGlobalValues(GlobalOrdinal globalRow,
                                     const ArrayView<const GlobalOrdinal> &cols,
                                     const ArrayView<const Scalar>        &vals) = 0;

    //! Replace matrix entries, using local IDs.
    /** All index values must be in the local space.
        Note that if a value is not already present for the specified location in the matrix, the input value will be ignored silently.
     */
    virtual void replaceLocalValues(LocalOrdinal localRow,
                                    const ArrayView<const LocalOrdinal> &cols,
                                    const ArrayView<const Scalar>       &vals) = 0;

    //! Set all matrix entries equal to scalar
    virtual void setAllToScalar(const Scalar &alpha)= 0;

    //! Scale the current values of a matrix, this = alpha*this.
    virtual void scale(const Scalar &alpha)= 0;

    //@}

    //! @name Transformational Methods
    //@{

    /*! Resume fill operations.
      After calling fillComplete(), resumeFill() must be called before initiating any changes to the matrix.

      resumeFill() may be called repeatedly.

      \post  <tt>isFillActive() == true<tt>
      \post  <tt>isFillComplete() == false<tt>
    */
    virtual void resumeFill(const RCP< ParameterList > &params=null) = 0;

    /*! \brief Signal that data entry is complete, specifying domain and range maps.

    Off-node indices are distributed (via globalAssemble()), indices are sorted, redundant indices are eliminated, and global indices are transformed to local indices.

    \pre  <tt>isFillActive() == true<tt>
    \pre <tt>isFillComplete()() == false<tt>

    \post <tt>isFillActive() == false<tt>
    \post <tt>isFillComplete() == true<tt>
    \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
    */
    virtual void fillComplete(const RCP<const Map> &domainMap, const RCP<const Map> &rangeMap, const RCP<ParameterList> &params = null) =0;

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
    virtual void fillComplete(const RCP<ParameterList> &params=null) =0;

    //@}

    //! @name Methods implementing RowMatrix
    //@{

    //! Returns the Map that describes the row distribution in this matrix.
    virtual const RCP<const Map> & getRowMap() const { return getRowMap(GetCurrentViewLabel()); }

    //! Returns the Map that describes the row distribution in this matrix.
    virtual const RCP<const Map> & getRowMap(viewLabel_t viewLabel) const {
      TEUCHOS_TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == false, Xpetra::Exceptions::RuntimeError, "Xpetra::Matrix.GetRowMap(): view '" + viewLabel + "' does not exist.");
      return operatorViewTable_.get(viewLabel)->GetRowMap();
    }

    //! \brief Returns the Map that describes the column distribution in this matrix.
    //! This might be <tt>null</tt> until fillComplete() is called.
    virtual const RCP<const Map> & getColMap() const { return getColMap(GetCurrentViewLabel()); }

    //! \brief Returns the Map that describes the column distribution in this matrix.
    virtual const RCP<const Map> & getColMap(viewLabel_t viewLabel) const {
      TEUCHOS_TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == false, Xpetra::Exceptions::RuntimeError, "Xpetra::Matrix.GetColMap(): view '" + viewLabel + "' does not exist.");
      return operatorViewTable_.get(viewLabel)->GetColMap();
    }

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

    //! Returns the global number of entries in this matrix.
    virtual global_size_t getGlobalNumEntries() const =0;

    //! Returns the local number of entries in this matrix.
    virtual size_t getNodeNumEntries() const =0;

    //! Returns the current number of entries on this node in the specified local row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
    virtual size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const =0;

    //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
    /** Undefined if isFillActive().
     */
    virtual size_t getGlobalMaxNumRowEntries() const =0;

    //! \brief Returns the maximum number of entries across all rows/columns on this node.
    /** Undefined if isFillActive().
     */
    virtual size_t getNodeMaxNumRowEntries() const =0;

    //! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
    virtual bool isLocallyIndexed() const =0;

    //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
    virtual bool isGloballyIndexed() const =0;

    //! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
    virtual bool isFillComplete() const =0;

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
      \param Indices  - (Out) Local column indices corresponding to values.
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

    //! Get Frobenius norm of the matrix
    virtual typename ScalarTraits<Scalar>::magnitudeType getFrobeniusNorm() const = 0;

    //! Left scale matrix using the given vector entries
    virtual void leftScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) = 0;

    //! Right scale matrix using the given vector entries
    virtual void rightScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) = 0;


    //! Returns true if globalConstants have been computed; false otherwise
    virtual bool haveGlobalConstants() const = 0;


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
//     virtual void multiply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans, Scalar alpha, Scalar beta) const=0;

    //@}

    //! Implements DistObject interface
    //{@

    //! Access function for the Tpetra::Map this DistObject was constructed with.
    virtual const Teuchos::RCP< const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > getMap() const = 0;

    // TODO: first argument of doImport/doExport should be a Xpetra::DistObject

    //! Import.
    virtual void doImport(const Matrix &source,
                          const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM) = 0;

    //! Export.
    virtual void doExport(const Matrix &dest,
                          const Import< LocalOrdinal, GlobalOrdinal, Node >& importer, CombineMode CM) = 0;

    //! Import (using an Exporter).
    virtual void doImport(const Matrix &source,
                          const Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM) = 0;

    //! Export (using an Importer).
    virtual void doExport(const Matrix &dest,
                          const Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM) = 0;

    // @}

    //! @name Overridden from Teuchos::Describable
    //@{

    // TODO: describe of views can be done here

    //   /** \brief Return a simple one-line description of this object. */
    //   virtual std::string description() const =0;

    //   /** \brief Print the object with some verbosity level to an FancyOStream object. */
    //   virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const =0;

    //@}


    //! @name Overridden from Teuchos::Describable
    //@{

    /** \brief Return a simple one-line description of this object. */
    virtual std::string description() const =0;

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const =0;
    //@}

    // JG: Added:

    //! Supports the getCrsGraph() call
    virtual bool hasCrsGraph() const =0;
    
    //! Returns the CrsGraph associated with this matrix.
    virtual RCP<const CrsGraph> getCrsGraph() const =0;

    // ----------------------------------------------------------------------------------
    // "TEMPORARY" VIEW MECHANISM
    /**
     * Set fixed block size of operator (e.g., 3 for 3 DOFs per node).
     *
     * @param blksize: block size denoting how many DOFs per node are used (LocalOrdinal)
     * @param offset:  global offset allows to define operators with global indices starting from a given value "offset" instead of 0. (GlobalOrdinal, default = 0)
     * */
    void SetFixedBlockSize(LocalOrdinal blksize, GlobalOrdinal offset=0) {

      TEUCHOS_TEST_FOR_EXCEPTION(isFillComplete() == false, Exceptions::RuntimeError, "Xpetra::Matrix::SetFixedBlockSize(): operator is not filled and completed."); // TODO: do we need this? we just wanna "copy" the domain and range map

      std::vector<size_t> stridingInfo;
      stridingInfo.push_back(Teuchos::as<size_t>(blksize));
      LocalOrdinal stridedBlockId = -1;

      RCP<const Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node> > stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                                    this->getRangeMap(),
                                                    stridingInfo,
                                                    stridedBlockId,
                                                    offset
                                                    );
      RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                              this->getDomainMap(),
                                              stridingInfo,
                                              stridedBlockId,
                                              offset
                                              );

      if(IsView("stridedMaps") == true) RemoveView("stridedMaps");
      CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);
    }

    //==========================================================================

    LocalOrdinal GetFixedBlockSize() const {
      if(IsView("stridedMaps")==true) {
        Teuchos::RCP<const StridedMap<LocalOrdinal, GlobalOrdinal, Node> > rangeMap = Teuchos::rcp_dynamic_cast<const StridedMap<LocalOrdinal, GlobalOrdinal, Node> >(getRowMap("stridedMaps"));
        Teuchos::RCP<const StridedMap<LocalOrdinal, GlobalOrdinal, Node> > domainMap = Teuchos::rcp_dynamic_cast<const StridedMap<LocalOrdinal, GlobalOrdinal, Node> >(getColMap("stridedMaps"));
        TEUCHOS_TEST_FOR_EXCEPTION(rangeMap  == Teuchos::null, Exceptions::BadCast, "Xpetra::Matrix::GetFixedBlockSize(): rangeMap is not of type StridedMap");
        TEUCHOS_TEST_FOR_EXCEPTION(domainMap == Teuchos::null, Exceptions::BadCast, "Xpetra::Matrix::GetFixedBlockSize(): domainMap is not of type StridedMap");
        TEUCHOS_TEST_FOR_EXCEPTION(domainMap->getFixedBlockSize() != rangeMap->getFixedBlockSize(), Exceptions::RuntimeError, "Xpetra::Matrix::GetFixedBlockSize(): block size of rangeMap and domainMap are different.");
        return Teuchos::as<LocalOrdinal>(domainMap->getFixedBlockSize()); // TODO: why LocalOrdinal?
      } else
        //TEUCHOS_TEST_FOR_EXCEPTION(false, Exceptions::RuntimeError, "Xpetra::Matrix::GetFixedBlockSize(): no strided maps available."); // TODO remove this
        return 1;
    }; //TODO: why LocalOrdinal?

    // ----------------------------------------------------------------------------------

    virtual void SetMaxEigenvalueEstimate(Scalar const &sigma) {
      operatorViewTable_.get(GetCurrentViewLabel())->SetMaxEigenvalueEstimate(sigma);
    }

    // ----------------------------------------------------------------------------------

    virtual Scalar GetMaxEigenvalueEstimate() const {
      return operatorViewTable_.get(GetCurrentViewLabel())->GetMaxEigenvalueEstimate();
    }

    // ----------------------------------------------------------------------------------
#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
#ifdef HAVE_XPETRA_TPETRA
    /// \brief Access the underlying local Kokkos::CrsMatrix object
    virtual local_matrix_type getLocalMatrix () const = 0;
#else
#ifdef __GNUC__
#warning "Xpetra Kokkos interface for CrsMatrix is enabled (HAVE_XPETRA_KOKKOS_REFACTOR) but Tpetra is disabled. The Kokkos interface needs Tpetra to be enabled, too."
#endif
#endif
#endif
    // ----------------------------------------------------------------------------------

    protected:
      Teuchos::Hashtable<viewLabel_t, RCP<MatrixView> > operatorViewTable_; // hashtable storing the operator views (keys = view names, values = views).

      viewLabel_t defaultViewLabel_;  // label of the view associated with inital Matrix construction
      viewLabel_t currentViewLabel_;  // label of the current view

  }; //class Matrix

} //namespace Xpetra

#define XPETRA_MATRIX_SHORT
#endif //XPETRA_MATRIX_DECL_HPP
