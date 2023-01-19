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

#ifndef XPETRA_CRSMATRIXWRAP_DEF_HPP
#define XPETRA_CRSMATRIXWRAP_DEF_HPP

#include <KokkosCompat_DefaultNode.hpp>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Hashtable.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_MultiVector.hpp"
#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_CrsMatrixFactory.hpp"

#include "Xpetra_Matrix.hpp"

#include "Xpetra_CrsMatrixWrap_decl.hpp"

namespace Xpetra {
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::CrsMatrixWrap (const RCP<const Map>& rowMap)
    : finalDefaultView_ (false)
  {
    matrixData_ = CrsMatrixFactory::Build (rowMap);
    CreateDefaultView ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::CrsMatrixWrap (const RCP<const Map>& rowMap,
                 size_t maxNumEntriesPerRow)
    : finalDefaultView_ (false)
  {
    matrixData_ = CrsMatrixFactory::Build (rowMap, maxNumEntriesPerRow);
    CreateDefaultView ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::CrsMatrixWrap (const RCP<const Map>& rowMap,
                 const ArrayRCP<const size_t>& NumEntriesPerRowToAlloc)
    : finalDefaultView_ (false)
  {
    matrixData_ = CrsMatrixFactory::Build(rowMap, NumEntriesPerRowToAlloc);
    CreateDefaultView ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::CrsMatrixWrap(const RCP<const Map> &rowMap, const RCP<const Map>& colMap, size_t maxNumEntriesPerRow)
    : finalDefaultView_(false)
  {
    // Set matrix data
    matrixData_ = CrsMatrixFactory::Build(rowMap, colMap, maxNumEntriesPerRow);

    // Default view
    CreateDefaultView();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::CrsMatrixWrap(const RCP<const Map> &rowMap, const RCP<const Map>& colMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc)
    : finalDefaultView_(false)
  {
    // Set matrix data
    matrixData_ = CrsMatrixFactory::Build(rowMap, colMap, NumEntriesPerRowToAlloc);

    // Default view
    CreateDefaultView();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::CrsMatrixWrap(const RCP<const Map> &rowMap, const RCP<const Map>& colMap, const local_matrix_type& lclMatrix, const Teuchos::RCP<Teuchos::ParameterList>& params)
    : finalDefaultView_(false)
  {
    // Set matrix data
    matrixData_ = CrsMatrixFactory::Build(rowMap, colMap, lclMatrix, params);

    // Default view
    CreateDefaultView();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::CrsMatrixWrap(const local_matrix_type& lclMatrix, const RCP<const Map> &rowMap, const RCP<const Map>& colMap,
                const RCP<const Map>& domainMap, const RCP<const Map>& rangeMap,
                const Teuchos::RCP<Teuchos::ParameterList>& params)
    : finalDefaultView_(false)
  {
    // Set matrix data
    matrixData_ = CrsMatrixFactory::Build(lclMatrix, rowMap, colMap, domainMap, rangeMap, params);

    // Default view
    CreateDefaultView();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::CrsMatrixWrap(RCP<CrsMatrix> matrix)
    : finalDefaultView_(matrix->isFillComplete())
  {
    // Set matrix data
    matrixData_ = matrix;

    // Default view
    CreateDefaultView();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::CrsMatrixWrap(const RCP<const CrsGraph>& graph, const RCP<ParameterList>& paramList)
    : finalDefaultView_(false)
  {
    // Set matrix data
    matrixData_ = CrsMatrixFactory::Build(graph, paramList);

    // Default view
    CreateDefaultView();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~CrsMatrixWrap() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &cols, const ArrayView<const Scalar> &vals) {
    matrixData_->insertGlobalValues(globalRow, cols, vals);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &cols, const ArrayView<const Scalar> &vals) {
    matrixData_->insertLocalValues(localRow, cols, vals);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceGlobalValues(GlobalOrdinal globalRow,
                           const ArrayView<const GlobalOrdinal> &cols,
                           const ArrayView<const Scalar>        &vals) { matrixData_->replaceGlobalValues(globalRow, cols, vals); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceLocalValues(LocalOrdinal localRow,
                          const ArrayView<const LocalOrdinal> &cols,
                          const ArrayView<const Scalar>       &vals) { matrixData_->replaceLocalValues(localRow, cols, vals); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setAllToScalar(const Scalar &alpha) { matrixData_->setAllToScalar(alpha); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::scale(const Scalar &alpha) {
    matrixData_->scale(alpha);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::resumeFill(const RCP< ParameterList > &params) {
    matrixData_->resumeFill(params);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::fillComplete(const RCP<const Map> &domainMap, const RCP<const Map> &rangeMap, const RCP<Teuchos::ParameterList> &params) {
    matrixData_->fillComplete(domainMap, rangeMap, params);

    // Update default view with the colMap because colMap can be <tt>null</tt> until fillComplete() is called.
    updateDefaultView();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::fillComplete(const RCP<ParameterList> &params) {
    matrixData_->fillComplete(params);

    // Update default view with the colMap because colMap can be <tt>null</tt> until fillComplete() is called.
    updateDefaultView();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumRows() const {
    return matrixData_->getGlobalNumRows();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumCols() const {
    return matrixData_->getGlobalNumCols();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalNumRows() const {
    return matrixData_->getLocalNumRows();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumEntries() const {
    return matrixData_->getGlobalNumEntries();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalNumEntries() const {
    return matrixData_->getLocalNumEntries();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNumEntriesInLocalRow(LocalOrdinal localRow) const {
    return matrixData_->getNumEntriesInLocalRow(localRow);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const {
    return matrixData_->getNumEntriesInGlobalRow(globalRow);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getGlobalMaxNumRowEntries() const {
    return matrixData_->getGlobalMaxNumRowEntries();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalMaxNumRowEntries() const {
    return matrixData_->getLocalMaxNumRowEntries();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isLocallyIndexed() const {
    return matrixData_->isLocallyIndexed();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isGloballyIndexed() const {
    return matrixData_->isGloballyIndexed();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isFillComplete() const {
    return matrixData_->isFillComplete();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalRowCopy(LocalOrdinal LocalRow,
                       const ArrayView<LocalOrdinal> &Indices,
                       const ArrayView<Scalar> &Values,
                       size_t &NumEntries
                       ) const {
    matrixData_->getLocalRowCopy(LocalRow, Indices, Values, NumEntries);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &indices, ArrayView<const Scalar> &values) const {
     matrixData_->getGlobalRowView(GlobalRow, indices, values);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices, ArrayView<const Scalar> &values) const {
     matrixData_->getLocalRowView(LocalRow, indices, values);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalDiagCopy(Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const {
    matrixData_->getLocalDiagCopy(diag);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalDiagOffsets(Teuchos::ArrayRCP<size_t> &offsets) const {
    matrixData_->getLocalDiagOffsets(offsets);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalDiagCopy(Xpetra::Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &diag, const Teuchos::ArrayView<const size_t> &offsets) const {
    matrixData_->getLocalDiagCopy(diag,offsets);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename ScalarTraits<Scalar>::magnitudeType CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getFrobeniusNorm() const {
    return matrixData_->getFrobeniusNorm();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::leftScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
    matrixData_->leftScale(x);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::rightScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
    matrixData_->rightScale(x);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::haveGlobalConstants() const {
    return matrixData_->haveGlobalConstants();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::apply(const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                   Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                   Teuchos::ETransp mode,
                   Scalar alpha,
                   Scalar beta) const {

    matrixData_->apply(X,Y,mode,alpha,beta);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::apply(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &X,
                  MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &Y,
                  Teuchos::ETransp mode,
                  Scalar alpha,
                  Scalar beta,
                  bool sumInterfaceValues,
                  const RCP<Import<LocalOrdinal, GlobalOrdinal, Node> >& regionInterfaceImporter,
                  const Teuchos::ArrayRCP<LocalOrdinal>& regionInterfaceLIDs
  ) const{
      matrixData_->apply(X,Y,mode,alpha,beta,sumInterfaceValues,regionInterfaceImporter,regionInterfaceLIDs);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getDomainMap() const {
    return matrixData_->getDomainMap();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getRangeMap() const {
    return matrixData_->getRangeMap();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>> & CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getColMap() const { return getColMap(Matrix::GetCurrentViewLabel()); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>> & CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getColMap(viewLabel_t viewLabel) const {
    TEUCHOS_TEST_FOR_EXCEPTION(Matrix::operatorViewTable_.containsKey(viewLabel) == false, Xpetra::Exceptions::RuntimeError, "Xpetra::Matrix.GetColMap(): view '" + viewLabel + "' does not exist.");
    updateDefaultView(); // If CrsMatrix::fillComplete() have been used instead of CrsMatrixWrap::fillComplete(), the default view is updated.
    return Matrix::operatorViewTable_.get(viewLabel)->GetColMap();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::removeEmptyProcessesInPlace(const Teuchos::RCP<const Map>& newMap) {
    matrixData_->removeEmptyProcessesInPlace(newMap);
    this->operatorViewTable_.get(this->GetCurrentViewLabel())->SetRowMap(matrixData_->getRowMap());
    this->operatorViewTable_.get(this->GetCurrentViewLabel())->SetColMap(matrixData_->getColMap());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const Teuchos::RCP< const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getMap() const {
    return matrixData_->getMap();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::doImport(const Matrix &source,
                const Xpetra::Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM) {
    const CrsMatrixWrap & sourceWrp = dynamic_cast<const CrsMatrixWrap &>(source);
    matrixData_->doImport(*sourceWrp.getCrsMatrix(), importer, CM);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::doExport(const Matrix &dest,
                const Xpetra::Import< LocalOrdinal, GlobalOrdinal, Node >& importer, CombineMode CM) {
    const CrsMatrixWrap & destWrp = dynamic_cast<const CrsMatrixWrap &>(dest);
    matrixData_->doExport(*destWrp.getCrsMatrix(), importer, CM);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::doImport(const Matrix &source,
                const Xpetra::Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM) {
    const CrsMatrixWrap & sourceWrp = dynamic_cast<const CrsMatrixWrap &>(source);
    matrixData_->doImport(*sourceWrp.getCrsMatrix(), exporter, CM);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::doExport(const Matrix &dest,
                const Xpetra::Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM) {
    const CrsMatrixWrap & destWrp = dynamic_cast<const CrsMatrixWrap &>(dest);
    matrixData_->doExport(*destWrp.getCrsMatrix(), exporter, CM);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::description() const {
    return "Xpetra::CrsMatrixWrap";
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    //     Teuchos::EVerbosityLevel vl = verbLevel;
    //     if (vl == VERB_DEFAULT) vl = VERB_LOW;
    //     RCP<const Comm<int> > comm = this->getComm();
    //     const int myImageID = comm->getRank(),
    //       numImages = comm->getSize();

    //     if (myImageID == 0) out << this->description() << std::endl;

    matrixData_->describe(out,verbLevel);

    // Teuchos::OSTab tab(out);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setObjectLabel( const std::string &objectLabel ) {
    Teuchos::LabeledObject::setObjectLabel(objectLabel);
    matrixData_->setObjectLabel(objectLabel);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type::HostMirror
  CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalMatrixHost () const {
    return matrixData_->getLocalMatrixHost();
  }
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type
  CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalMatrixDevice () const {
    return matrixData_->getLocalMatrixDevice();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::hasCrsGraph() const {return true;}


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getCrsGraph() const { return matrixData_->getCrsGraph(); }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getCrsMatrix() const {  return matrixData_; }

// Default view is created after fillComplete()
  // Because ColMap might not be available before fillComplete().
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::CreateDefaultView() {

    // Create default view
    this->defaultViewLabel_ = "point";
    this->CreateView(this->GetDefaultViewLabel(), matrixData_->getRowMap(), matrixData_->getColMap());

    // Set current view
    this->currentViewLabel_ = this->GetDefaultViewLabel();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::updateDefaultView() const {
    if ((finalDefaultView_ == false) &&  matrixData_->isFillComplete() ) {
      // Update default view with the colMap
      Matrix::operatorViewTable_.get(Matrix::GetDefaultViewLabel())->SetColMap(matrixData_->getColMap());
      finalDefaultView_ = true;
    }
  }


  // Expert only
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceCrsMatrix(RCP<CrsMatrix> & M) {
    // Clear the old view table
    Teuchos::Hashtable<viewLabel_t, RCP<MatrixView> > dummy_table;
    Matrix::operatorViewTable_ = dummy_table;

    finalDefaultView_ = M->isFillComplete();
    // Set matrix data
    matrixData_ = M;
    

    // Default view
    CreateDefaultView();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetStorageBlockSize() const {
    return matrixData_->GetStorageBlockSize();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>::residual(
            const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > & X, 
            const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > & B,
            MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > & R) const {
    matrixData_->residual(X,B,R);
  }


} //namespace Xpetra

#endif //ifndef XPETRA_CRSMATRIXWRAP_DEF_HPP
