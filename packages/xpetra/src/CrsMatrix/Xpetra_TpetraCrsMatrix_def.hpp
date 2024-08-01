// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_TPETRACRSMATRIX_DEF_HPP
#define XPETRA_TPETRACRSMATRIX_DEF_HPP

#include <Xpetra_MultiVectorFactory.hpp>
#include "Xpetra_TpetraCrsMatrix_decl.hpp"
#include "Tpetra_Details_residual.hpp"

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsMatrix(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &rowMap, size_t maxNumEntriesPerRow, const Teuchos::RCP<Teuchos::ParameterList> &params)
  : mtx_(Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rowMap), maxNumEntriesPerRow, params))) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsMatrix(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &rowMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, const Teuchos::RCP<Teuchos::ParameterList> &params)
  : mtx_(Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rowMap), NumEntriesPerRowToAlloc(), params))) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsMatrix(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &rowMap, const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &colMap, size_t maxNumEntriesPerRow, const Teuchos::RCP<Teuchos::ParameterList> &params)
  : mtx_(Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rowMap), toTpetra(colMap), maxNumEntriesPerRow, params))) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsMatrix(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &rowMap, const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &colMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, const Teuchos::RCP<Teuchos::ParameterList> &params)
  : mtx_(Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rowMap), toTpetra(colMap), NumEntriesPerRowToAlloc(), params))) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsMatrix(const Teuchos::RCP<const CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> &graph, const Teuchos::RCP<Teuchos::ParameterList> &params)
  : mtx_(Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(toTpetra(graph), params))) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsMatrix(const Teuchos::RCP<const CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> &graph, typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type::values_type &values, const Teuchos::RCP<Teuchos::ParameterList> &params)
  : mtx_(Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(toTpetra(graph), values, params))) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsMatrix(const Teuchos::RCP<const CrsMatrix> &sourceMatrix,
                                                                            const Import<LocalOrdinal, GlobalOrdinal, Node> &importer,
                                                                            const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &domainMap,
                                                                            const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &rangeMap,
                                                                            const Teuchos::RCP<Teuchos::ParameterList> &params) {
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> MyTpetraCrsMatrix;
  XPETRA_DYNAMIC_CAST(const TpetraCrsMatrixClass, *sourceMatrix, tSourceMatrix, "Xpetra::TpetraCrsMatrix constructor only accepts Xpetra::TpetraCrsMatrix as the input argument.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> v = tSourceMatrix.getTpetra_CrsMatrix();

  RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> myDomainMap = domainMap != Teuchos::null ? toTpetra(domainMap) : Teuchos::null;
  RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> myRangeMap  = rangeMap != Teuchos::null ? toTpetra(rangeMap) : Teuchos::null;
  mtx_                                                                  = Tpetra::importAndFillCompleteCrsMatrix<MyTpetraCrsMatrix>(tSourceMatrix.getTpetra_CrsMatrix(), toTpetra(importer), myDomainMap, myRangeMap, params);
  bool restrictComm                                                     = false;
  if (!params.is_null()) restrictComm = params->get("Restrict Communicator", restrictComm);
  if (restrictComm && mtx_->getRowMap().is_null()) mtx_ = Teuchos::null;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsMatrix(const Teuchos::RCP<const CrsMatrix> &sourceMatrix,
                                                                            const Export<LocalOrdinal, GlobalOrdinal, Node> &exporter,
                                                                            const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &domainMap,
                                                                            const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &rangeMap,
                                                                            const Teuchos::RCP<Teuchos::ParameterList> &params) {
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> MyTpetraCrsMatrix;
  XPETRA_DYNAMIC_CAST(const TpetraCrsMatrixClass, *sourceMatrix, tSourceMatrix, "Xpetra::TpetraCrsMatrix constructor only accepts Xpetra::TpetraCrsMatrix as the input argument.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> v = tSourceMatrix.getTpetra_CrsMatrix();

  RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> myDomainMap = domainMap != Teuchos::null ? toTpetra(domainMap) : Teuchos::null;
  RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> myRangeMap  = rangeMap != Teuchos::null ? toTpetra(rangeMap) : Teuchos::null;
  mtx_                                                                  = Tpetra::exportAndFillCompleteCrsMatrix<MyTpetraCrsMatrix>(tSourceMatrix.getTpetra_CrsMatrix(), toTpetra(exporter), myDomainMap, myRangeMap, params);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsMatrix(const Teuchos::RCP<const CrsMatrix> &sourceMatrix,
                                                                            const Import<LocalOrdinal, GlobalOrdinal, Node> &RowImporter,
                                                                            const Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node>> DomainImporter,
                                                                            const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &domainMap,
                                                                            const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &rangeMap,
                                                                            const Teuchos::RCP<Teuchos::ParameterList> &params) {
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> MyTpetraCrsMatrix;
  XPETRA_DYNAMIC_CAST(const TpetraCrsMatrixClass, *sourceMatrix, tSourceMatrix, "Xpetra::TpetraCrsMatrix constructor only accepts Xpetra::TpetraCrsMatrix as the input argument.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> v = tSourceMatrix.getTpetra_CrsMatrix();

  RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> myDomainMap = domainMap != Teuchos::null ? toTpetra(domainMap) : Teuchos::null;
  RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> myRangeMap  = rangeMap != Teuchos::null ? toTpetra(rangeMap) : Teuchos::null;

  mtx_              = Tpetra::importAndFillCompleteCrsMatrix<MyTpetraCrsMatrix>(tSourceMatrix.getTpetra_CrsMatrix(), toTpetra(RowImporter), toTpetra(*DomainImporter), myDomainMap, myRangeMap, params);
  bool restrictComm = false;
  if (!params.is_null()) restrictComm = params->get("Restrict Communicator", restrictComm);
  if (restrictComm && mtx_->getRowMap().is_null()) mtx_ = Teuchos::null;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsMatrix(const Teuchos::RCP<const CrsMatrix> &sourceMatrix,
                                                                            const Export<LocalOrdinal, GlobalOrdinal, Node> &RowExporter,
                                                                            const Teuchos::RCP<const Export<LocalOrdinal, GlobalOrdinal, Node>> DomainExporter,
                                                                            const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &domainMap,
                                                                            const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &rangeMap,
                                                                            const Teuchos::RCP<Teuchos::ParameterList> &params) {
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> MyTpetraCrsMatrix;
  XPETRA_DYNAMIC_CAST(const TpetraCrsMatrixClass, *sourceMatrix, tSourceMatrix, "Xpetra::TpetraCrsMatrix constructor only accepts Xpetra::TpetraCrsMatrix as the input argument.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> v = tSourceMatrix.getTpetra_CrsMatrix();

  RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> myDomainMap = domainMap != Teuchos::null ? toTpetra(domainMap) : Teuchos::null;
  RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> myRangeMap  = rangeMap != Teuchos::null ? toTpetra(rangeMap) : Teuchos::null;

  mtx_ = Tpetra::exportAndFillCompleteCrsMatrix<MyTpetraCrsMatrix>(tSourceMatrix.getTpetra_CrsMatrix(), toTpetra(RowExporter), toTpetra(*DomainExporter), myDomainMap, myRangeMap, params);
}

///////////////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_XPETRA_TPETRA
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsMatrix(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &rowMap,
                                                                            const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &colMap,
                                                                            const local_matrix_type &lclMatrix,
                                                                            const Teuchos::RCP<Teuchos::ParameterList> &params)
  : mtx_(Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rowMap), toTpetra(colMap), lclMatrix, params))) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsMatrix(
    const local_matrix_type &lclMatrix,
    const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &rowMap,
    const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &colMap,
    const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &domainMap,
    const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &rangeMap,
    const Teuchos::RCP<Teuchos::ParameterList> &params)
  : mtx_(Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(lclMatrix, toTpetra(rowMap), toTpetra(colMap), toTpetra(domainMap), toTpetra(rangeMap), params))) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsMatrix(
    const local_matrix_type &lclMatrix,
    const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &rowMap,
    const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &colMap,
    const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &domainMap,
    const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &rangeMap,
    const Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node>> &importer,
    const Teuchos::RCP<const Export<LocalOrdinal, GlobalOrdinal, Node>> &exporter,
    const Teuchos::RCP<Teuchos::ParameterList> &params)
  : mtx_(Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(lclMatrix, toTpetra(rowMap), toTpetra(colMap), toTpetra(domainMap), toTpetra(rangeMap), toTpetra(importer), toTpetra(exporter), params))) {}
#endif

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~TpetraCrsMatrix() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &cols, const ArrayView<const Scalar> &vals) {
  XPETRA_MONITOR("TpetraCrsMatrix::insertGlobalValues");
  mtx_->insertGlobalValues(globalRow, cols, vals);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &cols, const ArrayView<const Scalar> &vals) {
  XPETRA_MONITOR("TpetraCrsMatrix::insertLocalValues");
  mtx_->insertLocalValues(localRow, cols, vals);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &cols, const ArrayView<const Scalar> &vals) {
  XPETRA_MONITOR("TpetraCrsMatrix::replaceGlobalValues");
  mtx_->replaceGlobalValues(globalRow, cols, vals);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceLocalValues(LocalOrdinal localRow,
                                                                                    const ArrayView<const LocalOrdinal> &cols,
                                                                                    const ArrayView<const Scalar> &vals) {
  XPETRA_MONITOR("TpetraCrsMatrix::replaceLocalValues");
  typedef typename ArrayView<const LocalOrdinal>::size_type size_type;
  const LocalOrdinal numValid =
      mtx_->replaceLocalValues(localRow, cols, vals);
  TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_type>(numValid) != cols.size(), std::runtime_error,
      "replaceLocalValues returned " << numValid << " != cols.size() = " << cols.size() << ".");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setAllToScalar(const Scalar &alpha) {
  XPETRA_MONITOR("TpetraCrsMatrix::setAllToScalar");
  mtx_->setAllToScalar(alpha);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::scale(const Scalar &alpha) {
  XPETRA_MONITOR("TpetraCrsMatrix::scale");
  mtx_->scale(alpha);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::allocateAllValues(size_t numNonZeros, ArrayRCP<size_t> &rowptr, ArrayRCP<LocalOrdinal> &colind, ArrayRCP<Scalar> &values) {
  XPETRA_MONITOR("TpetraCrsMatrix::allocateAllValues");
  rowptr.resize(getLocalNumRows() + 1);
  colind.resize(numNonZeros);
  values.resize(numNonZeros);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setAllValues(const ArrayRCP<size_t> &rowptr, const ArrayRCP<LocalOrdinal> &colind, const ArrayRCP<Scalar> &values) {
  XPETRA_MONITOR("TpetraCrsMatrix::setAllValues");
  mtx_->setAllValues(rowptr, colind, values);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getAllValues(ArrayRCP<const size_t> &rowptr, ArrayRCP<const LocalOrdinal> &colind, ArrayRCP<const Scalar> &values) const {
  XPETRA_MONITOR("TpetraCrsMatrix::getAllValues");
  // TODO: Change getAllValues interface to return Kokkos views,
  // TODO: or add getLocalRowPtrsHost, getLocalIndicesHost, etc., interfaces
  // TODO: and use them in MueLu
  rowptr = Kokkos::Compat::persistingView(mtx_->getLocalRowPtrsHost());
  colind = Kokkos::Compat::persistingView(mtx_->getLocalIndicesHost());
  values = Teuchos::arcp_reinterpret_cast<const Scalar>(
      Kokkos::Compat::persistingView(mtx_->getLocalValuesHost(
          Tpetra::Access::ReadOnly)));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getAllValues(ArrayRCP<Scalar> &values) {
  XPETRA_MONITOR("TpetraCrsMatrix::getAllValues");
  // TODO: Change getAllValues interface to return Kokkos view,
  // TODO: or add getLocalValuesDevice() interfaces
  values = Teuchos::arcp_reinterpret_cast<Scalar>(
      Kokkos::Compat::persistingView(mtx_->getLocalValuesDevice(
          Tpetra::Access::ReadWrite)));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::haveGlobalConstants() const { return mtx_->haveGlobalConstants(); }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::resumeFill(const RCP<ParameterList> &params) {
  XPETRA_MONITOR("TpetraCrsMatrix::resumeFill");
  mtx_->resumeFill(params);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::fillComplete(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &domainMap, const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &rangeMap, const RCP<ParameterList> &params) {
  XPETRA_MONITOR("TpetraCrsMatrix::fillComplete");
  mtx_->fillComplete(toTpetra(domainMap), toTpetra(rangeMap), params);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::fillComplete(const RCP<ParameterList> &params) {
  XPETRA_MONITOR("TpetraCrsMatrix::fillComplete");
  mtx_->fillComplete(params);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceDomainMapAndImporter(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &newDomainMap, Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node>> &newImporter) {
  XPETRA_MONITOR("TpetraCrsMatrix::replaceDomainMapAndImporter");
  XPETRA_DYNAMIC_CAST(const TpetraImportClass, *newImporter, tImporter, "Xpetra::TpetraCrsMatrix::replaceDomainMapAndImporter only accepts Xpetra::TpetraImport.");
  RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>> myImport = tImporter.getTpetra_Import();
  mtx_->replaceDomainMapAndImporter(toTpetra(newDomainMap), myImport);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::expertStaticFillComplete(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &domainMap,
                                                                                          const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &rangeMap,
                                                                                          const RCP<const Import<LocalOrdinal, GlobalOrdinal, Node>> &importer,
                                                                                          const RCP<const Export<LocalOrdinal, GlobalOrdinal, Node>> &exporter,
                                                                                          const RCP<ParameterList> &params) {
  XPETRA_MONITOR("TpetraCrsMatrix::expertStaticFillComplete");
  RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>> myImport;
  RCP<const Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>> myExport;

  if (importer != Teuchos::null) {
    XPETRA_DYNAMIC_CAST(const TpetraImportClass, *importer, tImporter, "Xpetra::TpetraCrsMatrix::expertStaticFillComplete only accepts Xpetra::TpetraImport.");
    myImport = tImporter.getTpetra_Import();
  }
  if (exporter != Teuchos::null) {
    XPETRA_DYNAMIC_CAST(const TpetraExportClass, *exporter, tExporter, "Xpetra::TpetraCrsMatrix::expertStaticFillComplete only accepts Xpetra::TpetraExport.");
    myExport = tExporter.getTpetra_Export();
  }

  mtx_->expertStaticFillComplete(toTpetra(domainMap), toTpetra(rangeMap), myImport, myExport, params);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRowMap() const {
  XPETRA_MONITOR("TpetraCrsMatrix::getRowMap");
  return toXpetra(mtx_->getRowMap());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getColMap() const {
  XPETRA_MONITOR("TpetraCrsMatrix::getColMap");
  return toXpetra(mtx_->getColMap());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getCrsGraph() const {
  XPETRA_MONITOR("TpetraCrsMatrix::getCrsGraph");
  return toXpetra(mtx_->getCrsGraph());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalNumRows() const {
  XPETRA_MONITOR("TpetraCrsMatrix::getGlobalNumRows");
  return mtx_->getGlobalNumRows();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalNumCols() const {
  XPETRA_MONITOR("TpetraCrsMatrix::getGlobalNumCols");
  return mtx_->getGlobalNumCols();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalNumRows() const {
  XPETRA_MONITOR("TpetraCrsMatrix::getLocalNumRows");
  return mtx_->getLocalNumRows();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalNumCols() const {
  XPETRA_MONITOR("TpetraCrsMatrix::getLocalNumCols");
  return mtx_->getLocalNumCols();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalNumEntries() const {
  XPETRA_MONITOR("TpetraCrsMatrix::getGlobalNumEntries");
  return mtx_->getGlobalNumEntries();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalNumEntries() const {
  XPETRA_MONITOR("TpetraCrsMatrix::getLocalNumEntries");
  return mtx_->getLocalNumEntries();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNumEntriesInLocalRow(LocalOrdinal localRow) const {
  XPETRA_MONITOR("TpetraCrsMatrix::getNumEntriesInLocalRow");
  return mtx_->getNumEntriesInLocalRow(localRow);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const {
  XPETRA_MONITOR("TpetraCrsMatrix::getNumEntriesInGlobalRow");
  return mtx_->getNumEntriesInGlobalRow(globalRow);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalMaxNumRowEntries() const {
  XPETRA_MONITOR("TpetraCrsMatrix::getGlobalMaxNumRowEntries");
  return mtx_->getGlobalMaxNumRowEntries();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalMaxNumRowEntries() const {
  XPETRA_MONITOR("TpetraCrsMatrix::getLocalMaxNumRowEntries");
  return mtx_->getLocalMaxNumRowEntries();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isLocallyIndexed() const {
  XPETRA_MONITOR("TpetraCrsMatrix::isLocallyIndexed");
  return mtx_->isLocallyIndexed();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isGloballyIndexed() const {
  XPETRA_MONITOR("TpetraCrsMatrix::isGloballyIndexed");
  return mtx_->isGloballyIndexed();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isFillComplete() const {
  XPETRA_MONITOR("TpetraCrsMatrix::isFillComplete");
  return mtx_->isFillComplete();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isFillActive() const {
  XPETRA_MONITOR("TpetraCrsMatrix::isFillActive");
  return mtx_->isFillActive();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename ScalarTraits<Scalar>::magnitudeType TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getFrobeniusNorm() const {
  XPETRA_MONITOR("TpetraCrsMatrix::getFrobeniusNorm");
  return mtx_->getFrobeniusNorm();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::supportsRowViews() const {
  XPETRA_MONITOR("TpetraCrsMatrix::supportsRowViews");
  return mtx_->supportsRowViews();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalRowCopy(LocalOrdinal LocalRow, const ArrayView<LocalOrdinal> &Indices, const ArrayView<Scalar> &Values, size_t &NumEntries) const {
  XPETRA_MONITOR("TpetraCrsMatrix::getLocalRowCopy");
  typename Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::nonconst_local_inds_host_view_type indices("indices", Indices.size());
  typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::nonconst_values_host_view_type values("values", Values.size());

  mtx_->getLocalRowCopy(LocalRow, indices, values, NumEntries);
  for (size_t i = 0; i < NumEntries; ++i) {
    Indices[i] = indices(i);
    Values[i]  = values(i);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalRowCopy(GlobalOrdinal GlobalRow, const ArrayView<GlobalOrdinal> &Indices, const ArrayView<Scalar> &Values, size_t &NumEntries) const {
  XPETRA_MONITOR("TpetraCrsMatrix::getGlobalRowCopy");
  typename Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::nonconst_global_inds_host_view_type indices("indices", Indices.size());
  typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::nonconst_values_host_view_type values("values", Values.size());

  mtx_->getGlobalRowCopy(GlobalRow, indices, values, NumEntries);
  for (size_t i = 0; i < NumEntries; ++i) {
    Indices[i] = indices(i);
    Values[i]  = values(i);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &Indices, ArrayView<const Scalar> &Values) const {
  XPETRA_MONITOR("TpetraCrsMatrix::getGlobalRowView");
  typename Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::global_inds_host_view_type indices;
  typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::values_host_view_type values;

  mtx_->getGlobalRowView(GlobalRow, indices, values);
  Indices = ArrayView<const GlobalOrdinal>(indices.data(), indices.extent(0));
  Values  = ArrayView<const Scalar>(reinterpret_cast<const Scalar *>(values.data()), values.extent(0));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &Indices, ArrayView<const Scalar> &Values) const {
  XPETRA_MONITOR("TpetraCrsMatrix::getLocalRowView");
  typename Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_inds_host_view_type indices;
  typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::values_host_view_type values;

  mtx_->getLocalRowView(LocalRow, indices, values);
  Indices = ArrayView<const LocalOrdinal>(indices.data(), indices.extent(0));
  Values  = ArrayView<const Scalar>(reinterpret_cast<const Scalar *>(values.data()), values.extent(0));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::apply(const MultiVector &X, MultiVector &Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const {
  XPETRA_MONITOR("TpetraCrsMatrix::apply");
  mtx_->apply(toTpetra(X), toTpetra(Y), mode, alpha, beta);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::apply(const MultiVector &X, MultiVector &Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta, bool sumInterfaceValues, const RCP<Import<LocalOrdinal, GlobalOrdinal, Node>> &regionInterfaceImporter, const Teuchos::ArrayRCP<LocalOrdinal> &regionInterfaceLIDs) const {
  XPETRA_MONITOR("TpetraCrsMatrix::apply(region)");
  RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> regionInterfaceMap = regionInterfaceImporter->getTargetMap();
  mtx_->localApply(toTpetra(X), toTpetra(Y), mode, alpha, beta);
  if (sumInterfaceValues) {
    // preform communication to propagate local interface
    // values to all the processor that share interfaces.
    RCP<MultiVector> matvecInterfaceTmp = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(regionInterfaceMap, 1);
    matvecInterfaceTmp->doImport(Y, *regionInterfaceImporter, INSERT);

    // sum all contributions to interface values
    // on all ranks
    ArrayRCP<Scalar> YData         = Y.getDataNonConst(0);
    ArrayRCP<Scalar> interfaceData = matvecInterfaceTmp->getDataNonConst(0);
    for (LocalOrdinal interfaceIdx = 0; interfaceIdx < static_cast<LocalOrdinal>(interfaceData.size()); ++interfaceIdx) {
      YData[regionInterfaceLIDs[interfaceIdx]] += interfaceData[interfaceIdx];
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getDomainMap() const {
  XPETRA_MONITOR("TpetraCrsMatrix::getDomainMap");
  return toXpetra(mtx_->getDomainMap());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRangeMap() const {
  XPETRA_MONITOR("TpetraCrsMatrix::getRangeMap");
  return toXpetra(mtx_->getRangeMap());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
  XPETRA_MONITOR("TpetraCrsMatrix::description");
  return mtx_->description();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
  XPETRA_MONITOR("TpetraCrsMatrix::describe");
  mtx_->describe(out, verbLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setObjectLabel(const std::string &objectLabel) {
  XPETRA_MONITOR("TpetraCrsMatrix::setObjectLabel");
  Teuchos::LabeledObject::setObjectLabel(objectLabel);
  mtx_->setObjectLabel(objectLabel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsMatrix(const TpetraCrsMatrix &matrix)
  : mtx_(Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(*(matrix.mtx_), Teuchos::Copy))) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalDiagCopy(Vector &diag) const {
  XPETRA_MONITOR("TpetraCrsMatrix::getLocalDiagCopy");
  XPETRA_DYNAMIC_CAST(TpetraVectorClass, diag, tDiag, "Xpetra::TpetraCrsMatrix.getLocalDiagCopy() only accept Xpetra::TpetraVector as input arguments.");
  mtx_->getLocalDiagCopy(*tDiag.getTpetra_Vector());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalDiagOffsets(Teuchos::ArrayRCP<size_t> &offsets) const {
  XPETRA_MONITOR("TpetraCrsMatrix::getLocalDiagOffsets");
  mtx_->getLocalDiagOffsets(offsets);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalDiagCopy(Vector &diag, const Teuchos::ArrayView<const size_t> &offsets) const {
  XPETRA_MONITOR("TpetraCrsMatrix::getLocalDiagCopy");
  mtx_->getLocalDiagCopy(*(toTpetra(diag)), offsets);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalDiagCopy(Vector &diag, const Kokkos::View<const size_t *, typename Node::device_type, Kokkos::MemoryUnmanaged> &offsets) const {
  XPETRA_MONITOR("TpetraCrsMatrix::getLocalDiagCopy");
  mtx_->getLocalDiagCopy(*(toTpetra(diag)), offsets);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceDiag(const Vector &diag) {
  XPETRA_MONITOR("TpetraCrsMatrix::replaceDiag");
  Tpetra::replaceDiagonalCrsMatrix(*mtx_, *(toTpetra(diag)));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::leftScale(const Vector &x) {
  XPETRA_MONITOR("TpetraCrsMatrix::leftScale");
  mtx_->leftScale(*(toTpetra(x)));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::rightScale(const Vector &x) {
  XPETRA_MONITOR("TpetraCrsMatrix::rightScale");
  mtx_->rightScale(*(toTpetra(x)));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getMap() const {
  XPETRA_MONITOR("TpetraCrsMatrix::getMap");
  return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(mtx_->getMap()));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doImport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &source,
                                                                          const Import<LocalOrdinal, GlobalOrdinal, Node> &importer, CombineMode CM) {
  XPETRA_MONITOR("TpetraCrsMatrix::doImport");

  XPETRA_DYNAMIC_CAST(const TpetraCrsMatrixClass, source, tSource, "Xpetra::TpetraCrsMatrix::doImport only accept Xpetra::TpetraCrsMatrix as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> v = tSource.getTpetra_CrsMatrix();
  // mtx_->doImport(toTpetraCrsMatrix(source), *tImporter.getTpetra_Import(), toTpetra(CM));
  mtx_->doImport(*v, toTpetra(importer), toTpetra(CM));
}

//! Export.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doExport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &dest,
                                                                          const Import<LocalOrdinal, GlobalOrdinal, Node> &importer, CombineMode CM) {
  XPETRA_MONITOR("TpetraCrsMatrix::doExport");

  XPETRA_DYNAMIC_CAST(const TpetraCrsMatrixClass, dest, tDest, "Xpetra::TpetraCrsMatrix::doImport only accept Xpetra::TpetraCrsMatrix as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> v = tDest.getTpetra_CrsMatrix();
  mtx_->doExport(*v, toTpetra(importer), toTpetra(CM));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doImport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &source,
                                                                          const Export<LocalOrdinal, GlobalOrdinal, Node> &exporter, CombineMode CM) {
  XPETRA_MONITOR("TpetraCrsMatrix::doImport");

  XPETRA_DYNAMIC_CAST(const TpetraCrsMatrixClass, source, tSource, "Xpetra::TpetraCrsMatrix::doImport only accept Xpetra::TpetraCrsMatrix as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> v = tSource.getTpetra_CrsMatrix();
  mtx_->doImport(*v, toTpetra(exporter), toTpetra(CM));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doExport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &dest,
                                                                          const Export<LocalOrdinal, GlobalOrdinal, Node> &exporter, CombineMode CM) {
  XPETRA_MONITOR("TpetraCrsMatrix::doExport");

  XPETRA_DYNAMIC_CAST(const TpetraCrsMatrixClass, dest, tDest, "Xpetra::TpetraCrsMatrix::doImport only accept Xpetra::TpetraCrsMatrix as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> v = tDest.getTpetra_CrsMatrix();
  mtx_->doExport(*v, toTpetra(exporter), toTpetra(CM));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::removeEmptyProcessesInPlace(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &newMap) {
  XPETRA_MONITOR("TpetraCrsMatrix::removeEmptyProcessesInPlace");
  mtx_->removeEmptyProcessesInPlace(toTpetra(newMap));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::hasMatrix() const {
  return !mtx_.is_null();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsMatrix(const Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &mtx)
  : mtx_(mtx) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getTpetra_CrsMatrix() const { return mtx_; }

//! Get the underlying Tpetra matrix
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getTpetra_CrsMatrixNonConst() const { return mtx_; }  // TODO: remove

//! Compute a residual R = B - (*this) * X
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::residual(const MultiVector &X,
                                                                          const MultiVector &B,
                                                                          MultiVector &R) const {
  Tpetra::Details::residual(*mtx_, toTpetra(X), toTpetra(B), toTpetra(R));
}

////////////////////////////////////////////
////////////////////////////////////////////
// End of TpetrCrsMatrix class definition //
////////////////////////////////////////////
////////////////////////////////////////////

}  // namespace Xpetra

#endif  // XPETRA_TPETRACRSMATRIX_DEF_HPP
