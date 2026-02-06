// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_TPETRABLOCKCRSMATRIX_DEF_HPP
#define XPETRA_TPETRABLOCKCRSMATRIX_DEF_HPP

#include "Xpetra_TpetraBlockCrsMatrix_decl.hpp"
#include "Xpetra_TpetraCrsGraph.hpp"

namespace Xpetra {

//! Constructor specifying fixed number of entries for each row (not implemented)
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap,
                         size_t maxNumEntriesPerRow,
                         const Teuchos::RCP<Teuchos::ParameterList> &params) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in" + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Constructor specifying (possibly different) number of entries in each row (not implemented)
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap,
                         const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc,
                         const Teuchos::RCP<Teuchos::ParameterList> &params) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in" + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Constructor specifying column Map and fixed number of entries for each row (not implemented)
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap,
                         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &colMap,
                         size_t maxNumEntriesPerRow,
                         const Teuchos::RCP<Teuchos::ParameterList> &params) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in" + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Constructor specifying column Map and number of entries in each row (not implemented)
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap,
                         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &colMap,
                         const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc,
                         const Teuchos::RCP<Teuchos::ParameterList> &params) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in" + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Constructor specifying a previously constructed graph ( not implemented )
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP<const CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > &graph,
                         const Teuchos::RCP<Teuchos::ParameterList> &params)
// : mtx_(Teuchos::rcp(new Tpetra::BlockCrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node >(toTpetra(graph), params)))
// * there is no Tpetra::BlockCrsMatrix(graph, params) c'tor.  We throw anyways here so no need to set mtx_.
{
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in" + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Constructor specifying a previously constructed graph & blocksize
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP<const CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > &graph,
                         const LocalOrdinal blockSize)
  : mtx_(Teuchos::rcp(new Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(*toTpetra(graph), blockSize))) {}

//! Constructor specifying a previously constructed graph, point maps & blocksize
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP<const CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > &graph,
                         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &pointDomainMap,
                         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &pointRangeMap,
                         const LocalOrdinal blockSize)
  : mtx_(Teuchos::rcp(new Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(*toTpetra(graph), *toTpetra(pointDomainMap), *toTpetra(pointRangeMap), blockSize))) {}

//! Constructor for a fused import ( not implemented )
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &sourceMatrix,
                         const Import<LocalOrdinal, GlobalOrdinal, Node> &importer,
                         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &domainMap,
                         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rangeMap,
                         const Teuchos::RCP<Teuchos::ParameterList> &params) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in" + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Constructor for a fused export (not implemented(
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &sourceMatrix,
                         const Export<LocalOrdinal, GlobalOrdinal, Node> &exporter,
                         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &domainMap,
                         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rangeMap,
                         const Teuchos::RCP<Teuchos::ParameterList> &params) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in" + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Constructor for a fused import ( not implemented )
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &sourceMatrix,
                         const Import<LocalOrdinal, GlobalOrdinal, Node> &RowImporter,
                         const Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> > DomainImporter,
                         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &domainMap,
                         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rangeMap,
                         const Teuchos::RCP<Teuchos::ParameterList> &params) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in" + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Constructor for a fused export (not implemented(
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &sourceMatrix,
                         const Export<LocalOrdinal, GlobalOrdinal, Node> &RowExporter,
                         const Teuchos::RCP<const Export<LocalOrdinal, GlobalOrdinal, Node> > DomainExporter,
                         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &domainMap,
                         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rangeMap,
                         const Teuchos::RCP<Teuchos::ParameterList> &params) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in" + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Destructor.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ~TpetraBlockCrsMatrix() {}

//@}

//! Insert matrix entries, using global TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::IDs (not implemented)
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    insertGlobalValues(GlobalOrdinal globalRow,
                       const ArrayView<const GlobalOrdinal> &cols,
                       const ArrayView<const Scalar> &vals) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in" + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Insert matrix entries, using local IDs (not implemented)
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    insertLocalValues(LocalOrdinal localRow,
                      const ArrayView<const LocalOrdinal> &cols,
                      const ArrayView<const Scalar> &vals) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in" + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Replace matrix entries, using global IDs (not implemented)
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceGlobalValues(GlobalOrdinal globalRow,
                        const ArrayView<const GlobalOrdinal> &cols,
                        const ArrayView<const Scalar> &vals) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in" + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Replace matrix entries, using local IDs.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &cols, const ArrayView<const Scalar> &vals) {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::replaceLocalValues");
  mtx_->replaceLocalValues(localRow, cols.getRawPtr(), vals.getRawPtr(), cols.size());
}

//! Set all matrix entries equal to scalarThis.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    setAllToScalar(const Scalar &alpha) {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::setAllToScalar");
  mtx_->setAllToScalar(alpha);
}

//! Scale the current values of a matrix, this = alpha*this (not implemented)
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    scale(const Scalar &alpha) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Allocates and returns ArrayRCPs of the Crs arrays --- This is an Xpetra-only routine.
//** \warning This is an expert-only routine and should not be called from user code. (not implemented)
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    allocateAllValues(size_t numNonZeros, ArrayRCP<size_t> &rowptr, ArrayRCP<LocalOrdinal> &colind, ArrayRCP<Scalar> &values) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Sets the 1D pointer arrays of the graph (not impelmented)
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    setAllValues(const ArrayRCP<size_t> &rowptr, const ArrayRCP<LocalOrdinal> &colind, const ArrayRCP<Scalar> &values) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Gets the 1D pointer arrays of the graph (not implemented)
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getAllValues(ArrayRCP<const size_t> &rowptr,
                 ArrayRCP<const LocalOrdinal> &colind,
                 ArrayRCP<const Scalar> &values) const {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Gets the 1D pointer arrays of the graph (not implemented)
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getAllValues(ArrayRCP<Scalar> &values) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//@}

// Transformational Methods
//@{

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    resumeFill(const RCP<ParameterList> &params) {
  /*noop*/
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    fillComplete(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &domainMap,
                 const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rangeMap,
                 const RCP<ParameterList> &params) {
  /*noop*/
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    fillComplete(const RCP<ParameterList> &params) {
  /*noop*/
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceDomainMapAndImporter(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &newDomainMap,
                                Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> > &newImporter) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    expertStaticFillComplete(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &domainMap,
                             const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rangeMap,
                             const RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> > &importer,
                             const RCP<const Export<LocalOrdinal, GlobalOrdinal, Node> > &exporter,
                             const RCP<ParameterList> &params) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//@}

//! @name Methods implementing RowMatrix

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getRowMap() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getRowMap");
  return toXpetra(mtx_->getRowMap());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getColMap() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getColMap");
  return toXpetra(mtx_->getColMap());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getCrsGraph() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getCrsGraph");
  using G_t              = Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  using G_x              = TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  RCP<G_t> t_graph       = Teuchos::rcp_const_cast<G_t>(Teuchos::rcpFromRef(mtx_->getCrsGraph()));
  RCP<const G_x> x_graph = rcp(new G_x(t_graph));
  return x_graph;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getGlobalNumRows() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getGlobalNumRows");
  return mtx_->getGlobalNumRows();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getGlobalNumCols() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getGlobalNumCols");
  return mtx_->getGlobalNumCols();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalNumRows() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getLocalNumRows");
  return mtx_->getLocalNumRows();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalNumCols() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getLocalNumCols");
  return mtx_->getLocalNumCols();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getGlobalNumEntries() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getGlobalNumEntries");
  return mtx_->getGlobalNumEntries();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalNumEntries() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getLocalNumEntries");
  return mtx_->getLocalNumEntries();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getNumEntriesInLocalRow(LocalOrdinal localRow) const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getNumEntriesInLocalRow");
  return mtx_->getNumEntriesInLocalRow(localRow);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getNumEntriesInGlobalRow");
  return mtx_->getNumEntriesInGlobalRow(globalRow);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalMaxNumRowEntries() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getGlobalMaxNumRowEntries");
  return mtx_->getGlobalMaxNumRowEntries();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalMaxNumRowEntries() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getLocalMaxNumRowEntries");
  return mtx_->getLocalMaxNumRowEntries();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isLocallyIndexed() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::isLocallyIndexed");
  return mtx_->isLocallyIndexed();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isGloballyIndexed() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::isGloballyIndexed");
  return mtx_->isGloballyIndexed();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isFillComplete() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::isFillComplete");
  return mtx_->isFillComplete();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isFillActive() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::isFillActive");
  return false;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename ScalarTraits<Scalar>::magnitudeType TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getFrobeniusNorm() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getFrobeniusNorm");
  return mtx_->getFrobeniusNorm();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::supportsRowViews() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::supportsRowViews");
  return mtx_->supportsRowViews();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalRowCopy(LocalOrdinal LocalRow,
                    const ArrayView<LocalOrdinal> &Indices,
                    const ArrayView<Scalar> &Values,
                    size_t &NumEntries) const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getLocalRowCopy");
  typename Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::nonconst_local_inds_host_view_type indices("indices", Indices.size());
  typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::nonconst_values_host_view_type values("values", Values.size());

  mtx_->getLocalRowCopy(LocalRow, indices, values, NumEntries);
  for (size_t i = 0; i < NumEntries; ++i) {
    Indices[i] = indices(i);
    Values[i]  = values(i);
  }
}
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &Indices,
                    ArrayView<const Scalar> &Values) const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getLocalRowView");
  typename Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_inds_host_view_type indices;
  typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::values_host_view_type values;

  mtx_->getLocalRowView(LocalRow, indices, values);
  Indices = ArrayView<const LocalOrdinal>(indices.data(), indices.extent(0));
  Values  = ArrayView<const Scalar>(reinterpret_cast<const Scalar *>(values.data()), values.extent(0));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getGlobalRowView(GlobalOrdinal GlobalRow,
                     ArrayView<const GlobalOrdinal> &Indices,
                     ArrayView<const Scalar> &Values) const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getGlobalRowView");
  typename Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::global_inds_host_view_type indices;
  typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::values_host_view_type values;

  mtx_->getGlobalRowView(GlobalRow, indices, values);
  Indices = ArrayView<const GlobalOrdinal>(indices.data(), indices.extent(0));
  Values  = ArrayView<const Scalar>(reinterpret_cast<const Scalar *>(values.data()), values.extent(0));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getGlobalRowCopy(GlobalOrdinal GlobalRow,
                     const ArrayView<GlobalOrdinal> &Indices,
                     const ArrayView<Scalar> &Values,
                     size_t &NumEntries) const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getGlobalRowCopy");
  typename Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::nonconst_global_inds_host_view_type indices("indices", Indices.size());
  typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::nonconst_values_host_view_type values("values", Values.size());

  mtx_->getGlobalRowCopy(GlobalRow, indices, values, NumEntries);
  for (size_t i = 0; i < NumEntries; ++i) {
    Indices[i] = indices(i);
    Values[i]  = values(i);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    haveGlobalConstants() const { return true; }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    apply(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &X,
          MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &Y,
          Teuchos::ETransp mode,
          Scalar alpha,
          Scalar beta) const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::apply");
  mtx_->apply(toTpetra(X), toTpetra(Y), mode, alpha, beta);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    apply(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &X,
          MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &Y,
          Teuchos::ETransp mode,
          Scalar alpha,
          Scalar beta,
          bool sumInterfaceValues,
          const RCP<Import<LocalOrdinal, GlobalOrdinal, Node> > &regionInterfaceImporter,
          const Teuchos::ArrayRCP<LocalOrdinal> &regionInterfaceLIDs) const {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getDomainMap() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getDomainMap");
  return toXpetra(mtx_->getDomainMap());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getRangeMap() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getRangeMap");
  return toXpetra(mtx_->getRangeMap());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    description() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::description");
  return mtx_->description();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    describe(Teuchos::FancyOStream &out,
             const Teuchos::EVerbosityLevel verbLevel) const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::describe");
  mtx_->describe(out, verbLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    setObjectLabel(const std::string &objectLabel) {
  XPETRA_MONITOR("TpetraCrsMatrix::setObjectLabel");
  Teuchos::LabeledObject::setObjectLabel(objectLabel);
  mtx_->setObjectLabel(objectLabel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalDiagCopy(Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &diag) const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getLocalDiagCopy");
  XPETRA_DYNAMIC_CAST(TpetraVectorClass,
                      diag,
                      tDiag,
                      "Xpetra::TpetraBlockCrsMatrix.getLocalDiagCopy() only accept Xpetra::TpetraVector as input arguments.");
  mtx_->getLocalDiagCopy(*tDiag.getTpetra_Vector());
}

//! Get a copy of the diagonal entries owned by this node, with local row indices.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalDiagCopy(Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &diag,
                     const Teuchos::ArrayView<const size_t> &offsets) const {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Get a copy of the diagonal entries owned by this node, with local row indices.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalDiagCopy(Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &diag,
                     const Kokkos::View<const size_t *, typename Node::device_type, Kokkos::MemoryUnmanaged> &offsets) const {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalDiagOffsets(Teuchos::ArrayRCP<size_t> &offsets) const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getLocalDiagOffsets");

  const size_t lclNumRows = mtx_->getGraph()->getLocalNumRows();
  if (static_cast<size_t>(offsets.size()) < lclNumRows) {
    offsets.resize(lclNumRows);
  }

  // The input ArrayRCP must always be a host pointer.  Thus, if
  // device_type::memory_space is Kokkos::HostSpace, it's OK for us
  // to write to that allocation directly as a Kokkos::View.
  typedef typename Node::device_type device_type;
  typedef typename device_type::memory_space memory_space;
  if (std::is_same<memory_space, Kokkos::HostSpace>::value) {
    // It is always syntactically correct to assign a raw host
    // pointer to a device View, so this code will compile correctly
    // even if this branch never runs.
    typedef Kokkos::View<size_t *, device_type, Kokkos::MemoryUnmanaged> output_type;
    output_type offsetsOut(offsets.getRawPtr(), offsets.size());
    mtx_->getLocalDiagOffsets(offsetsOut);
  } else {
    Kokkos::View<size_t *, device_type> offsetsTmp("diagOffsets", offsets.size());
    mtx_->getLocalDiagOffsets(offsetsTmp);
    typedef Kokkos::View<size_t *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> output_type;
    output_type offsetsOut(offsets.getRawPtr(), offsets.size());
    Kokkos::deep_copy(offsetsOut, offsetsTmp);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceDiag(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &diag) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix::replaceDiag: function not implemented in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    leftScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &x) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    rightScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &x) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMap() const {
  XPETRA_MONITOR("TpetraBlockCrsMatrix::getMap");
  return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(mtx_->getMap()));
}

//! Import.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doImport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &source,
             const Import<LocalOrdinal, GlobalOrdinal, Node> &importer, CombineMode CM) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Export.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doExport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &dest,
             const Import<LocalOrdinal, GlobalOrdinal, Node> &importer, CombineMode CM) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Import (using an Exporter).
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doImport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &source,
             const Export<LocalOrdinal, GlobalOrdinal, Node> &exporter, CombineMode CM) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

//! Export (using an Importer).
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doExport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &dest,
             const Export<LocalOrdinal, GlobalOrdinal, Node> &exporter, CombineMode CM) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    removeEmptyProcessesInPlace(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &newMap) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    hasMatrix() const {
  return !mtx_.is_null();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP<Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &mtx)
  : mtx_(mtx) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getTpetra_BlockCrsMatrix() const {
  return mtx_;
}

// TODO: remove
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getTpetra_BlockCrsMatrixNonConst() const {
  return mtx_;
}

// was:     typedef typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type local_matrix_type;
// using local_matrix_type = typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalMatrixDevice() const {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix does not support getLocalMatrix due to missing Kokkos::CrsMatrix in Tpetra's experimental implementation in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));

#ifndef __NVCC__
  local_matrix_type ret;
#endif  // __NVCC__

  TEUCHOS_UNREACHABLE_RETURN(ret);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
#if KOKKOS_VERSION >= 40799
typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type::host_mirror_type
#else
typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type::HostMirror
#endif
TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalMatrixHost() const {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix does not support getLocalMatrix due to missing Kokkos::CrsMatrix in Tpetra's experimental implementation in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));

#ifndef __NVCC__
#if KOKKOS_VERSION >= 40799
  typename local_matrix_type::host_mirror_type ret;
#else
  typename local_matrix_type::HostMirror ret;
#endif
#endif  // __NVCC__

  TEUCHOS_UNREACHABLE_RETURN(ret);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    setAllValues(const typename local_matrix_type::row_map_type &ptr,
                 const typename local_matrix_type::StaticCrsGraphType::entries_type::non_const_type &ind,
                 const typename local_matrix_type::values_type &val) {
  throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix does not support setAllValues due to missing Kokkos::CrsMatrix in Tpetra's experimental implementation in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

}  // namespace Xpetra

#endif  // XPETRA_TPETRABLOCKCRSMATRIX_DEF_HPP
