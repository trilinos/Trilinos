// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_TPETRACRSGRAPH_DEF_HPP
#define XPETRA_TPETRACRSGRAPH_DEF_HPP
#include "Xpetra_TpetraConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Tpetra_CrsGraph.hpp"

#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_TpetraCrsGraph_decl.hpp"
#include "Xpetra_Utils.hpp"
#include "Xpetra_TpetraMap.hpp"
#include "Xpetra_TpetraImport.hpp"
#include "Xpetra_TpetraExport.hpp"

namespace Xpetra {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsGraph(const RCP<const map_type> &rowMap, size_t maxNumEntriesPerRow, const RCP<ParameterList> &params)
  : graph_(Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rowMap), maxNumEntriesPerRow, params))) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsGraph(const RCP<const Map> &rowMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, const RCP<ParameterList> &params)
  : graph_(Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rowMap), NumEntriesPerRowToAlloc(), params))) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsGraph(const RCP<const Map> &rowMap, const RCP<const Map> &colMap, size_t maxNumEntriesPerRow, const RCP<ParameterList> &params)
  : graph_(Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rowMap), toTpetra(colMap), maxNumEntriesPerRow, params))) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsGraph(const RCP<const Map> &rowMap, const RCP<const Map> &colMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, const RCP<ParameterList> &params)
  : graph_(Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rowMap), toTpetra(colMap), NumEntriesPerRowToAlloc(), params))) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
    TpetraCrsGraph(const Teuchos::RCP<const CrsGraph> &sourceGraph,
                   const Import &importer,
                   const Teuchos::RCP<const Map> &domainMap,
                   const Teuchos::RCP<const Map> &rangeMap,
                   const Teuchos::RCP<Teuchos::ParameterList> &params) {
  typedef Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> MyTpetraCrsGraph;
  XPETRA_DYNAMIC_CAST(const TpetraCrsGraphClass, *sourceGraph, tSourceGraph, "Xpetra::TpetraCrsMatrix constructor only accepts Xpetra::TpetraCrsMatrix as the input argument.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > v = tSourceGraph.getTpetra_CrsGraph();

  RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > myDomainMap = domainMap != Teuchos::null ? toTpetra(domainMap) : Teuchos::null;
  RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > myRangeMap  = rangeMap != Teuchos::null ? toTpetra(rangeMap) : Teuchos::null;
  graph_                                                                 = Tpetra::importAndFillCompleteCrsGraph<MyTpetraCrsGraph>(v, toTpetra(importer), myDomainMap, myRangeMap, params);
  bool restrictComm                                                      = false;
  if (!params.is_null()) restrictComm = params->get("Restrict Communicator", restrictComm);
  if (restrictComm && graph_->getRowMap().is_null()) graph_ = Teuchos::null;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
    TpetraCrsGraph(const Teuchos::RCP<const Map> &rowMap,
                   const Teuchos::RCP<const Map> &colMap,
                   const typename local_graph_type::row_map_type &rowPointers,
                   const typename local_graph_type::entries_type::non_const_type &columnIndices,
                   const Teuchos::RCP<Teuchos::ParameterList> &plist)
  : graph_(Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rowMap), toTpetra(colMap), rowPointers, columnIndices, plist))) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
    TpetraCrsGraph(const Teuchos::RCP<const map_type> &rowMap,
                   const Teuchos::RCP<const map_type> &colMap,
                   const local_graph_type &lclGraph,
                   const Teuchos::RCP<Teuchos::ParameterList> &params)
  : graph_(Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rowMap), toTpetra(colMap), lclGraph, params))) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
    TpetraCrsGraph(const local_graph_type &lclGraph,
                   const Teuchos::RCP<const map_type> &rowMap,
                   const Teuchos::RCP<const map_type> &colMap,
                   const Teuchos::RCP<const map_type> &domainMap,
                   const Teuchos::RCP<const map_type> &rangeMap,
                   const Teuchos::RCP<Teuchos::ParameterList> &params)
  : graph_(Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>(lclGraph, toTpetra(rowMap), toTpetra(colMap), toTpetra(domainMap), toTpetra(rangeMap), params))) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
    TpetraCrsGraph(const Teuchos::RCP<const Map> &rowMap,
                   const Teuchos::RCP<const Map> &colMap,
                   const Teuchos::ArrayRCP<size_t> &rowPointers,
                   const Teuchos::ArrayRCP<LocalOrdinal> &columnIndices,
                   const Teuchos::RCP<Teuchos::ParameterList> &params)
  : graph_(Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rowMap), toTpetra(colMap), rowPointers, columnIndices, params))) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::~TpetraCrsGraph() {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::insertGlobalIndices(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &indices) {
  XPETRA_MONITOR("TpetraCrsGraph::insertGlobalIndices");
  graph_->insertGlobalIndices(globalRow, indices);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::insertLocalIndices(const LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &indices) {
  XPETRA_MONITOR("TpetraCrsGraph::insertLocalIndices");
  graph_->insertLocalIndices(localRow, indices);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::removeLocalIndices(LocalOrdinal localRow) {
  XPETRA_MONITOR("TpetraCrsGraph::removeLocalIndices");
  graph_->removeLocalIndices(localRow);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
    allocateAllIndices(size_t numNonZeros, ArrayRCP<size_t> &rowptr, ArrayRCP<LocalOrdinal> &colind) {
  rowptr.resize(getLocalNumRows() + 1);
  colind.resize(numNonZeros);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
    setAllIndices(const ArrayRCP<size_t> &rowptr, const ArrayRCP<LocalOrdinal> &colind) {
  graph_->setAllIndices(rowptr, colind);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
    getAllIndices(ArrayRCP<const size_t> &rowptr, ArrayRCP<const LocalOrdinal> &colind) const {
  rowptr = Kokkos::Compat::persistingView(graph_->getLocalRowPtrsHost());
  colind = Kokkos::Compat::persistingView(graph_->getLocalIndicesHost());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::fillComplete(const RCP<const Map> &domainMap, const RCP<const Map> &rangeMap, const RCP<ParameterList> &params) {
  XPETRA_MONITOR("TpetraCrsGraph::fillComplete");
  graph_->fillComplete(toTpetra(domainMap), toTpetra(rangeMap), params);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::fillComplete(const RCP<ParameterList> &params) {
  XPETRA_MONITOR("TpetraCrsGraph::fillComplete");
  graph_->fillComplete(params);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
    expertStaticFillComplete(const Teuchos::RCP<const map_type> &domainMap,
                             const Teuchos::RCP<const map_type> &rangeMap,
                             const Teuchos::RCP<const Import> &importer,
                             const Teuchos::RCP<const Export> &exporter,
                             const Teuchos::RCP<Teuchos::ParameterList> &params) {
  XPETRA_MONITOR("TpetraCrsGraph::expertStaticFillComplete");
  RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > myImport;
  RCP<const Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node> > myExport;

  if (importer != Teuchos::null) {
    XPETRA_DYNAMIC_CAST(const TpetraImportClass, *importer, tImporter, "Xpetra::TpetraCrsGraph::expertStaticFillComplete only accepts Xpetra::TpetraImport.");
    myImport = tImporter.getTpetra_Import();
  }
  if (exporter != Teuchos::null) {
    XPETRA_DYNAMIC_CAST(const TpetraExportClass, *exporter, tExporter, "Xpetra::TpetraCrsGraph::expertStaticFillComplete only accepts Xpetra::TpetraExport.");
    myExport = tExporter.getTpetra_Export();
  }

  graph_->expertStaticFillComplete(toTpetra(domainMap), toTpetra(rangeMap), myImport, myExport, params);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Comm<int> > TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getComm() const {
  XPETRA_MONITOR("TpetraCrsGraph::getComm");
  return graph_->getComm();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getRowMap() const {
  XPETRA_MONITOR("TpetraCrsGraph::getRowMap");
  return toXpetra(graph_->getRowMap());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getColMap() const {
  XPETRA_MONITOR("TpetraCrsGraph::getColMap");
  return toXpetra(graph_->getColMap());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getDomainMap() const {
  XPETRA_MONITOR("TpetraCrsGraph::getDomainMap");
  return toXpetra(graph_->getDomainMap());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getRangeMap() const {
  XPETRA_MONITOR("TpetraCrsGraph::getRangeMap");
  return toXpetra(graph_->getRangeMap());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> > TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getImporter() const {
  XPETRA_MONITOR("TpetraCrsGraph::getImporter");
  return toXpetra(graph_->getImporter());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Export<LocalOrdinal, GlobalOrdinal, Node> > TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getExporter() const {
  XPETRA_MONITOR("TpetraCrsGraph::getExporter");
  return toXpetra(graph_->getExporter());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getGlobalNumRows() const {
  XPETRA_MONITOR("TpetraCrsGraph::getGlobalNumRows");
  return graph_->getGlobalNumRows();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getGlobalNumCols() const {
  XPETRA_MONITOR("TpetraCrsGraph::getGlobalNumCols");
  return graph_->getGlobalNumCols();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getLocalNumRows() const {
  XPETRA_MONITOR("TpetraCrsGraph::getLocalNumRows");
  return graph_->getLocalNumRows();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getLocalNumCols() const {
  XPETRA_MONITOR("TpetraCrsGraph::getLocalNumCols");
  return graph_->getLocalNumCols();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getIndexBase() const {
  XPETRA_MONITOR("TpetraCrsGraph::getIndexBase");
  return graph_->getIndexBase();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getGlobalNumEntries() const {
  XPETRA_MONITOR("TpetraCrsGraph::getGlobalNumEntries");
  return graph_->getGlobalNumEntries();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getLocalNumEntries() const {
  XPETRA_MONITOR("TpetraCrsGraph::getLocalNumEntries");
  return graph_->getLocalNumEntries();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const {
  XPETRA_MONITOR("TpetraCrsGraph::getNumEntriesInGlobalRow");
  return graph_->getNumEntriesInGlobalRow(globalRow);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getNumEntriesInLocalRow(LocalOrdinal localRow) const {
  XPETRA_MONITOR("TpetraCrsGraph::getNumEntriesInLocalRow");
  return graph_->getNumEntriesInLocalRow(localRow);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getNumAllocatedEntriesInGlobalRow(GlobalOrdinal globalRow) const {
  XPETRA_MONITOR("TpetraCrsGraph::getNumAllocatedEntriesInGlobalRow");
  return graph_->getNumAllocatedEntriesInGlobalRow(globalRow);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getNumAllocatedEntriesInLocalRow(LocalOrdinal localRow) const {
  XPETRA_MONITOR("TpetraCrsGraph::getNumAllocatedEntriesInLocalRow");
  return graph_->getNumAllocatedEntriesInLocalRow(localRow);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getGlobalMaxNumRowEntries() const {
  XPETRA_MONITOR("TpetraCrsGraph::getGlobalMaxNumRowEntries");
  return graph_->getGlobalMaxNumRowEntries();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getLocalMaxNumRowEntries() const {
  XPETRA_MONITOR("TpetraCrsGraph::getLocalMaxNumRowEntries");
  return graph_->getLocalMaxNumRowEntries();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::hasColMap() const {
  XPETRA_MONITOR("TpetraCrsGraph::hasColMap");
  return graph_->hasColMap();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::isLocallyIndexed() const {
  XPETRA_MONITOR("TpetraCrsGraph::isLocallyIndexed");
  return graph_->isLocallyIndexed();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::isGloballyIndexed() const {
  XPETRA_MONITOR("TpetraCrsGraph::isGloballyIndexed");
  return graph_->isGloballyIndexed();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::isFillComplete() const {
  XPETRA_MONITOR("TpetraCrsGraph::isFillComplete");
  return graph_->isFillComplete();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::isStorageOptimized() const {
  XPETRA_MONITOR("TpetraCrsGraph::isStorageOptimized");
  return graph_->isStorageOptimized();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &Indices) const {
  XPETRA_MONITOR("TpetraCrsGraph::getGlobalRowView");
  typename Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::global_inds_host_view_type indices;
  graph_->getGlobalRowView(GlobalRow, indices);
  Indices = ArrayView<const GlobalOrdinal>(indices.data(), indices.extent(0));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &Indices) const {
  XPETRA_MONITOR("TpetraCrsGraph::getLocalRowView");
  typename Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_inds_host_view_type indices;
  graph_->getLocalRowView(LocalRow, indices);
  Indices = ArrayView<const LocalOrdinal>(indices.data(), indices.extent(0));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
#if KOKKOS_VERSION >= 40799
typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type::host_mirror_type TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getLocalGraphHost() const {
#else
typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type::HostMirror TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getLocalGraphHost() const {
#endif
  return getTpetra_CrsGraph()->getLocalGraphHost();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getLocalGraphDevice() const {
  return getTpetra_CrsGraph()->getLocalGraphDevice();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getLocalDiagOffsets(const Kokkos::View<size_t *, typename Node::device_type, Kokkos::MemoryUnmanaged> &offsets) const {
  getTpetra_CrsGraph()->getLocalDiagOffsets(offsets);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::computeGlobalConstants() {
  // mfh 07 May 2018: See GitHub Issue #2565.
  graph_->computeGlobalConstants();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
std::string TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::description() const {
  XPETRA_MONITOR("TpetraCrsGraph::description");
  return graph_->description();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
  XPETRA_MONITOR("TpetraCrsGraph::describe");
  graph_->describe(out, verbLevel);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<const size_t> TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getNodeRowPtrs() const {
  XPETRA_MONITOR("TpetraCrsGraph::getNodeRowPtrs");
  return Kokkos::Compat::persistingView(graph_->getLocalRowPtrsHost());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getMap() const {
  XPETRA_MONITOR("TpetraCrsGraph::getMap");
  return rcp(new TpetraMap(graph_->getMap()));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::doImport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &source,
                                                                 const Import &importer, CombineMode CM) {
  XPETRA_MONITOR("TpetraCrsGraph::doImport");

  XPETRA_DYNAMIC_CAST(const TpetraCrsGraphClass, source, tSource, "Xpetra::TpetraCrsGraph::doImport only accept Xpetra::TpetraCrsGraph as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > v = tSource.getTpetra_CrsGraph();
  // graph_->doImport(toTpetraCrsGraph(source), *tImporter.getTpetra_Import(), toTpetra(CM));

  graph_->doImport(*v, toTpetra(importer), toTpetra(CM));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::doExport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &dest,
                                                                 const Import &importer, CombineMode CM) {
  XPETRA_MONITOR("TpetraCrsGraph::doExport");

  XPETRA_DYNAMIC_CAST(const TpetraCrsGraphClass, dest, tDest, "Xpetra::TpetraCrsGraph::doImport only accept Xpetra::TpetraCrsGraph as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > v = tDest.getTpetra_CrsGraph();
  graph_->doExport(*v, toTpetra(importer), toTpetra(CM));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::doImport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &source,
                                                                 const Export &exporter, CombineMode CM) {
  XPETRA_MONITOR("TpetraCrsGraph::doImport");

  XPETRA_DYNAMIC_CAST(const TpetraCrsGraphClass, source, tSource, "Xpetra::TpetraCrsGraph::doImport only accept Xpetra::TpetraCrsGraph as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > v = tSource.getTpetra_CrsGraph();

  graph_->doImport(*v, toTpetra(exporter), toTpetra(CM));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::doExport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &dest,
                                                                 const Export &exporter, CombineMode CM) {
  XPETRA_MONITOR("TpetraCrsGraph::doExport");

  XPETRA_DYNAMIC_CAST(const TpetraCrsGraphClass, dest, tDest, "Xpetra::TpetraCrsGraph::doImport only accept Xpetra::TpetraCrsGraph as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > v = tDest.getTpetra_CrsGraph();

  graph_->doExport(*v, toTpetra(exporter), toTpetra(CM));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::TpetraCrsGraph(const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > &graph)
  : graph_(graph) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getTpetra_CrsGraph() const { return graph_; }

}  // namespace Xpetra
#endif  // XPETRA_TPETRACRSGRAPH_DEF_HPP
