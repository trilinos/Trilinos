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


template<class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::TpetraCrsGraph(const RCP< const map_type > &rowMap, size_t maxNumEntriesPerRow, const RCP< ParameterList > &params)
: graph_(Teuchos::rcp(new Tpetra::CrsGraph< LocalOrdinal, GlobalOrdinal, Node >(toTpetra(rowMap), maxNumEntriesPerRow, params))) {  }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::TpetraCrsGraph(const RCP< const Map > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, const RCP< ParameterList > &params)
: graph_(Teuchos::rcp(new Tpetra::CrsGraph< LocalOrdinal, GlobalOrdinal, Node >(toTpetra(rowMap), NumEntriesPerRowToAlloc(), params))) {  }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::TpetraCrsGraph(const RCP< const Map>&rowMap, const RCP< const Map > &colMap, size_t maxNumEntriesPerRow, const RCP< ParameterList > &params)
: graph_(Teuchos::rcp(new Tpetra::CrsGraph< LocalOrdinal, GlobalOrdinal, Node >(toTpetra(rowMap), toTpetra(colMap), maxNumEntriesPerRow, params))) {  }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::TpetraCrsGraph(const RCP< const Map > &rowMap, const RCP< const Map > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, const RCP< ParameterList > &params)
: graph_(Teuchos::rcp(new Tpetra::CrsGraph< LocalOrdinal, GlobalOrdinal, Node >(toTpetra(rowMap), toTpetra(colMap), NumEntriesPerRowToAlloc(), params))) {  }


template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::
TpetraCrsGraph(const Teuchos::RCP<const CrsGraph >& sourceGraph,
               const Import & importer,
               const Teuchos::RCP<const Map >& domainMap,
               const Teuchos::RCP<const Map >& rangeMap,
               const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> MyTpetraCrsGraph;
  XPETRA_DYNAMIC_CAST(const TpetraCrsGraphClass, *sourceGraph, tSourceGraph, "Xpetra::TpetraCrsMatrix constructor only accepts Xpetra::TpetraCrsMatrix as the input argument.");//TODO: remove and use toTpetra()
  RCP< const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > v = tSourceGraph.getTpetra_CrsGraph();

  RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > myDomainMap = domainMap!=Teuchos::null ? toTpetra(domainMap) : Teuchos::null;
  RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > myRangeMap  = rangeMap!=Teuchos::null  ? toTpetra(rangeMap)  : Teuchos::null;
  graph_=Tpetra::importAndFillCompleteCrsGraph<MyTpetraCrsGraph>(v,toTpetra(importer),myDomainMap,myRangeMap,params);
  bool restrictComm=false;
  if(!params.is_null()) restrictComm = params->get("Restrict Communicator",restrictComm);
  if(restrictComm && graph_->getRowMap().is_null()) graph_=Teuchos::null;
  
}


template<class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::
TpetraCrsGraph(const Teuchos::RCP< const Map > &rowMap,
               const Teuchos::RCP< const Map > &colMap,
               const typename local_graph_type::row_map_type& rowPointers,
               const typename local_graph_type::entries_type::non_const_type& columnIndices,
               const Teuchos::RCP< Teuchos::ParameterList > &plist)
  : graph_(Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rowMap), toTpetra(colMap), rowPointers, columnIndices, plist))) {  }
  
  
template<class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::
TpetraCrsGraph(const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap,
               const local_graph_type& lclGraph,
               const Teuchos::RCP<Teuchos::ParameterList>& params)
  : graph_(Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rowMap), toTpetra(colMap), lclGraph, params))) {  }
  
template<class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::
TpetraCrsGraph(const local_graph_type& lclGraph,
               const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap,
               const Teuchos::RCP<const map_type>& domainMap,
               const Teuchos::RCP<const map_type>& rangeMap,
               const Teuchos::RCP<Teuchos::ParameterList>& params)
  : graph_(Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>(lclGraph, toTpetra(rowMap), toTpetra(colMap), toTpetra(domainMap), toTpetra(rangeMap), params))) {  }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::
TpetraCrsGraph(const Teuchos::RCP< const Map > &rowMap,
               const Teuchos::RCP< const Map > &colMap,
               const Teuchos::ArrayRCP<size_t>& rowPointers,
               const Teuchos::ArrayRCP<LocalOrdinal>& columnIndices,
               const Teuchos::RCP<Teuchos::ParameterList>& params)
  : graph_(Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rowMap), toTpetra(colMap), rowPointers, columnIndices, params))) {  }


template<class LocalOrdinal, class GlobalOrdinal, class Node>
 TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::~TpetraCrsGraph() {  }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::insertGlobalIndices(GlobalOrdinal globalRow, const ArrayView< const GlobalOrdinal > &indices)
{ XPETRA_MONITOR("TpetraCrsGraph::insertGlobalIndices"); graph_->insertGlobalIndices(globalRow, indices); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::insertLocalIndices(const LocalOrdinal localRow, const ArrayView< const LocalOrdinal > &indices)
{ XPETRA_MONITOR("TpetraCrsGraph::insertLocalIndices"); graph_->insertLocalIndices(localRow, indices); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::removeLocalIndices(LocalOrdinal localRow)
{ XPETRA_MONITOR("TpetraCrsGraph::removeLocalIndices"); graph_->removeLocalIndices(localRow); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::
allocateAllIndices(size_t numNonZeros,ArrayRCP<size_t> & rowptr, ArrayRCP<LocalOrdinal> & colind) {
  rowptr.resize(getLocalNumRows()+1); colind.resize(numNonZeros);
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::
setAllIndices(const ArrayRCP<size_t> & rowptr, const ArrayRCP<LocalOrdinal> & colind) {
  graph_->setAllIndices(rowptr,colind);
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::
getAllIndices(ArrayRCP<const size_t>& rowptr, ArrayRCP<const LocalOrdinal>& colind) const {
  rowptr = Kokkos::Compat::persistingView(graph_->getLocalRowPtrsHost());
  colind = Kokkos::Compat::persistingView(graph_->getLocalIndicesHost());
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::fillComplete(const RCP< const Map > &domainMap, const RCP< const Map > &rangeMap, const RCP< ParameterList > &params)
{ XPETRA_MONITOR("TpetraCrsGraph::fillComplete"); graph_->fillComplete(toTpetra(domainMap), toTpetra(rangeMap), params); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::fillComplete(const RCP< ParameterList > &params)
{ XPETRA_MONITOR("TpetraCrsGraph::fillComplete"); graph_->fillComplete(params); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::
expertStaticFillComplete (const Teuchos::RCP<const map_type>& domainMap,
                          const Teuchos::RCP<const map_type>& rangeMap,
                          const Teuchos::RCP<const Import>& importer,
                          const Teuchos::RCP<const Export>& exporter,                          
                          const Teuchos::RCP<Teuchos::ParameterList>& params) {
  XPETRA_MONITOR("TpetraCrsGraph::expertStaticFillComplete");
  RCP<const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > myImport;
  RCP<const Tpetra::Export<LocalOrdinal,GlobalOrdinal,Node> > myExport;
  
  if(importer!=Teuchos::null) {
    XPETRA_DYNAMIC_CAST( const TpetraImportClass , *importer, tImporter, "Xpetra::TpetraCrsGraph::expertStaticFillComplete only accepts Xpetra::TpetraImport.");
    myImport = tImporter.getTpetra_Import();
  }
  if(exporter!=Teuchos::null) {
    XPETRA_DYNAMIC_CAST( const TpetraExportClass , *exporter, tExporter, "Xpetra::TpetraCrsGraph::expertStaticFillComplete only accepts Xpetra::TpetraExport.");
    myExport = tExporter.getTpetra_Export();
  }
  
  graph_->expertStaticFillComplete(toTpetra(domainMap),toTpetra(rangeMap),myImport,myExport,params);
}


template<class LocalOrdinal, class GlobalOrdinal, class Node>
RCP< const Comm< int > > TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getComm() const
{ XPETRA_MONITOR("TpetraCrsGraph::getComm"); return graph_->getComm(); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
RCP< const Map<LocalOrdinal, GlobalOrdinal, Node> >  TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getRowMap() const
{ XPETRA_MONITOR("TpetraCrsGraph::getRowMap"); return toXpetra(graph_->getRowMap()); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
RCP< const Map<LocalOrdinal, GlobalOrdinal, Node> >  TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getColMap() const
{ XPETRA_MONITOR("TpetraCrsGraph::getColMap"); return toXpetra(graph_->getColMap()); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
RCP< const Map<LocalOrdinal, GlobalOrdinal, Node> >  TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getDomainMap() const
{ XPETRA_MONITOR("TpetraCrsGraph::getDomainMap"); return toXpetra(graph_->getDomainMap()); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
RCP< const Map<LocalOrdinal, GlobalOrdinal, Node> >  TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getRangeMap() const
{ XPETRA_MONITOR("TpetraCrsGraph::getRangeMap"); return toXpetra(graph_->getRangeMap()); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
RCP< const Import< LocalOrdinal, GlobalOrdinal, Node > > TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getImporter() const
{ XPETRA_MONITOR("TpetraCrsGraph::getImporter"); return toXpetra(graph_->getImporter()); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
RCP< const Export< LocalOrdinal, GlobalOrdinal, Node > > TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getExporter() const
{ XPETRA_MONITOR("TpetraCrsGraph::getExporter"); return toXpetra(graph_->getExporter()); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumRows() const
{ XPETRA_MONITOR("TpetraCrsGraph::getGlobalNumRows"); return graph_->getGlobalNumRows(); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumCols() const
{ XPETRA_MONITOR("TpetraCrsGraph::getGlobalNumCols"); return graph_->getGlobalNumCols(); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getLocalNumRows() const
{ XPETRA_MONITOR("TpetraCrsGraph::getLocalNumRows"); return graph_->getLocalNumRows(); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getLocalNumCols() const
{ XPETRA_MONITOR("TpetraCrsGraph::getLocalNumCols"); return graph_->getLocalNumCols(); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getIndexBase() const
{ XPETRA_MONITOR("TpetraCrsGraph::getIndexBase"); return graph_->getIndexBase(); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumEntries() const
{ XPETRA_MONITOR("TpetraCrsGraph::getGlobalNumEntries"); return graph_->getGlobalNumEntries(); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getLocalNumEntries() const
{ XPETRA_MONITOR("TpetraCrsGraph::getLocalNumEntries"); return graph_->getLocalNumEntries(); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const
{ XPETRA_MONITOR("TpetraCrsGraph::getNumEntriesInGlobalRow"); return graph_->getNumEntriesInGlobalRow(globalRow); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNumEntriesInLocalRow(LocalOrdinal localRow) const
{ XPETRA_MONITOR("TpetraCrsGraph::getNumEntriesInLocalRow"); return graph_->getNumEntriesInLocalRow(localRow); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNumAllocatedEntriesInGlobalRow(GlobalOrdinal globalRow) const
{ XPETRA_MONITOR("TpetraCrsGraph::getNumAllocatedEntriesInGlobalRow"); return graph_->getNumAllocatedEntriesInGlobalRow(globalRow); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNumAllocatedEntriesInLocalRow(LocalOrdinal localRow) const
{ XPETRA_MONITOR("TpetraCrsGraph::getNumAllocatedEntriesInLocalRow"); return graph_->getNumAllocatedEntriesInLocalRow(localRow); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalMaxNumRowEntries() const
{ XPETRA_MONITOR("TpetraCrsGraph::getGlobalMaxNumRowEntries"); return graph_->getGlobalMaxNumRowEntries(); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getLocalMaxNumRowEntries() const
{ XPETRA_MONITOR("TpetraCrsGraph::getLocalMaxNumRowEntries"); return graph_->getLocalMaxNumRowEntries(); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::hasColMap() const
{ XPETRA_MONITOR("TpetraCrsGraph::hasColMap"); return graph_->hasColMap(); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isLocallyIndexed() const
{ XPETRA_MONITOR("TpetraCrsGraph::isLocallyIndexed"); return graph_->isLocallyIndexed(); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isGloballyIndexed() const
{ XPETRA_MONITOR("TpetraCrsGraph::isGloballyIndexed"); return graph_->isGloballyIndexed(); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isFillComplete() const
{ XPETRA_MONITOR("TpetraCrsGraph::isFillComplete"); return graph_->isFillComplete(); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isStorageOptimized() const
{ XPETRA_MONITOR("TpetraCrsGraph::isStorageOptimized"); return graph_->isStorageOptimized(); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView< const GlobalOrdinal > &Indices) const
{
  XPETRA_MONITOR("TpetraCrsGraph::getGlobalRowView");
  typename Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::global_inds_host_view_type indices;
  graph_->getGlobalRowView(GlobalRow, indices);
  Indices = ArrayView<const GlobalOrdinal> (indices.data(), indices.extent(0));
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getLocalRowView(LocalOrdinal LocalRow, ArrayView< const LocalOrdinal > &Indices) const
{
  XPETRA_MONITOR("TpetraCrsGraph::getLocalRowView");
  typename Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::local_inds_host_view_type indices;
  graph_->getLocalRowView(LocalRow, indices);
  Indices = ArrayView<const LocalOrdinal> (indices.data(), indices.extent(0));
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
typename Xpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::local_graph_type::HostMirror TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getLocalGraphHost () const {
  return getTpetra_CrsGraph()->getLocalGraphHost();
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
typename Xpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::local_graph_type TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getLocalGraphDevice () const {
  return getTpetra_CrsGraph()->getLocalGraphDevice();
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::computeGlobalConstants() {
      // mfh 07 May 2018: See GitHub Issue #2565.
      graph_->computeGlobalConstants();
    }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
std::string TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::description() const
{ XPETRA_MONITOR("TpetraCrsGraph::description"); return graph_->description(); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
{ XPETRA_MONITOR("TpetraCrsGraph::describe"); graph_->describe(out, verbLevel); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP< const size_t > TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNodeRowPtrs() const
{ XPETRA_MONITOR("TpetraCrsGraph::getNodeRowPtrs"); return Kokkos::Compat::persistingView(graph_->getLocalRowPtrsHost()); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP< const Map<LocalOrdinal, GlobalOrdinal, Node> > TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getMap() const
{ XPETRA_MONITOR("TpetraCrsGraph::getMap"); return rcp( new TpetraMap(graph_->getMap()) ); }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::doImport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &source,
                                                               const Import &importer, CombineMode CM){
  XPETRA_MONITOR("TpetraCrsGraph::doImport");
  
  XPETRA_DYNAMIC_CAST(const TpetraCrsGraphClass, source, tSource, "Xpetra::TpetraCrsGraph::doImport only accept Xpetra::TpetraCrsGraph as input arguments.");//TODO: remove and use toTpetra()
  RCP< const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > v = tSource.getTpetra_CrsGraph();
  //graph_->doImport(toTpetraCrsGraph(source), *tImporter.getTpetra_Import(), toTpetra(CM));
  
  graph_->doImport(*v, toTpetra(importer), toTpetra(CM));
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::doExport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &dest,
                                                               const Import& importer, CombineMode CM) {
  XPETRA_MONITOR("TpetraCrsGraph::doExport");
  
  XPETRA_DYNAMIC_CAST(const TpetraCrsGraphClass, dest, tDest, "Xpetra::TpetraCrsGraph::doImport only accept Xpetra::TpetraCrsGraph as input arguments.");//TODO: remove and use toTpetra()
  RCP< const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > v = tDest.getTpetra_CrsGraph();
  graph_->doExport(*v, toTpetra(importer), toTpetra(CM));
  
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::doImport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &source,
                                                               const Export& exporter, CombineMode CM){
  XPETRA_MONITOR("TpetraCrsGraph::doImport");
  
  XPETRA_DYNAMIC_CAST(const TpetraCrsGraphClass, source, tSource, "Xpetra::TpetraCrsGraph::doImport only accept Xpetra::TpetraCrsGraph as input arguments.");//TODO: remove and use toTpetra()
  RCP< const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > v = tSource.getTpetra_CrsGraph();
  
  graph_->doImport(*v, toTpetra(exporter), toTpetra(CM));
  
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::doExport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &dest,
                                                               const Export& exporter, CombineMode CM) {
  XPETRA_MONITOR("TpetraCrsGraph::doExport");
  
  XPETRA_DYNAMIC_CAST(const TpetraCrsGraphClass, dest, tDest, "Xpetra::TpetraCrsGraph::doImport only accept Xpetra::TpetraCrsGraph as input arguments.");//TODO: remove and use toTpetra()
  RCP< const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > v = tDest.getTpetra_CrsGraph();
  
      graph_->doExport(*v, toTpetra(exporter), toTpetra(CM));
      
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::TpetraCrsGraph(const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > &graph) : graph_(graph)
{ }

template<class LocalOrdinal, class GlobalOrdinal, class Node>
RCP< const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getTpetra_CrsGraph() const
{ return graph_; }




#ifdef HAVE_XPETRA_EPETRA

#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))

  // specialization of TpetraCrsGraph for GO=LO=int
  template <>
  class TpetraCrsGraph<int,int,EpetraNode>
    : public CrsGraph<int,int,EpetraNode>
  {
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef EpetraNode Node;

    // The following typedef is used by the XPETRA_DYNAMIC_CAST() macro.
    typedef TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node> TpetraCrsGraphClass;
    typedef Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

  public:

    //! @name Constructor/Destructor Methods
    //@{

    //! Constructor specifying fixed number of entries for each row.
    TpetraCrsGraph(const RCP< const map_type > &rowMap, size_t maxNumEntriesPerRow, const RCP< ParameterList > &params=null) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "int", typeid(EpetraNode).name() );
    }

    //! Constructor specifying (possibly different) number of entries in each row.
    TpetraCrsGraph(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, const RCP< ParameterList > &params=null) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "int", typeid(EpetraNode).name() );
    }

    //! Constructor specifying column Map and fixed number of entries for each row.
    TpetraCrsGraph(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, const RCP< ParameterList > &params=null) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "int", typeid(EpetraNode).name() );
    }

    //! Constructor specifying column Map and number of entries in each row.
    TpetraCrsGraph(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, const RCP< ParameterList > &params=null) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "int", typeid(EpetraNode).name() );
    }

    /// \brief Constructor specifying column Map and arrays containing the graph in sorted, local ids.
    ///
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///
    /// \param rowPointers [in] The beginning of each row in the graph,
    ///   as in a CSR "rowptr" array.  The length of this vector should be
    ///   equal to the number of rows in the graph, plus one.  This last
    ///   entry should store the nunber of nonzeros in the graph.
    ///
    /// \param columnIndices [in] The local indices of the columns,
    ///   as in a CSR "colind" array.  The length of this vector
    ///   should be equal to the number of unknowns in the graph.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    TpetraCrsGraph(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap,
                   const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap,
                   const typename local_graph_type::row_map_type& rowPointers,
                   const typename local_graph_type::entries_type::non_const_type& columnIndices,
                   const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(),
                                   typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(),
                                   "int",
                                   typeid(EpetraNode).name());
    }

    /// \brief Constructor specifying column Map and a local (sorted)
    ///   graph, which the resulting CrsGraph views.
    ///
    /// Unlike most other CrsGraph constructors, successful completion
    /// of this constructor will result in a fill-complete graph.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///
    /// \param lclGraph [in] A locally indexed Kokkos::StaticCrsGraph
    ///   whose local row indices come from the specified row Map, and
    ///   whose local column indices come from the specified column
    ///   Map.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    TpetraCrsGraph(const Teuchos::RCP<const map_type>& rowMap,
                   const Teuchos::RCP<const map_type>& colMap,
                   const local_graph_type& lclGraph,
                   const Teuchos::RCP<Teuchos::ParameterList>& params) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(),
                                   typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(),
                                   "int",
                                   typeid(EpetraNode).name());
    }

    /// \brief Constructor specifying column, domain and range maps, and a
    ///   local (sorted) graph, which the resulting CrsGraph views.
    ///
    /// Unlike most other CrsGraph constructors, successful completion
    /// of this constructor will result in a fill-complete graph.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///
    /// \param domainMap [in] The graph's domain Map. MUST be one to
    ///   one!
    ///
    /// \param rangeMap [in] The graph's range Map.  MUST be one to
    ///   one!  May be, but need not be, the same as the domain Map.
    ///
    /// \param lclGraph [in] A locally indexed Kokkos::StaticCrsGraph
    ///   whose local row indices come from the specified row Map, and
    ///   whose local column indices come from the specified column
    ///   Map.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    TpetraCrsGraph(const local_graph_type& lclGraph,
                   const Teuchos::RCP<const map_type>& rowMap,
                   const Teuchos::RCP<const map_type>& colMap,
                   const Teuchos::RCP<const map_type>& domainMap = Teuchos::null,
                   const Teuchos::RCP<const map_type>& rangeMap = Teuchos::null,
                   const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(),
                                   typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(),
                                   "int",
                                   typeid(EpetraNode).name());
    }

    //! Destructor.
    virtual ~TpetraCrsGraph() {  }

    //@}

    //! @name Insertion/Removal Methods
    //@{

    //! Insert global indices into the graph.
    void insertGlobalIndices(GlobalOrdinal globalRow, const ArrayView< const GlobalOrdinal > &indices) { }

    //! Insert local indices into the graph.
    void insertLocalIndices(const LocalOrdinal localRow, const ArrayView< const LocalOrdinal > &indices) { }

    //! Remove all graph indices from the specified local row.
    void removeLocalIndices(LocalOrdinal localRow) { }

    //! Allocates the 1D pointer arrays of the graph
    void allocateAllIndices(size_t numNonZeros,ArrayRCP<size_t> & rowptr, ArrayRCP<LocalOrdinal> & colind){ }

    //! Sets the 1D pointer arrays of the graph.
    void setAllIndices(const ArrayRCP<size_t> & rowptr, const ArrayRCP<LocalOrdinal> & colind){ }

    //! Gets the 1D pointer arrays of the graph.
    void getAllIndices(ArrayRCP<const size_t>& rowptr, ArrayRCP<const LocalOrdinal>& colind) const{ }


    //@}

    //! @name Transformational Methods
    //@{

    //! Signal that data entry is complete, specifying domain and range maps.
    void fillComplete(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap, const RCP< ParameterList > &params=null) { }

    //! Signal that data entry is complete.
    void fillComplete(const RCP< ParameterList > &params=null) { }

    //! Expert version of fillComplete
    void expertStaticFillComplete (const Teuchos::RCP<const map_type>& domainMap,
                                   const Teuchos::RCP<const map_type>& rangeMap,
                                   const Teuchos::RCP<const Import< LocalOrdinal, GlobalOrdinal, Node > >& importer = null,
                                   const Teuchos::RCP<const Export< LocalOrdinal, GlobalOrdinal, Node > >& exporter = null,                          
                                   const Teuchos::RCP<Teuchos::ParameterList>& params=null){ } 

    //@}

    //! @name Methods implementing RowGraph.
    //@{

    //! Returns the communicator.
    RCP< const Comm< int > > getComm() const { return Teuchos::null; }

    //! Returns the Map that describes the row distribution in this graph.
    RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getRowMap() const { return Teuchos::null; }

    //! Returns the Map that describes the column distribution in this graph.
    RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getColMap() const { return Teuchos::null; }

    //! Returns the Map associated with the domain of this graph.
    RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getDomainMap() const { return Teuchos::null; }

    //! Returns the Map associated with the domain of this graph.
    RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getRangeMap() const { return Teuchos::null; }

    //! Returns the importer associated with this graph.
    RCP< const Import< LocalOrdinal, GlobalOrdinal, Node > > getImporter() const { return Teuchos::null; }

    //! Returns the exporter associated with this graph.
    RCP< const Export< LocalOrdinal, GlobalOrdinal, Node > > getExporter() const { return Teuchos::null; }

    //! Returns the number of global rows in the graph.
    global_size_t getGlobalNumRows() const { return 0; }

    //! Returns the number of global columns in the graph.
    global_size_t getGlobalNumCols() const { return 0; }

    //! Returns the number of graph rows owned on the calling node.
    size_t getLocalNumRows() const { return 0; }

    //! Returns the number of columns connected to the locally owned rows of this graph.
    size_t getLocalNumCols() const { return 0; }

    //! Returns the index base for global indices for this graph.
    GlobalOrdinal getIndexBase() const { return 0; }

    //! Returns the global number of entries in the graph.
    global_size_t getGlobalNumEntries() const { return 0; }

    //! Returns the local number of entries in the graph.
    size_t getLocalNumEntries() const { return 0; }

    //! Returns the current number of entries on this node in the specified global row.
    size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const { return 0; }

    //! Returns the current number of entries on this node in the specified local row.
    size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const { return 0; }

    //! Returns the current number of allocated entries for this node in the specified global row .
    size_t getNumAllocatedEntriesInGlobalRow(GlobalOrdinal globalRow) const { return 0; }

    //! Returns the current number of allocated entries on this node in the specified local row.
    size_t getNumAllocatedEntriesInLocalRow(LocalOrdinal localRow) const { return 0; }

    //! Maximum number of entries in all rows over all processes.
    size_t getGlobalMaxNumRowEntries() const { return 0; }

    //! Maximum number of entries in all rows owned by the calling process.
    size_t getLocalMaxNumRowEntries() const { return 0; }

    //! Whether the graph has a column Map.
    bool hasColMap() const { return false; }

    //! Whether column indices are stored using local indices on the calling process.
    bool isLocallyIndexed() const { return false; }

    //! Whether column indices are stored using global indices on the calling process.
    bool isGloballyIndexed() const { return false; }

    //! Whether fillComplete() has been called and the graph is in compute mode.
    bool isFillComplete() const { return false; }

    //! Returns true if storage has been optimized.
    bool isStorageOptimized() const { return false; }

    //! Return a const, nonpersisting view of global indices in the given row.
    void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView< const GlobalOrdinal > &Indices) const {  }

    //! Return a const, nonpersisting view of local indices in the given row.
    void getLocalRowView(LocalOrdinal LocalRow, ArrayView< const LocalOrdinal > &indices) const {  }

    /// \brief Access the local KokkosSparse::StaticCrsGraph data
    local_graph_type getLocalGraph () const {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                                 "Epetra does not support Kokkos::StaticCrsGraph!");
      TEUCHOS_UNREACHABLE_RETURN((local_graph_type()));
    }

    //! Dummy implementation for computeGlobalConstants
    void computeGlobalConstants() { }

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const { return std::string(""); }

    //! Print the object with some verbosity level to an FancyOStream object.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const {  }

    //@}

    //! @name Advanced methods, at increased risk of deprecation.
    //@{

    //! Get an ArrayRCP of the row-offsets.
    ArrayRCP< const size_t > getNodeRowPtrs() const { return Teuchos::ArrayRCP< const size_t>(); }

    //@}

    //! Implements DistObject interface
    //{@

    //! Access function for the Tpetra::Map this DistObject was constructed with.
    Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > getMap() const { return Teuchos::null; }

    //! Import.
    void doImport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &source,
                  const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM) { }

    //! Export.
    void doExport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &dest,
                  const Import< LocalOrdinal, GlobalOrdinal, Node >& importer, CombineMode CM) { }

    //! Import (using an Exporter).
    void doImport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &source,
                  const Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM) { }

    //! Export (using an Importer).
    void doExport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &dest,
                  const Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM) { }

    // @}

    //! @name Xpetra specific
    //@{

    //! TpetraCrsGraph constructor to wrap a Tpetra::CrsGraph object
    TpetraCrsGraph(const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > &graph)  {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "int", typeid(EpetraNode).name() );
    }

    //! Get the underlying Tpetra graph
    RCP< const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > getTpetra_CrsGraph() const { return Teuchos::null; }

    //@}
  }; // TpetraCrsGraph class (specialization for LO=GO=int and NO=EpetraNode)
#endif

#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))))

  // specialization of TpetraCrsGraph for GO=long long and NO=EpetraNode
  template <>
  class TpetraCrsGraph<int,long long,EpetraNode>
    : public CrsGraph<int,long long,EpetraNode>
  {
    typedef int LocalOrdinal;
    typedef long long GlobalOrdinal;
    typedef EpetraNode Node;

    // The following typedef is used by the XPETRA_DYNAMIC_CAST() macro.
    typedef TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node> TpetraCrsGraphClass;
    typedef Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

  public:

    //! @name Constructor/Destructor Methods
    //@{

    //! Constructor specifying fixed number of entries for each row.
    TpetraCrsGraph(const RCP< const map_type > &rowMap, size_t maxNumEntriesPerRow, const RCP< ParameterList > &params=null) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "long long", typeid(EpetraNode).name() );
    }

    //! Constructor specifying (possibly different) number of entries in each row.
    TpetraCrsGraph(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, const RCP< ParameterList > &params=null) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "long long", typeid(EpetraNode).name() );
    }

    //! Constructor specifying column Map and fixed number of entries for each row.
    TpetraCrsGraph(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, const RCP< ParameterList > &params=null) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "long long", typeid(EpetraNode).name() );
    }

    //! Constructor specifying column Map and number of entries in each row.
    TpetraCrsGraph(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, const RCP< ParameterList > &params=null) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "long long", typeid(EpetraNode).name() );
    }

    /// \brief Constructor specifying column Map and arrays containing the graph in sorted, local ids.
    ///
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///
    /// \param rowPointers [in] The beginning of each row in the graph,
    ///   as in a CSR "rowptr" array.  The length of this vector should be
    ///   equal to the number of rows in the graph, plus one.  This last
    ///   entry should store the nunber of nonzeros in the graph.
    ///
    /// \param columnIndices [in] The local indices of the columns,
    ///   as in a CSR "colind" array.  The length of this vector
    ///   should be equal to the number of unknowns in the graph.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    TpetraCrsGraph(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap,
                   const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap,
                   const typename local_graph_type::row_map_type& rowPointers,
                   const typename local_graph_type::entries_type::non_const_type& columnIndices,
                   const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(),
                                   typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(),
                                   "int",
                                   typeid(EpetraNode).name());
    }

    /// \brief Constructor specifying column Map and a local (sorted)
    ///   graph, which the resulting CrsGraph views.
    ///
    /// Unlike most other CrsGraph constructors, successful completion
    /// of this constructor will result in a fill-complete graph.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///
    /// \param lclGraph [in] A locally indexed Kokkos::StaticCrsGraph
    ///   whose local row indices come from the specified row Map, and
    ///   whose local column indices come from the specified column
    ///   Map.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    TpetraCrsGraph(const Teuchos::RCP<const map_type>& rowMap,
                   const Teuchos::RCP<const map_type>& colMap,
                   const local_graph_type& lclGraph,
                   const Teuchos::RCP<Teuchos::ParameterList>& params) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(),
                                   typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(),
                                   "int",
                                   typeid(EpetraNode).name());
    }

    /// \brief Constructor specifying column, domain and range maps, and a
    ///   local (sorted) graph, which the resulting CrsGraph views.
    ///
    /// Unlike most other CrsGraph constructors, successful completion
    /// of this constructor will result in a fill-complete graph.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///
    /// \param domainMap [in] The graph's domain Map. MUST be one to
    ///   one!
    ///
    /// \param rangeMap [in] The graph's range Map.  MUST be one to
    ///   one!  May be, but need not be, the same as the domain Map.
    ///
    /// \param lclGraph [in] A locally indexed Kokkos::StaticCrsGraph
    ///   whose local row indices come from the specified row Map, and
    ///   whose local column indices come from the specified column
    ///   Map.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    TpetraCrsGraph(const local_graph_type& lclGraph,
                   const Teuchos::RCP<const map_type>& rowMap,
                   const Teuchos::RCP<const map_type>& colMap,
                   const Teuchos::RCP<const map_type>& domainMap = Teuchos::null,
                   const Teuchos::RCP<const map_type>& rangeMap = Teuchos::null,
                   const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(),
                                   typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(),
                                   "int",
                                   typeid(EpetraNode).name());
    }

    /// \brief Constructor specifying column Map and arrays containing the graph in sorted, local ids.
    ///
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///
    /// \param rowPointers [in] The beginning of each row in the graph,
    ///   as in a CSR "rowptr" array.  The length of this vector should be
    ///   equal to the number of rows in the graph, plus one.  This last
    ///   entry should store the nunber of nonzeros in the graph.
    ///
    /// \param columnIndices [in] The local indices of the columns,
    ///   as in a CSR "colind" array.  The length of this vector
    ///   should be equal to the number of unknowns in the graph.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    TpetraCrsGraph(const Teuchos::RCP< const map_type > &rowMap,
                   const Teuchos::RCP< const map_type > &colMap,
                   const Teuchos::ArrayRCP<size_t>& rowPointers,
                   const Teuchos::ArrayRCP<LocalOrdinal>& columnIndices,
                   const Teuchos::RCP<Teuchos::ParameterList>& params) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(),
                                   typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(),
                                   "int",
                                   typeid(EpetraNode).name())
      }


    //! Destructor.
    virtual ~TpetraCrsGraph() {  }

    //@}

    //! @name Insertion/Removal Methods
    //@{

    //! Insert global indices into the graph.
    void insertGlobalIndices(GlobalOrdinal globalRow, const ArrayView< const GlobalOrdinal > &indices) { }

    //! Insert local indices into the graph.
    void insertLocalIndices(const LocalOrdinal localRow, const ArrayView< const LocalOrdinal > &indices) { }

    //! Remove all graph indices from the specified local row.
    void removeLocalIndices(LocalOrdinal localRow) { }

    
    //! Allocates the 1D pointer arrays of the graph
    void allocateAllIndices(size_t numNonZeros,ArrayRCP<size_t> & rowptr, ArrayRCP<LocalOrdinal> & colind){ }

    //! Sets the 1D pointer arrays of the graph.
    void setAllIndices(const ArrayRCP<size_t> & rowptr, const ArrayRCP<LocalOrdinal> & colind){ } 

    //! Gets the 1D pointer arrays of the graph.
    void getAllIndices(ArrayRCP<const size_t>& rowptr, ArrayRCP<const LocalOrdinal>& colind) const { } 


    //@}

    //! @name Transformational Methods
    //@{

    //! Signal that data entry is complete, specifying domain and range maps.
    void fillComplete(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap, const RCP< ParameterList > &params=null) { }

    //! Signal that data entry is complete.
    void fillComplete(const RCP< ParameterList > &params=null) { }

    //! Expert version of fillComplete
    void expertStaticFillComplete (const Teuchos::RCP<const map_type>& domainMap,
                                   const Teuchos::RCP<const map_type>& rangeMap,
                                   const Teuchos::RCP<const Import< LocalOrdinal, GlobalOrdinal, Node > >& importer = null,
                                   const Teuchos::RCP<const Export< LocalOrdinal, GlobalOrdinal, Node > >& exporter = null,                          
                                   const Teuchos::RCP<Teuchos::ParameterList>& params=null){ } 

    //@}

    //! @name Methods implementing RowGraph.
    //@{

    //! Returns the communicator.
    RCP< const Comm< int > > getComm() const { return Teuchos::null; }

    //! Returns the Map that describes the row distribution in this graph.
    RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getRowMap() const { return Teuchos::null; }

    //! Returns the Map that describes the column distribution in this graph.
    RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getColMap() const { return Teuchos::null; }

    //! Returns the Map associated with the domain of this graph.
    RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getDomainMap() const { return Teuchos::null; }

    //! Returns the Map associated with the domain of this graph.
    RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getRangeMap() const { return Teuchos::null; }

    //! Returns the importer associated with this graph.
    RCP< const Import< LocalOrdinal, GlobalOrdinal, Node > > getImporter() const { return Teuchos::null; }

    //! Returns the exporter associated with this graph.
    RCP< const Export< LocalOrdinal, GlobalOrdinal, Node > > getExporter() const { return Teuchos::null; }

    //! Returns the number of global rows in the graph.
    global_size_t getGlobalNumRows() const { return 0; }

    //! Returns the number of global columns in the graph.
    global_size_t getGlobalNumCols() const { return 0; }

    //! Returns the number of graph rows owned on the calling node.
    size_t getLocalNumRows() const { return 0; }

    //! Returns the number of columns connected to the locally owned rows of this graph.
    size_t getLocalNumCols() const { return 0; }

    //! Returns the index base for global indices for this graph.
    GlobalOrdinal getIndexBase() const { return 0; }

    //! Returns the global number of entries in the graph.
    global_size_t getGlobalNumEntries() const { return 0; }

    //! Returns the local number of entries in the graph.
    size_t getLocalNumEntries() const { return 0; }

    //! Returns the current number of entries on this node in the specified global row.
    size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const { return 0; }

    //! Returns the current number of entries on this node in the specified local row.
    size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const { return 0; }

    //! Returns the current number of allocated entries for this node in the specified global row .
    size_t getNumAllocatedEntriesInGlobalRow(GlobalOrdinal globalRow) const { return 0; }

    //! Returns the current number of allocated entries on this node in the specified local row.
    size_t getNumAllocatedEntriesInLocalRow(LocalOrdinal localRow) const { return 0; }

    //! Maximum number of entries in all rows over all processes.
    size_t getGlobalMaxNumRowEntries() const { return 0; }

    //! Maximum number of entries in all rows owned by the calling process.
    size_t getLocalMaxNumRowEntries() const { return 0; }

    //! Whether the graph has a column Map.
    bool hasColMap() const { return false; }

    //! Whether column indices are stored using local indices on the calling process.
    bool isLocallyIndexed() const { return false; }

    //! Whether column indices are stored using global indices on the calling process.
    bool isGloballyIndexed() const { return false; }

    //! Whether fillComplete() has been called and the graph is in compute mode.
    bool isFillComplete() const { return false; }

    //! Returns true if storage has been optimized.
    bool isStorageOptimized() const { return false; }

    //! Return a const, nonpersisting view of global indices in the given row.
    void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView< const GlobalOrdinal > &Indices) const {  }

    //! Return a const, nonpersisting view of local indices in the given row.
    void getLocalRowView(LocalOrdinal LocalRow, ArrayView< const LocalOrdinal > &indices) const {  }

    local_graph_type getLocalGraphDevice () const {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                                 "Epetra does not support Kokkos::StaticCrsGraph!");
      TEUCHOS_UNREACHABLE_RETURN((local_graph_type()));
    }

    typename local_graph_type::HostMirror getLocalGraphHost () const {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                                 "Epetra does not support Kokkos::StaticCrsGraph!");
      TEUCHOS_UNREACHABLE_RETURN((local_graph_type::HostMirror()));
    }

    //! Dummy implementation for computeGlobalConstants
    void computeGlobalConstants() { }

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const { return std::string(""); }

    //! Print the object with some verbosity level to an FancyOStream object.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const {  }

    //@}

    //! @name Advanced methods, at increased risk of deprecation.
    //@{

    //! Get an ArrayRCP of the row-offsets.
    ArrayRCP< const size_t > getNodeRowPtrs() const { return Teuchos::ArrayRCP< const size_t>(); }

    //@}

    //! Implements DistObject interface
    //{@

    //! Access function for the Tpetra::Map this DistObject was constructed with.
    Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > getMap() const { return Teuchos::null; }

    //! Import.
    void doImport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &source,
                  const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM) { }

    //! Export.
    void doExport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &dest,
                  const Import< LocalOrdinal, GlobalOrdinal, Node >& importer, CombineMode CM) { }

    //! Import (using an Exporter).
    void doImport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &source,
                  const Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM) { }

    //! Export (using an Importer).
    void doExport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &dest,
                  const Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM) { }

    // @}

    //! @name Xpetra specific
    //@{

    //! TpetraCrsGraph constructor to wrap a Tpetra::CrsGraph object
    TpetraCrsGraph(const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > &graph)  {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "long long", typeid(EpetraNode).name() );
    }

    //! Get the underlying Tpetra graph
    RCP< const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > getTpetra_CrsGraph() const { return Teuchos::null; }

    //@}
  }; // TpetraCrsGraph class (specialization for GO=long long and NO=EpetraNode)
#endif

#endif // HAVE_XPETRA_EPETRA


} // Xpetra namespace
#endif //XPETRA_TPETRACRSGRAPH_DEF_HPP

