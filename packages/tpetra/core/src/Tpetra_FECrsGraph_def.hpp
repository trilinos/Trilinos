// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_FECRSGRAPH_DEF_HPP
#define TPETRA_FECRSGRAPH_DEF_HPP

#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_getEntryOnHost.hpp"

#include <type_traits>
#include <set>

namespace Tpetra {

template<class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
FECrsGraph(const Teuchos::RCP<const map_type> & ownedRowMap,
           const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,
           const size_t maxNumEntriesPerRow,
           const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter,
           const Teuchos::RCP<const map_type> & domainMap,
           const Teuchos::RCP<const map_type> & ownedRangeMap,
           const Teuchos::RCP<Teuchos::ParameterList>& params):
  FECrsGraph(ownedRowMap,ownedPlusSharedRowMap,maxNumEntriesPerRow,
             domainMap.is_null() ? ownedRowMap : domainMap,
             ownedPlusSharedToOwnedimporter,domainMap,ownedRangeMap,params)
{
  // Nothing else to do here
}


template<class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
FECrsGraph(const Teuchos::RCP<const map_type> & ownedRowMap,
           const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,
           const size_t maxNumEntriesPerRow,
           const Teuchos::RCP<const map_type> & ownedPlusSharedDomainMap,
           const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter,
           const Teuchos::RCP<const map_type> & ownedDomainMap,
           const Teuchos::RCP<const map_type> & ownedRangeMap,
           const Teuchos::RCP<Teuchos::ParameterList>& params):
  crs_graph_type(ownedPlusSharedRowMap, maxNumEntriesPerRow, params),
  ownedRowsImporter_(ownedPlusSharedToOwnedimporter),
  ownedDomainMap_(ownedDomainMap.is_null() ? ownedRowMap : ownedDomainMap),
  ownedRangeMap_(ownedRangeMap.is_null() ? ownedRowMap : ownedRangeMap)
{
  this->domainMap_  = ownedPlusSharedDomainMap.is_null() ? ownedPlusSharedRowMap : ownedPlusSharedDomainMap;
  Teuchos::RCP<const map_type> dummy;
  setup(ownedRowMap,ownedPlusSharedRowMap,dummy,params);
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
FECrsGraph (const Teuchos::RCP<const map_type> & ownedRowMap,
            const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,
            const Kokkos::DualView<const size_t*, device_type>& numEntPerRow,
            const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter,
            const Teuchos::RCP<const map_type> & domainMap,
            const Teuchos::RCP<const map_type> & ownedRangeMap,
            const Teuchos::RCP<Teuchos::ParameterList>& params):
  FECrsGraph(ownedRowMap,ownedPlusSharedRowMap,numEntPerRow,
             domainMap.is_null() ? ownedRowMap : domainMap,
             ownedPlusSharedToOwnedimporter,domainMap,ownedRangeMap,params)
{
  // Nothing else to do here
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
FECrsGraph (const Teuchos::RCP<const map_type> & ownedRowMap,
            const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,
            const Kokkos::DualView<const size_t*, device_type>& numEntPerRow,
            const Teuchos::RCP<const map_type> & ownedPlusSharedDomainMap,
            const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter,
            const Teuchos::RCP<const map_type> & ownedDomainMap,
            const Teuchos::RCP<const map_type> & ownedRangeMap,
            const Teuchos::RCP<Teuchos::ParameterList>& params):
  crs_graph_type( ownedPlusSharedRowMap, numEntPerRow, params),
  ownedRowsImporter_(ownedPlusSharedToOwnedimporter),
  ownedDomainMap_(ownedDomainMap.is_null() ? ownedRowMap : ownedDomainMap),
  ownedRangeMap_(ownedRangeMap.is_null() ? ownedRowMap : ownedRangeMap)
{
  this->domainMap_ = ownedPlusSharedDomainMap.is_null() ? ownedPlusSharedRowMap : ownedPlusSharedDomainMap;
  Teuchos::RCP<const map_type> dummy;
  setup(ownedRowMap,ownedPlusSharedRowMap,dummy,params);
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
FECrsGraph(const Teuchos::RCP<const map_type> & ownedRowMap,
           const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,
           const Teuchos::RCP<const map_type> & ownedPlusSharedColMap,
           const size_t maxNumEntriesPerRow,
           const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter,
           const Teuchos::RCP<const map_type> & domainMap,
           const Teuchos::RCP<const map_type> & ownedRangeMap,
           const Teuchos::RCP<Teuchos::ParameterList>& params):
  FECrsGraph(ownedRowMap,ownedPlusSharedRowMap,ownedPlusSharedColMap,maxNumEntriesPerRow,
             domainMap.is_null() ? ownedRowMap : domainMap,
             ownedPlusSharedToOwnedimporter, domainMap, ownedRangeMap, params)
{
  // Nothing else to do here
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
FECrsGraph(const Teuchos::RCP<const map_type> & ownedRowMap,
           const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,
           const Teuchos::RCP<const map_type> & ownedPlusSharedColMap,
           const size_t maxNumEntriesPerRow,
           const Teuchos::RCP<const map_type> & ownedPlusSharedDomainMap,
           const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter,
           const Teuchos::RCP<const map_type> & ownedDomainMap,
           const Teuchos::RCP<const map_type> & ownedRangeMap,
           const Teuchos::RCP<Teuchos::ParameterList>& params):
  crs_graph_type(ownedPlusSharedRowMap, ownedPlusSharedColMap,maxNumEntriesPerRow, params),
  ownedRowsImporter_(ownedPlusSharedToOwnedimporter),
  ownedDomainMap_(ownedDomainMap.is_null() ? ownedRowMap : ownedDomainMap),
  ownedRangeMap_(ownedRangeMap.is_null() ? ownedRowMap : ownedRangeMap)
{
  this->domainMap_ = ownedPlusSharedDomainMap.is_null() ? ownedPlusSharedRowMap : ownedPlusSharedDomainMap;
  setup(ownedRowMap,ownedPlusSharedRowMap, ownedPlusSharedColMap,params);
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
FECrsGraph (const Teuchos::RCP<const map_type> & ownedRowMap,
            const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,
            const Teuchos::RCP<const map_type> & ownedPlusSharedColMap,
            const Kokkos::DualView<const size_t*, device_type>& numEntPerRow,
            const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter,
            const Teuchos::RCP<const map_type> & domainMap,
            const Teuchos::RCP<const map_type> & ownedRangeMap,
            const Teuchos::RCP<Teuchos::ParameterList>& params):
  FECrsGraph(ownedRowMap,ownedPlusSharedRowMap,ownedPlusSharedColMap,numEntPerRow,
             domainMap.is_null() ? ownedRowMap : domainMap,
             ownedPlusSharedToOwnedimporter, domainMap, ownedRangeMap, params)
{
  // Nothing else to do here
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
FECrsGraph (const Teuchos::RCP<const map_type> & ownedRowMap,
            const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,
            const Teuchos::RCP<const map_type> & ownedPlusSharedColMap,
            const Kokkos::DualView<const size_t*, device_type>& numEntPerRow,
            const Teuchos::RCP<const map_type> & ownedPlusSharedDomainMap,
            const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter,
            const Teuchos::RCP<const map_type> & ownedDomainMap,
            const Teuchos::RCP<const map_type> & ownedRangeMap,
            const Teuchos::RCP<Teuchos::ParameterList>& params):
  crs_graph_type(ownedPlusSharedRowMap, ownedPlusSharedColMap, numEntPerRow, params),
  ownedRowsImporter_(ownedPlusSharedToOwnedimporter),
  ownedDomainMap_(ownedDomainMap.is_null() ? ownedRowMap : ownedDomainMap),
  ownedRangeMap_(ownedRangeMap.is_null() ? ownedRowMap : ownedRangeMap)
{
  this->domainMap_ = ownedPlusSharedDomainMap.is_null() ? ownedPlusSharedRowMap : ownedPlusSharedDomainMap;
  setup(ownedRowMap,ownedPlusSharedRowMap, ownedPlusSharedColMap,params);
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
setup(const Teuchos::RCP<const map_type>  & ownedRowMap,
      const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,
      const Teuchos::RCP<const map_type> & ownedPlusSharedColMap,
      const Teuchos::RCP<Teuchos::ParameterList>& params)
{
 const char tfecfFuncName[] = "FECrsGraph::setup(): ";

 TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(ownedRowMap.is_null (), std::runtime_error, "ownedRowMap is null.");
 TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(ownedPlusSharedRowMap.is_null (), std::runtime_error, "ownedPlusSharedRowMap is null.");

 // If we have a colMap, we're local, otherwise global
 if(ownedPlusSharedColMap.is_null()) this->allocateIndices(GlobalIndices);
 else this->allocateIndices(LocalIndices);

 activeCrsGraph_     = Teuchos::rcp(new FEWhichActive(FE_ACTIVE_OWNED_PLUS_SHARED));
 fillState_ = Teuchos::rcp(new FillState(FillState::closed));

 // Use a very strong map equivalence check
 bool maps_are_the_same = ownedRowMap->isSameAs(*ownedPlusSharedRowMap);
 if(!maps_are_the_same) {
   // Make an importer if we need to, check map compatability if we don't
   if(ownedRowsImporter_.is_null()) {
       ownedRowsImporter_ = Teuchos::rcp(new import_type(ownedRowMap,ownedPlusSharedRowMap));
   }
   else {
     TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(!ownedRowMap->isSameAs(*ownedRowsImporter_->getSourceMap()), std::runtime_error, "ownedRowMap does not match importer source map.");
     TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(!ownedPlusSharedRowMap->isSameAs(*ownedRowsImporter_->getTargetMap()), std::runtime_error, "ownedPlusSharedRowMap does not match importer target map.");
   }

   // Make sure the ownedPlusSharedRowMap has at least as many entries at the ownedRowMap (due to our superset requriement)
   TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( ownedRowsImporter_->getNumSameIDs() != ownedRowsImporter_->getSourceMap()->getLocalNumElements(),
                                          std::runtime_error,"ownedRowMap contains entries which are not in the ownedPlusSharedRowMap.");

   TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( ownedRowMap->getLocalNumElements() > ownedPlusSharedRowMap->getLocalNumElements(),
                                          std::runtime_error,"ownedRowMap more entries than the ownedPlusSharedRowMap.");

   // The locallyFitted check is debug mode only since it is more expensive
   const bool debug = ::Tpetra::Details::Behavior::debug ();
   if(debug) {
     TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( !ownedPlusSharedRowMap->isLocallyFitted(*ownedRowMap),
                                            std::runtime_error,"ownedPlusSharedRowMap must be locally fitted to the ownedRowMap");

   }
 } else {
   // We don't need the importer in this case, since the two row maps are identical (e.g., serial mode).
   // Setting this to null helps later to detect whether there is a need for the second graph (the owned one).
   ownedRowsImporter_ = Teuchos::null;
 }
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::doOwnedPlusSharedToOwned(const CombineMode CM) {
  const char tfecfFuncName[] = "FECrsGraph::doOwnedPlusSharedToOwned(CombineMode): ";
  if(!ownedRowsImporter_.is_null() && *activeCrsGraph_ == FE_ACTIVE_OWNED_PLUS_SHARED) {
    Teuchos::RCP<const map_type> ownedRowMap = ownedRowsImporter_->getSourceMap();

    // Do a self-export in "restricted mode"
    this->doExport(*this,*ownedRowsImporter_,CM,true);

    // Under the "if you own an element, you own at least one of its nodes" assumption,
    // we can start by making a columnmap for ownedPlusShared
    if(!this->hasColMap()) {
      Teuchos::Array<int> remotePIDs (0);
      this->makeColMap(remotePIDs);
    }

    // Now run CrsGraph's fillComplete to get the final importer
    crs_graph_type::fillComplete(this->domainMap_,this->getRowMap());

    // In debug mode, we check to make sure the "if you own an element, you own at least one of its nodes"
    // However, we give the user the ability to bypass this check. This can be useful when dealing with mesh
    // that are generated and partitioned by a TPL, and/or where it is so complicated to rearrange the
    // dofs ownership that the customer app might be willing to accept a small performance hit.
    const bool debug = ::Tpetra::Details::Behavior::debug ();
    const bool checkColGIDsInAtLeastOneOwnedRow =
        this->getMyNonconstParamList().is_null() ? true :
                  this->getMyNonconstParamList()->get("Check Col GIDs In At Least One Owned Row",true);
    if (debug && checkColGIDsInAtLeastOneOwnedRow) {
      Teuchos::RCP<const map_type> colmap = this->getColMap();
      Teuchos::Array<bool> flag(colmap->getLocalNumElements(),false);
      typename crs_graph_type::nonconst_local_inds_host_view_type indices("indices",this->getLocalMaxNumRowEntries());

      for(size_t i=0; i<ownedRowMap->getLocalNumElements(); i++)  {
        size_t NumEntries=0;
        this->getLocalRowCopy(i,indices,NumEntries);
        for(size_t j=0; j<NumEntries; j++)
          flag[indices[j]] = true;
      }

      int lclCount = 0;
      for(size_t i=0; i<(size_t)flag.size(); i++)
        if(!flag[i])
          ++lclCount;

      // Perform a reduction over the input comm, so that ranks with success=true won't be hanging.
      // Note: this only ensures things don't hang
      int gblCount = lclCount;
      auto comm = this->getComm();
      if (!comm.is_null()) {
        Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &lclCount, &gblCount);
      }
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (gblCount > 0,
         std::invalid_argument, "if you own an element (in the finite element sense) you "
         "must also own one of the attached nodes.  This assumption has been violated in "
         "your matrix fill on at least one MPI rank:\n"
         "  locally,  there are " + std::to_string(lclCount) + " col gids not connected to any owned gid.\n"
         "  globally, there are " + std::to_string(gblCount) + " col gids not connected to any owned gid.\n"
         "NOTE: you can disable this check by setting a parameter list with the option\n"
         "      'Check Col GIDs In At Least One Owned Row' set to false.\n"
         "NOTE: the parameter list must be set AFTER construction, since it would not be recognized as valid"
         "by the base class CrsGraph.\n");
    }

    // Time to build an owned localGraph via subviews
    inactiveCrsGraph_ = Teuchos::rcp(new crs_graph_type(*this, ownedRowMap));
    inactiveCrsGraph_->fillComplete(ownedDomainMap_,ownedRangeMap_);
  }
}//end doOverlapToLocal


template<class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::doOwnedToOwnedPlusShared(const CombineMode /* CM */) {
  // This should be a no-op for all of our purposes
}//end doLocalToOverlap

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::switchActiveCrsGraph() {
  if(*activeCrsGraph_ == FE_ACTIVE_OWNED_PLUS_SHARED)
    *activeCrsGraph_ = FE_ACTIVE_OWNED;
  else
    *activeCrsGraph_ = FE_ACTIVE_OWNED_PLUS_SHARED;

  if(inactiveCrsGraph_.is_null()) return;

  this->swap(*inactiveCrsGraph_);

}//end switchActiveCrsGraph


template<class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::endFill(
  const Teuchos::RCP<const map_type>& domainMap,
  const Teuchos::RCP<const map_type>& rangeMap)
{
  /* What has to go on here is complicated.
     First off, if we don't really have two graphs (e.g. the rowMaps are the same, because we're in serial or
     doing finite differences, things are easy --- just call fillComplete().

     If, we are in the parallel FE case, then:
     Precondition: FE_ACTIVE_OWNED_PLUS_SHARED mode

     Postconditions:
     1) FE_ACTIVE_OWNED mode
     2) The OWNED graph has been fillCompleted with an Aztec-compatible column map
     3) rowptr & (local) colinds are aliased between the two graphs
     4) The OWNED_PLUS_SHARED graph has been fillCompleted with a column map whose first chunk
        is the column map for the OWNED graph.
        If we assume that (a) if you own an element, you also own at least one of the connected nodes and (b) elements are cliques, then
        the columnMap is the same for both graphs!!!

     5) The OWNED_PLUS_SHARED graph has neither an importer nor exporter.  Making these is expensive and we don't need them.
   */
  // Precondition
  const char tfecfFuncName[] = "FECrsGraph::endFill(domainMap, rangeMap): ";
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(*activeCrsGraph_ != FE_ACTIVE_OWNED_PLUS_SHARED,std::runtime_error, "must be in owned+shared mode.");
  if(ownedRowsImporter_.is_null()) {
    // The easy case: One graph
    switchActiveCrsGraph();
    crs_graph_type::fillComplete(domainMap, rangeMap);
  }
  else {
    // The hard case: Two graphs

    // fillComplete the owned+shared graph in a way that generates the owned+shared grep w/o an importer or exporter
    // Migrate data to the owned graph
    doOwnedPlusSharedToOwned(Tpetra::ADD);

    // Load up the owned graph
    switchActiveCrsGraph();

  }
}


template<class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::endFill() {
  this->endFill(ownedDomainMap_, ownedRangeMap_);
}


template<class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::beginFill() {
  const char tfecfFuncName[] = "FECrsGraph::beginFill(): ";

  // Unlike FECrsMatrix and FEMultiVector, we do not allow you to call beginFill() after calling endFill()
  // So we throw an exception if you're in owned mode
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(*activeCrsGraph_ == FE_ACTIVE_OWNED,std::runtime_error, "can only be called once.");

}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::beginAssembly() {
  const char tfecfFuncName[] = "FECrsGraph::beginAssembly: ";
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    *fillState_ != FillState::closed,
    std::runtime_error,
    "Cannot beginAssembly, matrix is not in a closed state"
  );
  *fillState_ = FillState::open;
  this->beginFill();
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::endAssembly() {
  const char tfecfFuncName[] = "FECrsGraph::endAssembly: ";
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    *fillState_ != FillState::open,
    std::logic_error,
    "Cannot endAssembly, matrix is not open to fill but is closed."
  );
  *fillState_ = FillState::closed;
  this->endFill();
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::endAssembly(
  const Teuchos::RCP<const map_type>& domainMap,
  const Teuchos::RCP<const map_type>& rangeMap)
{
  const char tfecfFuncName[] = "FECrsGraph::endAssembly: ";
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    *fillState_ != FillState::open,
    std::logic_error,
    "Cannot endAssembly, matrix is not open to fill but is closed."
  );
  *fillState_ = FillState::closed;
  this->endFill(domainMap, rangeMap);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Teuchos::ParameterList>
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::getValidParameters () const
{
  auto valid_pl = Teuchos::rcp(new Teuchos::ParameterList("Tpetra::FECrsGraph"));
  valid_pl->validateParametersAndSetDefaults(*crs_graph_type::getValidParameters());
  valid_pl->set("Check Col GIDs In At Least One Owned Row",true);

  return valid_pl;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::insertGlobalIndicesImpl (
  const local_ordinal_type lclRow,
  const global_ordinal_type inputGblColInds[],
  const size_t numInputInds
){
  const char tfecfFuncName[] = "FECrsGraph::insertGlobalIndices: ";
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    *fillState_ != FillState::open,
    std::logic_error,
    "Cannot replace global values, matrix is not open to fill but is closed."
  );
  return crs_graph_type::insertGlobalIndicesImpl(lclRow, inputGblColInds, numInputInds);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::insertGlobalIndicesImpl (
  const RowInfo& rowInfo,
  const global_ordinal_type inputGblColInds[],
  const size_t numInputInds,
  std::function<void(const size_t, const size_t, const size_t)> fun
){
  const char tfecfFuncName[] = "FECrsGraph::insertGlobalIndices: ";
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    *fillState_ != FillState::open,
    std::logic_error,
    "Cannot replace global values, matrix is not open to fill but is closed."
  );
  return crs_graph_type::insertGlobalIndicesImpl(rowInfo, inputGblColInds, numInputInds, fun);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::insertLocalIndicesImpl (
  const local_ordinal_type lclRow,
  const Teuchos::ArrayView<const local_ordinal_type>& gblColInds,
  std::function<void(const size_t, const size_t, const size_t)> fun
){
  const char tfecfFuncName[] = "FECrsGraph::insertLocalIndices: ";
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    *fillState_ != FillState::open,
    std::logic_error,
    "Cannot replace global values, matrix is not open to fill but is closed."
  );
  return crs_graph_type::insertLocalIndicesImpl(lclRow, gblColInds, fun);
}

}  // end namespace Tpetra


//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_FECRSGRAPH_INSTANT(LO,GO,NODE) \
  template class FECrsGraph<LO, GO, NODE>;



#endif // TPETRA_FECRSGRAPH_DEF
