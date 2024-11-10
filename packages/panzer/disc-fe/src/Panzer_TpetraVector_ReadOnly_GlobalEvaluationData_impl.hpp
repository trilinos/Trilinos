// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_TpetraVector_ReadOnly_GlobalEvaluationData_impl_hpp__
#define __Panzer_TpetraVector_ReadOnly_GlobalEvaluationData_impl_hpp__

#include "Thyra_TpetraThyraWrappers.hpp"

namespace panzer {

using Teuchos::RCP;

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void 
TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
useConstantValues(const std::vector<GlobalOrdinalT> & indices,double value)
{
  TEUCHOS_TEST_FOR_EXCEPTION(isInitialized_,std::logic_error,
                             "TpetraVector_ReadOnly_GED has been initialized, cannot call \"useConstantValues\"!");

  // add this specification to the filtered pairs vector
  FilteredGlobalPair pair;
  pair.first = indices;
  pair.second = value;
  globalFilteredPairs_.push_back(pair);
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void
TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
initialize(const RCP<const ImportType>& importer,
           const RCP<const MapType>&    ghostedMap,
           const RCP<const MapType>&    ownedMap)
{
  importer_   = importer;
  ghostedMap_ = ghostedMap;
  ownedMap_   = ownedMap;

  // allocate the ghosted vector
  ghostedVector_ = Teuchos::rcp(new VectorType(ghostedMap_));

  // build up the thyra conversion data structures
  ghostedSpace_ = Thyra::tpetraVectorSpace<ScalarT, LocalOrdinalT,
    GlobalOrdinalT, NodeT>(ghostedMap_);
  ownedSpace_   = Thyra::tpetraVectorSpace<ScalarT, LocalOrdinalT,
    GlobalOrdinalT, NodeT>(ownedMap_);


   // translate filtered pair GIDs to LIDs
   // initialize some ghosted values to the user specified values
   filteredPairs_.resize(globalFilteredPairs_.size());
   for(std::size_t i=0;i<globalFilteredPairs_.size();i++) {
     std::vector<LocalOrdinalT> lids;
     const std::vector<GlobalOrdinalT> & gids = globalFilteredPairs_[i].first;
     for(std::size_t j=0;j<gids.size();j++) {
       LocalOrdinalT lid = ghostedMap->getLocalElement(gids[j]);

       // add legit LIDs to list
       if(lid>=0)
         lids.push_back(lid);
     }

     // convert original GID vector to LID vector, store value as well
     filteredPairs_[i].first = lids;
     filteredPairs_[i].second = globalFilteredPairs_[i].second;
   }
  
  isInitialized_ = true;
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void 
TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
globalToGhost(int /* mem */)
{
  TEUCHOS_TEST_FOR_EXCEPTION(ownedVector_ == Teuchos::null, std::logic_error,
    "Owned vector has not been set, can't perform the halo exchange!");

  // Initialize the ghosted data, zeroing out things, and filling in specified
  // constants.
  initializeData();

  // Do the global distribution.
  ghostedVector_->doImport(*ownedVector_, *importer_, Tpetra::INSERT);
  PHX::ExecSpace().fence();
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void 
TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
initializeData()
{
   TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_,std::logic_error,
                              "TpetraVector_ReadOnly_GED has not been initialized, cannot call \"initializeData\"!");

   ghostedVector_->putScalar(0.0);
   PHX::ExecSpace().fence();

   auto values = ghostedVector_->getLocalViewHost(Tpetra::Access::OverwriteAll);

   // initialize some ghosted values to the user specified values
   for(std::size_t i=0;i<filteredPairs_.size();i++) {
     const std::vector<int> & lids = filteredPairs_[i].first;
     double value = filteredPairs_[i].second;
     for(std::size_t j=0;j<lids.size();j++)
       values(lids[j], 0) = value;
   }
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void 
TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
setOwnedVector_Tpetra(const Teuchos::RCP<const VectorType>& ownedVector)
{
  TEUCHOS_ASSERT(isInitialized_);
  ownedVector_ = ownedVector;
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const typename TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::VectorType>
TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getOwnedVector_Tpetra() const
{
  TEUCHOS_ASSERT(isInitialized_);
  return ownedVector_;
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<typename TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::VectorType>
TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedVector_Tpetra() const
{
  TEUCHOS_ASSERT(isInitialized_);
  TEUCHOS_ASSERT(ghostedVector_!=Teuchos::null);

  return ghostedVector_;
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void 
TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
setOwnedVector(const Teuchos::RCP<const Thyra::VectorBase<double> >&
	ownedVector)
{
  typedef Thyra::TpetraOperatorVectorExtraction<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> TOE;
  TEUCHOS_ASSERT(isInitialized_);
  ownedVector_ = TOE::getConstTpetraVector(ownedVector);
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Thyra::VectorBase<double> > 
TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getOwnedVector() const
{
  TEUCHOS_ASSERT(isInitialized_);
  return (ownedVector_ == Teuchos::null) ? Teuchos::null :
    Thyra::createConstVector(ownedVector_, ownedSpace_);
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Thyra::VectorBase<double> > 
TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedVector() const
{
  TEUCHOS_ASSERT(isInitialized_);
  TEUCHOS_ASSERT(ghostedVector_!=Teuchos::null);

  return Thyra::createVector(ghostedVector_,ghostedSpace_);
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void
TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
print(std::ostream & os) const
{
  const std::string tab = "    ";
  os << "\n";
  os << tab << "TpetraVector_ReadOnly_GlobalEvaluationData\n"
     << tab << "  init    = " << isInitialized_ << "\n"
     << tab << "  owned   = " << ownedVector_   << "\n"
     << tab << "  ghosted = " << ghostedVector_ << "\n";
}


}

#endif
