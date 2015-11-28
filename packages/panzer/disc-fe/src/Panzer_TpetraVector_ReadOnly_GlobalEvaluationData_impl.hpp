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
initialize(const RCP<const ImportType> & importer,
           const RCP<const MapType> & ghostedMap,
           const RCP<const MapType> & uniqueMap)
{
  importer_ = importer;
  ghostedMap_ = ghostedMap;
  uniqueMap_ = uniqueMap;

  // allocate the ghosted vector
  ghostedVector_ = Teuchos::rcp(new VectorType(ghostedMap_));

  // build up the thyra conversion data structures
  ghostedSpace_ = Thyra::tpetraVectorSpace<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(ghostedMap_);
  uniqueSpace_ = Thyra::tpetraVectorSpace<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(uniqueMap_);


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
globalToGhost(int mem)
{
   TEUCHOS_TEST_FOR_EXCEPTION(uniqueVector_==Teuchos::null,std::logic_error,
                              "Unique vector has not been set, can't perform the halo exchange!");

   // initialize the ghosted data, zeroing out things, and filling in specified constants
   initializeData();

   // do the global distribution
   ghostedVector_->doImport(*uniqueVector_,*importer_,Tpetra::INSERT);
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void 
TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
initializeData()
{
   TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_,std::logic_error,
                              "TpetraVector_ReadOnly_GED has not been initialized, cannot call \"initializeData\"!");

   ghostedVector_->putScalar(0.0);

   typedef typename VectorType::dual_view_type::t_dev::memory_space DMS;
   auto values_2d = ghostedVector_->template getLocalView<DMS>();
   auto values = Kokkos::subview(values_2d, Kokkos::ALL (), 0);

   // initialize some ghosted values to the user specified values
   for(std::size_t i=0;i<filteredPairs_.size();i++) {
     const std::vector<int> & lids = filteredPairs_[i].first;
     double value = filteredPairs_[i].second;
     for(std::size_t j=0;j<lids.size();j++)
       values(lids[j]) = value;
   }
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void 
TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
setUniqueVector_Tpetra(const Teuchos::RCP<const VectorType> & uniqueVector)
{
  TEUCHOS_ASSERT(isInitialized_);

  uniqueVector_ = uniqueVector;
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const typename TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::VectorType>
TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getUniqueVector_Tpetra() const
{
  TEUCHOS_ASSERT(isInitialized_);

  return uniqueVector_;
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
setUniqueVector(const Teuchos::RCP<const Thyra::VectorBase<double> > & uniqueVector)
{
  typedef Thyra::TpetraOperatorVectorExtraction<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> TOE;

  TEUCHOS_ASSERT(isInitialized_);
  uniqueVector_ = TOE::getConstTpetraVector(uniqueVector);
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Thyra::VectorBase<double> > 
TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getUniqueVector() const
{
  TEUCHOS_ASSERT(isInitialized_);

  return (uniqueVector_==Teuchos::null) ? Teuchos::null : Thyra::createConstVector(uniqueVector_,uniqueSpace_);
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
     << tab << "  unique  = " << uniqueVector_ << "\n"
     << tab << "  ghosted = " << ghostedVector_ << "\n";
}


}

#endif
