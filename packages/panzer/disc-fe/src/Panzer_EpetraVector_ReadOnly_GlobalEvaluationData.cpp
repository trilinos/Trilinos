#include "Panzer_EpetraVector_ReadOnly_GlobalEvaluationData.hpp"

#include "Epetra_Import.h"

#include "Thyra_VectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

namespace panzer {

using Teuchos::RCP;

void 
EpetraVector_ReadOnly_GlobalEvaluationData::
useConstantValues(const std::vector<int> & indices,double value)
{
  TEUCHOS_TEST_FOR_EXCEPTION(isInitialized_,std::logic_error,
                             "EpetraVector_ReadOnly_GED has been initialized, cannot call \"useConstantValues\"!");

  // add this specification to the filetered pairs vector
  FilteredPair pair;
  pair.first = indices;
  pair.second = value;
  filteredPairs_.push_back(pair);
}

void
EpetraVector_ReadOnly_GlobalEvaluationData::
initialize(const RCP<const Epetra_Import> & importer,
           const RCP<const Epetra_Map> & ghostedMap,
           const RCP<const Epetra_Map> & uniqueMap)
{
  importer_ = importer;
  ghostedMap_ = ghostedMap;
  uniqueMap_ = uniqueMap;

  // allocate the ghosted vector
  ghostedVector_ = Teuchos::rcp(new Epetra_Vector(*ghostedMap_));

  // build up the thyra conversion data structures
  ghostedSpace_ = Thyra::create_VectorSpace(ghostedMap_);
  uniqueSpace_ = Thyra::create_VectorSpace(uniqueMap_);

  // translate filtered pair GIDs to LIDs
   // initialize some ghosted values to the user specified values
   for(std::size_t i=0;i<filteredPairs_.size();i++) {
     std::vector<int> lids;
     const std::vector<int> & gids = filteredPairs_[i].first;
     for(std::size_t j=0;j<gids.size();j++) {
       int lid = ghostedMap->LID(gids[j]);

       // add legit LIDs to list
       if(lid>=0)
         lids.push_back(lid);
     }

     // overwrite original GID vector with new LID vector
     filteredPairs_[i].first = lids;
   }
  
  isInitialized_ = true;
}

void 
EpetraVector_ReadOnly_GlobalEvaluationData::
globalToGhost(int mem)
{
   TEUCHOS_TEST_FOR_EXCEPTION(uniqueVector_==Teuchos::null,std::logic_error,
                              "Unique vector has not been set, can't perform the halo exchange!");

   // initialize the ghosted data, zeroing out things, and filling in specified constants
   initializeData();

   Teuchos::RCP<const Epetra_Vector> uniqueVector_ep = Thyra::get_Epetra_Vector(*uniqueMap_,uniqueVector_);
  
   // do the global distribution
   ghostedVector_->Import(*uniqueVector_ep,*importer_,Insert);
}

void 
EpetraVector_ReadOnly_GlobalEvaluationData::
initializeData()
{
   TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_,std::logic_error,
                              "EpetraVector_ReadOnly_GED has not been initialized, cannot call \"initializeData\"!");

   ghostedVector_->PutScalar(0.0);

   // initialize some ghosted values to the user specified values
   for(std::size_t i=0;i<filteredPairs_.size();i++) {
     const std::vector<int> & lids = filteredPairs_[i].first;
     double value = filteredPairs_[i].second;
     for(std::size_t j=0;j<lids.size();j++)
       (*ghostedVector_)[lids[j]] = value;
   }
}

void 
EpetraVector_ReadOnly_GlobalEvaluationData::
setUniqueVector_Epetra(const Teuchos::RCP<const Epetra_Vector> & uniqueVector)
{
  TEUCHOS_ASSERT(isInitialized_);

  uniqueVector_ = Thyra::create_Vector(uniqueVector,uniqueSpace_);
}

Teuchos::RCP<Epetra_Vector> 
EpetraVector_ReadOnly_GlobalEvaluationData::
getGhostedVector_Epetra() const
{
  TEUCHOS_ASSERT(isInitialized_);
  TEUCHOS_ASSERT(ghostedVector_!=Teuchos::null);

  return ghostedVector_;
}

void 
EpetraVector_ReadOnly_GlobalEvaluationData::
setUniqueVector(const Teuchos::RCP<const Thyra::VectorBase<double> > & uniqueVector)
{
  TEUCHOS_ASSERT(isInitialized_);

  // uniqueVector_ep_ = Thyra::get_Epetra_Vector(*uniqueMap_,uniqueVector);
  uniqueVector_ = uniqueVector;
/*
  std::cout << "SETTING UNIQUE" << std::endl;
  std::cout << Teuchos::describe(*uniqueVector,Teuchos::VERB_EXTREME) << std::endl;
  uniqueVector_->Print(std::cout);
  std::cout << std::endl;
  */
}

Teuchos::RCP<const Thyra::VectorBase<double> > 
EpetraVector_ReadOnly_GlobalEvaluationData::
getUniqueVector() const
{
  TEUCHOS_ASSERT(isInitialized_);

  // return (uniqueVector_==Teuchos::null) ? Teuchos::null : Thyra::create_Vector(uniqueVector_,uniqueSpace_);
  return uniqueVector_;
}

Teuchos::RCP<Thyra::VectorBase<double> > 
EpetraVector_ReadOnly_GlobalEvaluationData::
getGhostedVector() const
{
  TEUCHOS_ASSERT(isInitialized_);
  TEUCHOS_ASSERT(ghostedVector_!=Teuchos::null);

  return Thyra::create_Vector(ghostedVector_,ghostedSpace_);
}

void
EpetraVector_ReadOnly_GlobalEvaluationData::
print(std::ostream & os) const
{
  const std::string tab = "    ";
  os << "\n";
  os << tab << "EpetraVector_ReadOnly_GlobalEvaluationData\n"
     << tab << "  init    = " << isInitialized_ << "\n"
     << tab << "  unique  = " << uniqueVector_ << "\n"
     << tab << "  ghosted = " << ghostedVector_ << "\n";

  /*
  os << "GHOSTED MAP\n";
  ghostedMap_->Print(os);
  os << "\n\n";

  os << "GHOSTED Vector\n";
  ghostedVector_->Print(os);
  os << "\n\n";

  os << "UNIQUE MAP\n";
  uniqueMap_->Print(os);
  os << "\n\n";

  os << "UNIQUE Vector\n";
  uniqueVector_->Print(os);
  os << "\n\n";
  */
}


}
