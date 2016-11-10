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
initialize(const RCP<const Epetra_Import>& importer,
           const RCP<const Epetra_Map>&    ghostedMap,
           const RCP<const Epetra_Map>&    ownedMap)
{
  importer_   = importer;
  ghostedMap_ = ghostedMap;
  ownedMap_   = ownedMap;

  // allocate the ghosted vector
  ghostedVector_ = Teuchos::rcp(new Epetra_Vector(*ghostedMap_));

  // build up the thyra conversion data structures
  ghostedSpace_ = Thyra::create_VectorSpace(ghostedMap_);
  ownedSpace_   = Thyra::create_VectorSpace(ownedMap_);

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
  TEUCHOS_TEST_FOR_EXCEPTION(ownedVector_ == Teuchos::null,std::logic_error,
    "Owned vector has not been set, can't perform the halo exchange!");

  // Initialize the ghosted data, zeroing out things, and filling in specified
  // constants.
  initializeData();
  Teuchos::RCP<const Epetra_Vector> ownedVector_ep = Thyra::get_Epetra_Vector(
    *ownedMap_, ownedVector_);
  
  // Do the global distribution.
  ghostedVector_->Import(*ownedVector_ep, *importer_, Insert);
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
setOwnedVector_Epetra(const Teuchos::RCP<const Epetra_Vector>& ownedVector)
{
  TEUCHOS_ASSERT(isInitialized_);
  ownedVector_ = Thyra::create_Vector(ownedVector, ownedSpace_);
}

///////////////////////////////////////////////////////////////////////////////
//
//  getOwnedVector_Epetra()
//
///////////////////////////////////////////////////////////////////////////////
Teuchos::RCP<const Epetra_Vector> EpetraVector_ReadOnly_GlobalEvaluationData::
getOwnedVector_Epetra() const
{
  TEUCHOS_ASSERT(isInitialized_);
  TEUCHOS_ASSERT(ownedVector_ != Teuchos::null);
  return Thyra::get_Epetra_Vector(*ownedMap_, ownedVector_);
} // end of getOwnedVector_Epetra()

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
setOwnedVector(const Teuchos::RCP<const Thyra::VectorBase<double> >&
	ownedVector)
{
  TEUCHOS_ASSERT(isInitialized_);
  // ownedVector_ep_ = Thyra::get_Epetra_Vector(*ownedMap_, ownedVector);
  ownedVector_ = ownedVector;
/*
  std::cout << "SETTING OWNED" << std::endl;
  std::cout << Teuchos::describe(*ownedVector, Teuchos::VERB_EXTREME) << std::endl;
  ownedVector_->Print(std::cout);
  std::cout << std::endl;
*/
}

Teuchos::RCP<const Thyra::VectorBase<double> > 
EpetraVector_ReadOnly_GlobalEvaluationData::
getOwnedVector() const
{
  TEUCHOS_ASSERT(isInitialized_);
  // return (ownedVector_==Teuchos::null) ? Teuchos::null :
  // 	 Thyra::create_Vector(ownedVector_, ownedSpace_);
  return ownedVector_;
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
     << tab << "  owned   = " << ownedVector_   << "\n"
     << tab << "  ghosted = " << ghostedVector_ << "\n";

  /*
  os << "GHOSTED MAP\n";
  ghostedMap_->Print(os);
  os << "\n\n";

  os << "GHOSTED Vector\n";
  ghostedVector_->Print(os);
  os << "\n\n";

  os << "OWNED MAP\n";
  ownedMap_->Print(os);
  os << "\n\n";

  os << "OWNED Vector\n";
  ownedVector_->Print(os);
  os << "\n\n";
  */
}


}
