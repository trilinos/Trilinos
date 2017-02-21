#include "Panzer_EpetraVector_Write_GlobalEvaluationData.hpp"

#include "Epetra_Export.h"

#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"

namespace panzer {

using Teuchos::RCP;

void
EpetraVector_Write_GlobalEvaluationData::
initialize(const RCP<const Epetra_Export>& exporter,
           const RCP<const Epetra_Map>&    ghostedMap,
           const RCP<const Epetra_Map>&    ownedMap)
{
  exporter_   = exporter;
  ghostedMap_ = ghostedMap;
  ownedMap_   = ownedMap;

  // allocate the ghosted vector
  ghostedVector_ = Teuchos::rcp(new Epetra_Vector(*ghostedMap_));

  // build up the thyra conversion data structures
  ghostedSpace_ = Thyra::create_VectorSpace(ghostedMap_);
  ownedSpace_   = Thyra::create_VectorSpace(ownedMap_);

  isInitialized_ = true;
}

void 
EpetraVector_Write_GlobalEvaluationData::
ghostToGlobal(int mem)
{
  TEUCHOS_TEST_FOR_EXCEPTION(ownedVector_ == Teuchos::null,std::logic_error,
    "Owned vector has not been set, can't perform the halo exchange!");

  // Initialize the ghosted data, zeroing out things, and filling in specified
  // constants.
  initializeData();
  Teuchos::RCP<Epetra_Vector> ownedVector_ep = Thyra::get_Epetra_Vector(
    *ownedMap_, ownedVector_);

  // set different combine modes
  Epetra_CombineMode cm = Add;
  switch(getCombineMode()) {
    case CM_Sum:
      cm = Add;
      break;
    case CM_Min:
      cm = Epetra_Min;
      break;
    case CM_Max:
      cm = Epetra_Max;
      break;
    case CM_Insert:
      cm = Insert;
      break;
  };
  
  // Do the global distribution.
  ownedVector_ep->Export(*ghostedVector_,*exporter_, cm);
}

void 
EpetraVector_Write_GlobalEvaluationData::
initializeData()
{
   TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_,std::logic_error,
                              "EpetraVector_Write_GED has not been initialized, cannot call \"initializeData\"!");

   Thyra::put_scalar(0.0,ownedVector_.ptr());
}

void 
EpetraVector_Write_GlobalEvaluationData::
setOwnedVector_Epetra(const Teuchos::RCP<Epetra_Vector>& ownedVector)
{
  TEUCHOS_ASSERT(isInitialized_);
  ownedVector_ = Thyra::create_Vector(ownedVector, ownedSpace_);
}

Teuchos::RCP<Epetra_Vector> 
EpetraVector_Write_GlobalEvaluationData::
getGhostedVector_Epetra() const
{
  TEUCHOS_ASSERT(isInitialized_);
  TEUCHOS_ASSERT(ghostedVector_!=Teuchos::null);

  return ghostedVector_;
}

void 
EpetraVector_Write_GlobalEvaluationData::
setOwnedVector(const Teuchos::RCP<Thyra::VectorBase<double> >&
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

Teuchos::RCP<Thyra::VectorBase<double> > 
EpetraVector_Write_GlobalEvaluationData::
getOwnedVector() const
{
  TEUCHOS_ASSERT(isInitialized_);
  // return (ownedVector_==Teuchos::null) ? Teuchos::null :
  // 	 Thyra::create_Vector(ownedVector_, ownedSpace_);
  return ownedVector_;
}

Teuchos::RCP<Thyra::VectorBase<double> > 
EpetraVector_Write_GlobalEvaluationData::
getGhostedVector() const
{
  TEUCHOS_ASSERT(isInitialized_);
  TEUCHOS_ASSERT(ghostedVector_!=Teuchos::null);

  return Thyra::create_Vector(ghostedVector_,ghostedSpace_);
}

void
EpetraVector_Write_GlobalEvaluationData::
print(std::ostream & os) const
{
  const std::string tab = "    ";
  os << "\n";
  os << tab << "EpetraVector_Write_GlobalEvaluationData\n"
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
