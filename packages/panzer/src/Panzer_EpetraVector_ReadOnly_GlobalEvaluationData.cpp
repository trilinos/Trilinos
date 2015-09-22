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
  
  isInitialized_ = true;
}

void 
EpetraVector_ReadOnly_GlobalEvaluationData::
globalToGhost(int mem)
{
   TEUCHOS_TEST_FOR_EXCEPTION(uniqueVector_==Teuchos::null,std::logic_error,
                              "Unique vector has not been set, can't perform the halo exchange!");

   // do the global distribution
   ghostedVector_->PutScalar(0.0);
   ghostedVector_->Import(*uniqueVector_,*importer_,Insert);
}

void 
EpetraVector_ReadOnly_GlobalEvaluationData::
initializeData()
{
   TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_,std::logic_error,
                              "EpetraVector_ReadOnly_GED has not been initialized, cannot call \"initializeData\"!");

   ghostedVector_->PutScalar(0.0);
}

void 
EpetraVector_ReadOnly_GlobalEvaluationData::
setUniqueVector_Epetra(const Teuchos::RCP<const Epetra_Vector> & uniqueVector)
{
  TEUCHOS_ASSERT(isInitialized_);

  uniqueVector_ = uniqueVector;
}

Teuchos::RCP<const Epetra_Vector> 
EpetraVector_ReadOnly_GlobalEvaluationData::
getUniqueVector_Epetra() const
{
  TEUCHOS_ASSERT(isInitialized_);

  return uniqueVector_;
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

  uniqueVector_ = Thyra::get_Epetra_Vector(*uniqueMap_,uniqueVector);
}

Teuchos::RCP<const Thyra::VectorBase<double> > 
EpetraVector_ReadOnly_GlobalEvaluationData::
getUniqueVector() const
{
  TEUCHOS_ASSERT(isInitialized_);

  return (uniqueVector_==Teuchos::null) ? Teuchos::null : Thyra::create_Vector(uniqueVector_,uniqueSpace_);
}

Teuchos::RCP<Thyra::VectorBase<double> > 
EpetraVector_ReadOnly_GlobalEvaluationData::
getGhostedVector() const
{
  TEUCHOS_ASSERT(isInitialized_);
  TEUCHOS_ASSERT(ghostedVector_!=Teuchos::null);

  return Thyra::create_Vector(ghostedVector_,ghostedSpace_);
}

}
