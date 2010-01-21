#include "Epetra/Teko_ReorderedMappingStrategy.hpp"

#include "PB_BlockedReordering.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcpFromRef;

namespace Teko {
namespace Epetra {

ReorderedMappingStrategy::ReorderedMappingStrategy(const BlockReorderManager & brm,const Teuchos::RCP<const MappingStrategy> & map)
   : reorderManager_(brm), mapStrategy_(map)
{
   rangeMap_ = mapStrategy_->rangeMap();
   domainMap_ = mapStrategy_->domainMap();
}

void ReorderedMappingStrategy::copyEpetraIntoThyra(const Epetra_MultiVector& X,
                                                 const Teuchos::Ptr<Thyra::MultiVectorBase<double> > & thyra_X,
                                                 const Teko::Epetra::EpetraOperatorWrapper & eow) const
{
   using Teuchos::ptr_const_cast;
   using Teuchos::rcp_const_cast;

   // first flatten the vector: notice this just works on the block structure
   RCP<Thyra::ProductMultiVectorBase<double> > prod_X = rcp_dynamic_cast<Thyra::ProductMultiVectorBase<double> >(rcpFromRef(*thyra_X));
   RCP<Thyra::MultiVectorBase<double> > flat_X = buildFlatMultiVector(reorderManager_,prod_X);

   // now use the underlying mapping strategy to copy the flat vector 
   mapStrategy_->copyEpetraIntoThyra(X,flat_X.ptr(),eow);
}

void ReorderedMappingStrategy::copyThyraIntoEpetra(const RCP<const Thyra::MultiVectorBase<double> > & thyra_Y,
                                                 Epetra_MultiVector& Y,
                                                 const Teko::Epetra::EpetraOperatorWrapper & eow) const
{
   // first flatten the vector: notice this just works on the block structure
   RCP<const Thyra::ProductMultiVectorBase<double> > prod_Y = rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<double> >(rcpFromRef(*thyra_Y));
   RCP<const Thyra::MultiVectorBase<double> > flat_Y = buildFlatMultiVector(reorderManager_,prod_Y);

   // now use the underlying mapping strategy to copy the flat vector 
   mapStrategy_->copyThyraIntoEpetra(flat_Y,Y,eow);
}

} // end namespace Epetra
} // end namespace Teko
