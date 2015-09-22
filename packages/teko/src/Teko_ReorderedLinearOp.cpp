#include "Teko_ReorderedLinearOp.hpp"

namespace Teko {

ReorderedLinearOp::ReorderedLinearOp(const Teuchos::RCP<const BlockReorderManager> & mgr,
                                     const Teuchos::RCP<Thyra::LinearOpBase<double> > & blockedOp)
   : mgr_(mgr), blockedOp_(blockedOp)
{
   
   range_ = buildFlatVectorSpace(*mgr_,blockedOp_->range());
   domain_ = buildFlatVectorSpace(*mgr_,blockedOp_->domain());
}

VectorSpace ReorderedLinearOp::range() const
{
   return range_;
}

VectorSpace ReorderedLinearOp::domain() const
{
   return domain_;
}

void ReorderedLinearOp::implicitApply(const MultiVector & x, MultiVector & y,
                                      const double alpha, const double beta) const
{
   using Teuchos::rcp_dynamic_cast;

   Teuchos::RCP<const Thyra::MultiVectorBase<double> > reorderX 
      = Teko::buildReorderedMultiVector(*mgr_,rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<double> >(x));
   MultiVector reorderY = Teko::buildReorderedMultiVector(*mgr_,rcp_dynamic_cast<Thyra::ProductMultiVectorBase<double> >(y));

   // this will automatically fill the right data
   Thyra::apply(*blockedOp_,Thyra::NOTRANS,*reorderX,reorderY.ptr(),alpha,beta);
}

} // end namespace Teko
