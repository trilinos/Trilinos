#include "Teko_BlockImplicitLinearOp.hpp"

namespace Teko {

using Teuchos::rcpFromRef;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcp_const_cast;
using Teuchos::RCP;

using Thyra::ProductMultiVectorBase;

bool BlockImplicitLinearOp::opSupportedImpl(const Thyra::EOpTransp M_trans) const
{
  return (M_trans == Thyra::NOTRANS);
}

void BlockImplicitLinearOp::applyImpl(
  const Thyra::EOpTransp M_trans,
  const Thyra::MultiVectorBase<double> & x,
  const Teuchos::Ptr<Thyra::MultiVectorBase<double> > & y,
  const double alpha,
  const double beta
  ) const
{
   TEST_FOR_EXCEPTION(M_trans!=Thyra::NOTRANS, std::runtime_error,
     "Linear operators of inherited type BlockImplicitLinearOp "
     "cannot handle conjugation (yet!)");

   // cast source vector
   RCP<const ProductMultiVectorBase<double> > src =
     rcp_dynamic_cast<const ProductMultiVectorBase<double> >(rcpFromRef(x));
   BlockedMultiVector srcX = rcp_const_cast<ProductMultiVectorBase<double> >(src);

   // cast destination vector
   BlockedMultiVector destY =
     rcp_dynamic_cast<ProductMultiVectorBase<double> >(rcpFromPtr(y));

   // call apply
   implicitApply(srcX,destY,alpha,beta);
}
 
} // end namespace Teko
