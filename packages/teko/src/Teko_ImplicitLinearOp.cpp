#include "Teko_ImplicitLinearOp.hpp"

namespace Teko {

using Teuchos::rcpFromRef;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcp_const_cast;
using Teuchos::RCP;

using Thyra::MultiVectorBase;

void ImplicitLinearOp::apply(const Thyra::EConj conj, 
           const Thyra::MultiVectorBase<double> & x, Thyra::MultiVectorBase<double> * y,
           const double alpha, const double beta) const
{
   TEST_FOR_EXCEPTION(conj!=Thyra::NONCONJ_ELE,std::runtime_error,
                           "Linear operators of inherited type Teko::ImplicitLinearOp "
                           "cannot handle conjugation (yet!)");

   MultiVector srcX = rcp_const_cast<MultiVectorBase<double> >(rcpFromRef(x));
   MultiVector destY = rcp_dynamic_cast<MultiVectorBase<double> >(rcpFromRef(*y));

   // call apply
   implicitApply(srcX,destY,alpha,beta);
}
 
} // end namespace Teko
