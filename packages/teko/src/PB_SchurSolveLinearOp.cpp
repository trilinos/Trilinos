#include "PB_SchurSolveLinearOp.hpp"

#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

using namespace Thyra;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::dyn_cast;

namespace PB {

// Thyra::LinearOpBase requirements
////////////////////////////////////////////////////////////////////////
SchurSolveLinearOp::SchurSolveLinearOp(const Teuchos::RCP<const Thyra::BlockedLinearOpBase<double> > & A,
                                       const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invA00,
                                       const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invS)
   : A_(A), invA00_(invA00), invS_(invS), 
     A10_(A->getBlock(1,0)), A01_(A->getBlock(0,1))
{
   RCP<const VectorSpaceBase<double> > productArray[2];
 
   // create and store range space
   productArray[0] = invA00_->range();
   productArray[1] = invS_->range();
   productRange_   = rcp(new DefaultProductVectorSpace<double>(2,&productArray[0]));

   // create and store range space
   productArray[0] = invA00_->domain();
   productArray[1] = invS_->domain(); 
   productDomain_  = rcp(new DefaultProductVectorSpace<double>(2,&productArray[0]));
      // "productVectorSpace" did not yield to me quick enough
}

void SchurSolveLinearOp::apply( const EConj conj, const MultiVectorBase<double> &X,
      MultiVectorBase<double> *Y, const double alpha, const double beta) const 
{
   typedef RCP<MultiVectorBase<double> > MultiVector;
   typedef RCP<const MultiVectorBase<double> > ConstMultiVector;

   TEUCHOS_ASSERT(conj==Thyra::NONCONJ_ELE);

   // convert source and destination vectors to "ProductMultiVectorBase"s
   const ProductMultiVectorBase<double> & src = dyn_cast<const ProductMultiVectorBase<double> >(X);
   ProductMultiVectorBase<double> & dst       = dyn_cast<ProductMultiVectorBase<double> >(*Y);

   // get src blocks
   ConstMultiVector f = src.getMultiVectorBlock(0); // f
   ConstMultiVector g = src.getMultiVectorBlock(1); // g
   
   // get destination blocks
   MultiVector u = dst.getNonconstMultiVectorBlock(0); // u (u^)
   MultiVector p = dst.getNonconstMultiVectorBlock(1); // p (p^)
   MultiVector ps = g->clone_mv(); // this is need b/c of p = -inv(S)*p

   // set temporary operator for performing inv(A_00)*A_01
   const Teuchos::RCP<const Thyra::LinearOpBase<double> > invA00_A01 = Thyra::multiply(invA00_,A01_);

   // compute actual product
   Thyra::apply(*invA00_, conj,  *f, &*u);                 // u   = inv(A_00) * f
   Thyra::apply(*A10_, conj,  *u, &*ps, -1.0, 1.0);        // ps += -A_10*u
   Thyra::apply(*invS_, conj, *ps,  &*p, -1.0);            // p   = -inv(S)*ps
   Thyra::apply( *invA00_A01, conj,  *p,  &*u, -1.0, 1.0); // u  += -inv(A_00)*A_01*p
}

} // end namespace PB
