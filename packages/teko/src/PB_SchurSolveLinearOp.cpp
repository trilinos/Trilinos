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
SchurSolveLinearOp::SchurSolveLinearOp(const Teuchos::RCP<const Thyra::BlockedLinearOpBase<double> > & fwdOp,
                                       const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invF,
                                       const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invS)
   : fwdOp_(fwdOp), invF_(invF), invS_(invS), 
     fwdL_(fwdOp->getBlock(1,0)), fwdU_(fwdOp->getBlock(0,1))
{
   RCP<const VectorSpaceBase<double> > productArray[2];
 
   // create and store range space
   productArray[0] = invF_->range();
   productArray[1] = invS_->range();
   productRange_   = rcp(new DefaultProductVectorSpace<double>(2,&productArray[0]));

   // create and store range space
   productArray[0] = invF_->domain();
   productArray[1] = invS_->domain(); 
   productDomain_  = rcp(new DefaultProductVectorSpace<double>(2,&productArray[0]));
      // "productVectorSpace" did not yield to me quick enough
}

void SchurSolveLinearOp::apply( const EConj conj, const MultiVectorBase<double> &X,
      MultiVectorBase<double> *Y, const double alpha, const double beta) const 
{
   typedef RCP<MultiVectorBase<double> > MultiVector;
   typedef RCP<const MultiVectorBase<double> > ConstMultiVector;

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

   // set temporary operator for performing inv(F)*G
   const Teuchos::RCP<const Thyra::LinearOpBase<double> > temp = Thyra::multiply(invF_,fwdU_);

   // compute actual product
   Thyra::apply(*invF_, conj,  *f, &*u);               // u   = inv(F) * f
   Thyra::apply(*fwdL_, conj,  *u, &*ps, -1.0, 1.0);   // ps += -L*u
   Thyra::apply(*invS_, conj, *ps,  &*p, -1.0);        // p   = -inv(S)*ps
   Thyra::apply( *temp, conj,  *p,  &*u, -1.0, 1.0);   // u  += -inv(F)*U*p
}

#if 0
// Thyra::BlockedLinearOpBase requirements
////////////////////////////////////////////////////////////////////////

/** \brief Return if the block <tt>(i,j)</tt> exists or not.
 *
 * \param  i  [in] Zero-based index for the block row.
 * \param  j  [in] Zero-based index for the block column.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>0 <= i && i < this->productRange()->numBlocks()</tt>
 * <li><tt>0 <= j && j < this->productDomain()->numBlocks()</tt>
 * </ul>
 */
bool SchurSolveLinearOp::blockExists(const int i, const int j) const 
{
   TEST_FOR_EXCEPT_MSG(true,"SchurSolveLinearOp::blockExists not implemented, this is a \"matrix-free\" operation");
   return false;
}
 
/** \brief Return if the block <tt>(i,j)</tt> is const only or not.
 *
 * \param  i  [in] Zero-based index for the block row.
 * \param  j  [in] Zero-based index for the block column.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>0 <= i && i < this->productRange()->numBlocks()</tt>
 * <li><tt>0 <= j && j < this->productDomain()->numBlocks()</tt>
 * </ul>
 */
bool SchurSolveLinearOp::blockIsConst(const int i, const int j) const
{
   TEST_FOR_EXCEPT_MSG(true,"SchurSolveLinearOp::blockIsConst exists not implemented, this is a \"matrix-free\" operation");
   return false;
}
 
/** \brief Return a non-const view of the block <tt>(i,j)</tt> if it exists.
 *
 * \param  i  [in] Zero-based index for the block row.
 * \param  j  [in] Zero-based index for the block column.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>this->blockIsConst(i,j)==false</tt>
 * <li><tt>0 <= i && i < this->productRange()->numBlocks()</tt>
 * <li><tt>0 <= j && j < this->productDomain()->numBlocks()</tt>
 * </ul>
 *
 * <b>Postconditions:</b><ul>
 * <li>[<tt>this->blockExists(i,j)==true</tt>] <tt>return.get()!=NULL</tt>
 * <li>[<tt>this->blockExists(i,j)==false</tt>] <tt>return.get()==NULL</tt>
 * </ul>
 */
Teuchos::RCP<LinearOpBase<double> > SchurSolveLinearOp::getNonconstBlock(const int i, const int j)
{
   TEST_FOR_EXCEPT_MSG(true,"SchurSolveLinearOp::getNonconstBlock not implemented, this is a \"matrix-free\" operation");
   return Teuchos::null;
}
 
/** \brief Return a const view of the block <tt>(i,j)</tt> if it exists.
 *
 * \param  i  [in] Zero-based index for the block row.
 * \param  j  [in] Zero-based index for the block column.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>0 <= i && i < this->productRange()->numBlocks()</tt>
 * <li><tt>0 <= j && j < this->productDomain()->numBlocks()</tt>
 * </ul>
 *
 * <b>Postconditions:</b><ul>
 * <li>[<tt>this->blockExists(i,j)==true</tt>] <tt>return.get()!=NULL</tt>
 * <li>[<tt>this->blockExists(i,j)==false</tt>] <tt>return.get()==NULL</tt>
 * </ul>
 */
Teuchos::RCP<const LinearOpBase<double> > SchurSolveLinearOp::getBlock(const int i, const int j) const
{
   TEST_FOR_EXCEPT_MSG(true,"SchurSolveLinearOp::getBlock not implemented, this is a \"matrix-free\" operation");
   return Teuchos::null;
}
#endif

} // end namespace PB
