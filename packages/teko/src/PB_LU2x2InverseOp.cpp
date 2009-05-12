#include "PB_LU2x2InverseOp.hpp"

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
LU2x2InverseOp::LU2x2InverseOp(const BlockedLinearOp & A,
                                         const LinearOp & invA00,
                                         const LinearOp & invS)
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

void LU2x2InverseOp::implicitApply(const BlockedMultiVector & x, BlockedMultiVector & y,
                                const double alpha, const double beta) const
{
   // get src blocks
   MultiVector f = getBlock(0,x); // f
   MultiVector g = getBlock(1,x); // g
  
   // get extra storage
   MultiVector ps = deepcopy(g); // this is need b/c of p = -inv(S)*p
   
   // get destination blocks
   MultiVector u = getBlock(0,y); // u (u^)
   MultiVector p = getBlock(1,y); // p (p^)

   // for efficiency make copies of u and p
   MultiVector uc,pc; // copies of u and p
   if(beta!=0) {
      // perform a deep copy
      uc = deepcopy(u);
      pc = deepcopy(p);
   } else {
      // perform a shallow copy

      // uc and pc point 
      // to the data of u and p
      uc = u;
      pc = p;
   }

   // set temporary operator for performing inv(A_00)*A_01
   LinearOp invA00_A01 = Thyra::multiply(invA00_,A01_);

   // compute actual product
   applyOp(invA00_,     f,  uc);            // u   = inv(A_00) * f
   applyOp(A10_,       uc,  ps, -1.0, 1.0); // ps += -A_10*u
   applyOp(invS_,      ps,  pc, -1.0);      // p   = -inv(S)*ps
   applyOp(invA00_A01, pc,  uc, -1.0, 1.0); // u  += -inv(A_00)*A_01*p

   // scale result by alpha
   if(beta!=0) {
      update(alpha,uc,beta,u); // u = alpha * uc + beta * u
      update(alpha,pc,beta,p); // p = alpha * pc + beta * p
   } 
   else if(alpha!=1.0) {  
      scale(alpha,u); // u = alpha * u
      scale(alpha,p); // p = alpha * p
   }
}

} // end namespace PB
