#include "PB_SIMPLEPreconditionerFactory.hpp"

#include "PB_Utilities.hpp"
#include "PB_InverseFactory.hpp"
#include "PB_BlockLowerTriInverseOp.hpp"
#include "PB_BlockUpperTriInverseOp.hpp"

using Teuchos::RCP;

namespace PB {
namespace NS {

// Constructor definition
SIMPLEPreconditionerFactory
   ::SIMPLEPreconditionerFactory(const RCP<const InverseFactory> & inverse,
                                 double alpha)
   : inverse_(inverse), alpha_(alpha)
{ }

// Use the factory to build the preconditioner (this is where the work goes)
LinearOp SIMPLEPreconditionerFactory
   ::buildPreconditionerOperator(BlockedLinearOp & blockOp,
                                 BlockPreconditionerState & state) const
{
   int rows = blockRowCount(blockOp);
   int cols = blockColCount(blockOp);
 
   TEUCHOS_ASSERT(rows==2); // sanity checks
   TEUCHOS_ASSERT(cols==2);

   // extract subblocks
   const LinearOp F  = getBlock(0,0,blockOp);
   const LinearOp Bt = getBlock(0,1,blockOp);
   const LinearOp B  = getBlock(1,0,blockOp);
   const LinearOp C  = getBlock(1,1,blockOp);

   // get inverse of diag(F)
   const LinearOp H = getInvDiagonalOp(F);

   // build approximate Schur complement: hatS = -C + B*H*Bt
   const LinearOp hatS = explicitAdd(scale(-1.0,C), 
                                     explicitMultiply(B,H,Bt));

   // build inverses for F and the approximate Schur complement
   const LinearOp invF = buildInverse(*inverse_,F);
   const LinearOp invS = buildInverse(*inverse_,hatS);

   std::vector<LinearOp> invDiag(2); // vector storing inverses

   // build lower triangular inverse matrix
   BlockedLinearOp L = zeroBlockedOp(blockOp);
   setBlock(1,0,L,B);
   L->endBlockFill();

   invDiag[0] = invF;
   invDiag[1] = scale(-1.0,invS);
   LinearOp invL = createBlockLowerTriInverseOp(L,invDiag);

   // build upper triangular matrix
   BlockedLinearOp U = zeroBlockedOp(blockOp);
   setBlock(0,1,U,scale(1.0/alpha_,multiply(H,Bt)));
   U->endBlockFill();

   invDiag[0] = identity(rangeSpace(invF));
   invDiag[1] = scale(alpha_,identity(rangeSpace(invS)));
   LinearOp invU = createBlockUpperTriInverseOp(U,invDiag);

   // std::stringstream ss;
   // ss << "SIMPLE Preconditioner: ( inv(F) = " << invF->description()
   //    << ", inv(S) = " << invS->description() << " )";

   // return implicit product operator
   return multiply(invU,invL,"SIMPLE");
}

} // end namespace NS
} // end namespace PB
