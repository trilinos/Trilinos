#include "PB_JacobiPreconditionerFactory.hpp"

using Teuchos::rcp;

namespace PB {

JacobiPreconditionerFactory::JacobiPreconditionerFactory(const LinearOp & invD0,const LinearOp & invD1)
      : invOpsStrategy_(rcp(new StaticInvDiagStrategy(invD0,invD1)))
{ }

JacobiPreconditionerFactory::JacobiPreconditionerFactory(const RCP<const BlockInvDiagonalStrategy> & strategy)
         : invOpsStrategy_(strategy)
{ }

LinearOp JacobiPreconditionerFactory::buildPreconditionerOperator(BlockedLinearOp & blo,BlockPreconditionerState & state) const
{
   int rows = blo->productRange()->numBlocks();
   int cols = blo->productDomain()->numBlocks();
 
   TEUCHOS_ASSERT(rows==cols);

   // get diagonal blocks
   std::vector<LinearOp> invDiag;
   invOpsStrategy_->getInvD(blo,state,invDiag);
   TEUCHOS_ASSERT(rows==invDiag.size());

   // create a blocked linear operator
   BlockedLinearOp precond = createBlockedOp();
   std::stringstream ss;
   ss << "Jacobi Preconditioner ( ";

   // start filling the blocked operator
   beginBlockFill(precond,rows,rows); // this is assuming the matrix is square

   // build blocked diagonal matrix
   for(int i=0;i<rows;i++) {
      ss << " op" << i << " = " << invDiag[i]->description() << ", ";
      precond->setBlock(i,i,invDiag[i]);
   }
   ss << " )";
   
   endBlockFill(precond);
   // done filling the blocked operator

   precond->setObjectLabel(ss.str());
   
   return precond; 
}

} // end namspace PB
