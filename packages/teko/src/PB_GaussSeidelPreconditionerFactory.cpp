#include "PB_GaussSeidelPreconditionerFactory.hpp"

#include "PB_BlockUpperTriInverseOp.hpp"
#include "PB_BlockLowerTriInverseOp.hpp"

using Teuchos::rcp;
using Teuchos::RCP;

namespace PB {

GaussSeidelPreconditionerFactory::GaussSeidelPreconditionerFactory(TriSolveType solveType,const LinearOp & invD0,const LinearOp & invD1)
      : invOpsStrategy_(rcp(new StaticInvDiagStrategy(invD0,invD1))), solveType_(solveType)
{ }

GaussSeidelPreconditionerFactory::GaussSeidelPreconditionerFactory(TriSolveType solveType,const RCP<const BlockInvDiagonalStrategy> & strategy)
         : invOpsStrategy_(strategy), solveType_(solveType)
{ }

LinearOp GaussSeidelPreconditionerFactory::buildPreconditionerOperator(BlockedLinearOp & blo) const
{
   int rows = blockRowCount(blo);
   int cols = blockColCount(blo);
 
   TEUCHOS_ASSERT(rows==cols);

   // get diagonal blocks
   std::vector<LinearOp> invDiag;
   invOpsStrategy_->getInvD(blo,invDiag);
   TEUCHOS_ASSERT(rows==invDiag.size());

   if(solveType_==GS_UseUpperTriangle) {
      // create a blocked linear operator
      BlockedLinearOp U = getUpperTriBlocks(blo);

      return createNewBlockUpperTriInverseOp(U,invDiag);
   } 
   else if(solveType_==GS_UseLowerTriangle) {
      // create a blocked linear operator
      BlockedLinearOp L = getLowerTriBlocks(blo);

      return createNewBlockLowerTriInverseOp(L,invDiag);
   }

   TEUCHOS_ASSERT(false); // we should never make it here!
}

} // end namspace PB
