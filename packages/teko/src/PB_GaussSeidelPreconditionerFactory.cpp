#include "PB_GaussSeidelPreconditionerFactory.hpp"

#include "PB_BlockUpperTriInverseOp.hpp"

using Teuchos::rcp;
using Teuchos::RCP;

namespace PB {

GaussSeidelPreconditionerFactory::GaussSeidelPreconditionerFactory(const LinearOp & invD0,const LinearOp & invD1)
      : invOpsStrategy_(rcp(new StaticInvDiagStrategy(invD0,invD1)))
{ }

GaussSeidelPreconditionerFactory::GaussSeidelPreconditionerFactory(const RCP<const BlockInvDiagonalStrategy> & strategy)
         : invOpsStrategy_(strategy)
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

   // create a blocked linear operator
   BlockedLinearOp U = getUpperTriBlocks(blo);

   return createNewBlockUpperTriInverseOp(U,invDiag);
}

} // end namspace PB
