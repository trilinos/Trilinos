#include "PB_BlockInvDiagonalStrategy.hpp"

namespace PB {

InvFactoryDiagStrategy::InvFactoryDiagStrategy(const Teuchos::RCP<const InverseFactory> & factory)
{
   // only one factory to use!
   invDiagFact_ = factory;
}

/** returns an (approximate) inverse of the diagonal blocks of A
  * where A is closely related to the original source for invD0 and invD1
  */
void InvFactoryDiagStrategy::getInvD(const BlockedLinearOp & A,BlockPreconditionerState & state,std::vector<LinearOp> & invDiag) const
{ 
   // loop over diagonals, build an inverse operator for each
   int diagCnt = A->productRange()->numBlocks();
   for(int i=0;i<diagCnt;i++) 
      invDiag.push_back(buildInverse(*invDiagFact_,getBlock(i,i,A)));
}

} // end namespace PB
