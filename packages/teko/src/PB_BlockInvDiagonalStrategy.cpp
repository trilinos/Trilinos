#include "PB_BlockInvDiagonalStrategy.hpp"

namespace PB {

InvFactoryDiagStrategy::InvFactoryDiagStrategy(const Teuchos::RCP<InverseFactory> & factory)
{
   // only one factory to use!
   invDiagFact_.resize(1,factory);
}

InvFactoryDiagStrategy::InvFactoryDiagStrategy(const std::vector<Teuchos::RCP<InverseFactory> > & factories)
{
   invDiagFact_ = factories;
}

/** returns an (approximate) inverse of the diagonal blocks of A
  * where A is closely related to the original source for invD0 and invD1
  */
void InvFactoryDiagStrategy::getInvD(const BlockedLinearOp & A,BlockPreconditionerState & state,std::vector<LinearOp> & invDiag) const
{ 
   // loop over diagonals, build an inverse operator for each
   int diagCnt = A->productRange()->numBlocks();
   int invCnt = invDiagFact_.size();

   // make sure correct number of inverse factories exist
   TEUCHOS_ASSERT(invCnt==diagCnt || invCnt==1);
   
   if(invCnt==1) {
      for(int i=0;i<diagCnt;i++) 
         invDiag.push_back(buildInverse(*invDiagFact_[0],getBlock(i,i,A)));
   }
   else {
      for(int i=0;i<diagCnt;i++) 
         invDiag.push_back(buildInverse(*invDiagFact_[i],getBlock(i,i,A)));
   }
}

} // end namespace PB
