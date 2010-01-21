#include "Teko_BlockInvDiagonalStrategy.hpp"

namespace Teko {

InvFactoryDiagStrategy::InvFactoryDiagStrategy(const Teuchos::RCP<InverseFactory> & factory)
{
   // only one factory to use!
   invDiagFact_.resize(1,factory);
   defaultInvFact_ = factory;
}

InvFactoryDiagStrategy::InvFactoryDiagStrategy(const std::vector<Teuchos::RCP<InverseFactory> > & factories,
                                               const Teuchos::RCP<InverseFactory> & defaultFact)
{
   invDiagFact_ = factories;

   if(defaultFact==Teuchos::null)
      defaultInvFact_ = invDiagFact_[0];
   else
      defaultInvFact_ = defaultFact;
}

/** returns an (approximate) inverse of the diagonal blocks of A
  * where A is closely related to the original source for invD0 and invD1
  */
void InvFactoryDiagStrategy::getInvD(const BlockedLinearOp & A,BlockPreconditionerState & state,std::vector<LinearOp> & invDiag) const
{ 
   Teko_DEBUG_SCOPE("InvFactoryDiagStrategy::getInvD",10);

   // loop over diagonals, build an inverse operator for each
   int diagCnt = A->productRange()->numBlocks();
   int invCnt = invDiagFact_.size();

   Teko_DEBUG_MSG("# diags = " << diagCnt << ", # inverses = " << invCnt,6);

   if(diagCnt<=invCnt) {
      for(int i=0;i<diagCnt;i++) 
         invDiag.push_back(buildInverse(*invDiagFact_[i],getBlock(i,i,A)));
   }
   else {
      for(int i=0;i<invCnt;i++) 
         invDiag.push_back(buildInverse(*invDiagFact_[i],getBlock(i,i,A)));

      for(int i=invCnt;i<diagCnt;i++) 
         invDiag.push_back(buildInverse(*defaultInvFact_,getBlock(i,i,A)));
   }
}

} // end namespace Teko
