#include "PB_LSCPreconditionerFactory.hpp"

#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"

#include "PB_LU2x2InverseOp.hpp"
#include "PB_Utilities.hpp"
#include "PB_BlockUpperTriInverseOp.hpp"

namespace PB {
namespace NS {

using Teuchos::rcp;
using Teuchos::RCP;

using Thyra::multiply;
using Thyra::add;
using Thyra::identity;

// Stabilized constructor
LSCPreconditionerFactory::LSCPreconditionerFactory(const LinearOp & invF,const LinearOp & invBQBtmC,
                                                   const LinearOp & invD,const LinearOp & invMass)
      : invOpsStrategy_(rcp(new StaticLSCStrategy(invF,invBQBtmC,invD,invMass)))
{ }

// Stable constructor
LSCPreconditionerFactory::LSCPreconditionerFactory(const LinearOp & invF, const LinearOp & invBQBtmC,
                                                   const LinearOp & invMass)
      : invOpsStrategy_(rcp(new StaticLSCStrategy(invF,invBQBtmC,invMass)))
{ }

// fully generic constructor
LSCPreconditionerFactory::LSCPreconditionerFactory(const RCP<const LSCStrategy> & strategy)
   : invOpsStrategy_(strategy)
{ }

// for PreconditionerFactoryBase
///////////////////////////////////////////////////////////////////////

// initialize a newly created preconditioner object
LinearOp LSCPreconditionerFactory::buildPreconditionerOperator(BlockedLinearOp & blockOp,BlockPreconditionerState & state) const
{
   // extract sub-matrices from source operator 
   LinearOp F  = blockOp->getBlock(0,0);
   LinearOp B  = blockOp->getBlock(1,0);
   LinearOp Bt = blockOp->getBlock(0,1);

   // build what is neccessary for the state object
   invOpsStrategy_->buildState(blockOp,state);

   // extract operators from strategy
   LinearOp invF      = invOpsStrategy_->getInvF(blockOp,state);
   LinearOp invBQBtmC = invOpsStrategy_->getInvBQBt(blockOp,state);
   LinearOp invD      = invOpsStrategy_->getInvD(blockOp,state);

   // if necessary build an identity mass matrix
   LinearOp invMass   = invOpsStrategy_->getInvMass(blockOp,state);
   if(invMass==Teuchos::null)
      invMass = identity<double>(F->range());

   // need to build Schur complement,  inv(P) = inv(B*Bt)*(B*F*Bt)*inv(B*Bt)

   // first construct middle operator: M = B * inv(Mass) * F * inv(Mass) * Bt
   LinearOp M = 
      //          (B * inv(Mass) ) * F * (inv(Mass) * Bt)
      multiply( multiply(B,invMass), F , multiply(invMass,Bt));
      
   // now construct a linear operator schur complement
   LinearOp invPschur; 
   if(invD!=Teuchos::null)
      invPschur = add(multiply(invBQBtmC, M , invBQBtmC), invD);
   else
      invPschur = multiply(invBQBtmC, M , invBQBtmC);

   // build diagonal operations
   std::vector<LinearOp> invDiag(2);
   invDiag[0] = invF;
   invDiag[1] = Thyra::scale(-1.0,invPschur);

   // get upper triangular matrix
   BlockedLinearOp U = getUpperTriBlocks(blockOp); 


   // label the preconditioner
   std::stringstream ss;
   ss << "LSC Preconditioner " << (invD!=Teuchos::null ? "(stabilized)" : "(stable)" ) 
      << ": ( inv(F) = " << invF->description()
      << ", inv(BQBt) = " << invBQBtmC->description() << " )";

   // build the preconditioner operator
   return createBlockUpperTriInverseOp(U,invDiag,ss.str());
}

} // end namespace NS
} // end namespace PB
