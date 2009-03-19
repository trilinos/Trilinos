#include "PB_LSCPreconditionerFactory.hpp"

// Thyra includes
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"

using namespace Thyra;
using namespace Teuchos;

namespace PB {
namespace NS {

// Stabilized constructor
LSCPreconditionerFactory::LSCPreconditionerFactory(const RCP<const LinearOpBase<double> > & invF,
                               const RCP<const LinearOpBase<double> > & invBQBtmC,
                               const RCP<const LinearOpBase<double> > & invD,
                               const RCP<const LinearOpBase<double> > & invMass)
      : Block2x2PreconditionerFactory(), invOpsStrategy_(Teuchos::null) 
{
   invOpsStrategy_ = rcp(new StaticLSCStrategy(invF,invBQBtmC,invD,invMass));
}

// Stable constructor
LSCPreconditionerFactory::LSCPreconditionerFactory(const RCP<const LinearOpBase<double> > & invF,
                               const RCP<const LinearOpBase<double> > & invBQBtmC,
                               const RCP<const LinearOpBase<double> > & invMass)
      : Block2x2PreconditionerFactory(), invOpsStrategy_(Teuchos::null) 
{
   invOpsStrategy_ = rcp(new StaticLSCStrategy(invF,invBQBtmC,invMass));
}

// fully generic constructor
LSCPreconditionerFactory::LSCPreconditionerFactory(const Teuchos::RCP<const LSCStrategy> & strategy)
   : invOpsStrategy_(strategy)
{ }

// for PreconditionerFactoryBase
///////////////////////////////////////////////////////////////////////

// initialize a newly created preconditioner object
void PB::NS::LSCPreconditionerFactory::initializePrec(
                    const RCP<const LinearOpSourceBase<double> > & opSrc,
                    PreconditionerBase<double> * prec,
                    const ESupportSolveUse supportSolveUse) const
{
   RCP<const BlockedLinearOpBase<double> > blockOp = rcp_dynamic_cast<const BlockedLinearOpBase<double> >(opSrc->getOp());
   // here we need to check that the block ranges match! (or maybe just let the runtime system do it! :) )

   // extract sub-matrices from source operator 
   const RCP<const LinearOpBase<double> > F  = blockOp->getBlock(0,0);
   const RCP<const LinearOpBase<double> > B  = blockOp->getBlock(1,0);
   const RCP<const LinearOpBase<double> > Bt = blockOp->getBlock(0,1);

   // extract operators from strategy
   const RCP<const LinearOpBase<double> > invF      = invOpsStrategy_->getInvF(blockOp);
   const RCP<const LinearOpBase<double> > invBQBtmC = invOpsStrategy_->getInvBQBt(blockOp);
   const RCP<const LinearOpBase<double> > invD      = invOpsStrategy_->getInvD(blockOp);

   // if necessary build an identity mass matrix
   RCP<const LinearOpBase<double> > invMass   = invOpsStrategy_->getInvMass(blockOp);
   if(invMass==Teuchos::null)
      invMass = Thyra::identity<double>(F->range());

   // need to build Schur complement,  inv(P) = inv(B*Bt)*(B*F*Bt)*inv(B*Bt)

   // first construct middle operator: M = B * inv(Mass) * F * inv(Mass) * Bt
   RCP<const LinearOpBase<double> > M = 
      //          (B * inv(Mass) ) * F * (inv(Mass) * Bt)
      multiply( multiply(B,invMass), F , multiply(invMass,Bt),"inv(Mass)*F*inv(Mass)");
      
   // now construct a linear operator schur complement
   RCP<const LinearOpBase<double> > invPschur; 
   if(invD!=Teuchos::null)
      invPschur = add(multiply(invBQBtmC, M , invBQBtmC,"inv(B*Bt)*(B*F*Bt)*inv(B*Bt)"),invD);
   else
      invPschur = multiply(invBQBtmC, M , invBQBtmC,"inv(B*Bt)*(B*F*Bt)*inv(B*Bt)");

   // build a preconditioner operator using the parent classes utility function
   const RCP<const LinearOpBase<double> > precOp = build2x2InverseOperator(blockOp,invF,invPschur);

   // must first cast that to be initialized
   DefaultPreconditioner<double> & dPrec = dyn_cast<DefaultPreconditioner<double> >(*prec);
   dPrec.initializeUnspecified(precOp); // should this be initializeUnspecified?
}

} // end namespace NS
} // end namespace PB
