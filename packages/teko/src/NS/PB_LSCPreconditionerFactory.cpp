#include "PB_LSCPreconditionerFactory.hpp"

#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

#include "PB_LU2x2InverseOp.hpp"
#include "PB_Utilities.hpp"
#include "PB_BlockUpperTriInverseOp.hpp"

#include "EpetraExt_RowMatrixOut.h"

namespace PB {
namespace NS {

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
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

LSCPreconditionerFactory::LSCPreconditionerFactory()
{ }

// for PreconditionerFactoryBase
///////////////////////////////////////////////////////////////////////

// initialize a newly created preconditioner object
LinearOp LSCPreconditionerFactory::buildPreconditionerOperator(BlockedLinearOp & blockOp,BlockPreconditionerState & state) const
{
   PB_DEBUG_MSG("BEGIN LSCPreconditionerFactory::buildPreconditionerOperator",10);

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

   // build the preconditioner operator: Use LDU or upper triangular approximation
   if(invOpsStrategy_->useFullLDU()) { 
      PB_DEBUG_MSG("END LSCPreconditionerFactory::buildPreconditionerOperator (Full LDU)",10);

      // solve using a full LDU decomposition
      return createLU2x2InverseOp(blockOp,invF,invPschur,"LSC-LDU");
   } else {
      // build diagonal operations
      std::vector<LinearOp> invDiag(2);
      invDiag[0] = invF;
      invDiag[1] = Thyra::scale(-1.0,invPschur);

      // get upper triangular matrix
      BlockedLinearOp U = getUpperTriBlocks(blockOp); 

      PB_DEBUG_MSG("END LSCPreconditionerFactory::buildPreconditionerOperator (Upper only)",10);

      // solve using only one inversion of F
      return createBlockUpperTriInverseOp(U,invDiag,"LSC-Upper");
   }
}

//! Initialize from a parameter list
void LSCPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & pl)
{
   RCP<const InverseLibrary> invLib = getInverseLibrary();

   // get string specifying inverse
   std::string invStr="", invVStr="", invPStr="";
   bool rowZeroing = true;

   // "parse" the parameter list
   if(pl.isParameter("Inverse Type"))
      invStr = pl.get<std::string>("Inverse Type");
   if(pl.isParameter("Inverse Velocity Type"))
     invVStr = pl.get<std::string>("Inverse Velocity Type");
   if(pl.isParameter("Inverse Pressure Type")) 
     invPStr = pl.get<std::string>("Inverse Pressure Type");
   if(pl.isParameter("Ignore Boundary Rows"))
     rowZeroing = pl.get<bool>("Ignore Boundary Rows");

   // set defaults as needed
   if(invStr=="") invVStr = "Amesos";
   if(invVStr=="") invVStr = invStr;
   if(invPStr=="") invPStr = invStr;

   //  two inverse factory objects
   RCP<const InverseFactory> invVFact, invPFact;

   // build velocity inverse factory
   invVFact = invLib->getInverseFactory(invVStr);
   invPFact = invVFact; // by default these are the same
   if(invVStr!=invPStr) // if different, build pressure inverse factory
      invPFact = invLib->getInverseFactory(invPStr);

   // based on parameter type build a strategy
   invOpsStrategy_ = rcp(new InvLSCStrategy(invVFact,invPFact,rowZeroing));
}

} // end namespace NS
} // end namespace PB
