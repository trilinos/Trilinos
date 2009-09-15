#include "PB_SIMPLEPreconditionerFactory.hpp"

#include "PB_Utilities.hpp"
#include "PB_InverseFactory.hpp"
#include "PB_BlockLowerTriInverseOp.hpp"
#include "PB_BlockUpperTriInverseOp.hpp"

#include "Teuchos_Time.hpp"

using Teuchos::RCP;

namespace PB {
namespace NS {

// Constructor definition
SIMPLEPreconditionerFactory
   ::SIMPLEPreconditionerFactory(const RCP<InverseFactory> & inverse,
                                 double alpha)
   : invVelFactory_(inverse), invPrsFactory_(inverse), alpha_(alpha)
{ }

SIMPLEPreconditionerFactory
   ::SIMPLEPreconditionerFactory(const RCP<InverseFactory> & invVFact,
                                 const RCP<InverseFactory> & invPFact,
                                 double alpha)
   : invVelFactory_(invVFact), invPrsFactory_(invPFact), alpha_(alpha)
{ }

SIMPLEPreconditionerFactory::SIMPLEPreconditionerFactory()
{ }

// Use the factory to build the preconditioner (this is where the work goes)
LinearOp SIMPLEPreconditionerFactory
   ::buildPreconditionerOperator(BlockedLinearOp & blockOp,
                                 BlockPreconditionerState & state) const
{
   PB_DEBUG_SCOPE("SIMPLEPreconditionerFactory::buildPreconditionerOperator",10);
   PB_DEBUG_EXPR(Teuchos::Time timer(""));

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
   PB_DEBUG_EXPR(timer.start(true));
   const LinearOp HBt = explicitMultiply(H,Bt);
   const LinearOp hatS = explicitAdd(C,scale(-1.0,explicitMultiply(B,HBt)));
   PB_DEBUG_EXPR(timer.stop());
   PB_DEBUG_MSG("SIMPLEPrecFact::buildPO Schur ConstTime = " << timer.totalElapsedTime(),2);

   // build the inverse for F 
   PB_DEBUG_EXPR(timer.start(true));
   InverseLinearOp invF = state.getInverse("invF");
   if(invF==Teuchos::null) {
      invF = buildInverse(*invVelFactory_,F);
      state.addInverse("invF",invF); 
   } else {
      rebuildInverse(*invVelFactory_,F,invF);
   }
   PB_DEBUG_EXPR(timer.stop());
   PB_DEBUG_MSG("SIMPLEPrecFact::buildPO GetInvFTime = " << timer.totalElapsedTime(),2);

   // build the approximate Schur complement: This is inefficient! FIXME
   PB_DEBUG_EXPR(timer.start(true));
   const LinearOp invS = buildInverse(*invPrsFactory_,hatS);
   PB_DEBUG_EXPR(timer.stop());
   PB_DEBUG_MSG("SIMPLEPrecFact::buildPO GetInvSTime = " << timer.totalElapsedTime(),2);

   std::vector<LinearOp> invDiag(2); // vector storing inverses

   // build lower triangular inverse matrix
   BlockedLinearOp L = zeroBlockedOp(blockOp);
   setBlock(1,0,L,B);
   endBlockFill(L);

   invDiag[0] = invF;
   invDiag[1] = invS;
   LinearOp invL = createBlockLowerTriInverseOp(L,invDiag);

   // build upper triangular matrix
   BlockedLinearOp U = zeroBlockedOp(blockOp);
   setBlock(0,1,U,scale(1.0/alpha_,HBt));
   endBlockFill(U);

   invDiag[0] = identity(rangeSpace(invF));
   invDiag[1] = scale(alpha_,identity(rangeSpace(invS)));
   LinearOp invU = createBlockUpperTriInverseOp(U,invDiag);

   // return implicit product operator
   return multiply(invU,invL,"SIMPLE");
}

//! Initialize from a parameter list
void SIMPLEPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & pl)
{
   RCP<const InverseLibrary> invLib = getInverseLibrary();

   // get string specifying inverse
   std::string invStr="", invVStr="", invPStr="";
   double alpha = 1.0;

   // "parse" the parameter list
   if(pl.isParameter("Inverse Type"))
      invStr = pl.get<std::string>("Inverse Type");
   if(pl.isParameter("Inverse Velocity Type"))
     invVStr = pl.get<std::string>("Inverse Velocity Type");
   if(pl.isParameter("Inverse Pressure Type")) 
     invPStr = pl.get<std::string>("Inverse Pressure Type");
   if(pl.isParameter("Alpha"))
     alpha = pl.get<double>("Alpha");

   // set defaults as needed
   if(invStr=="") invVStr = "Amesos";
   if(invVStr=="") invVStr = invStr;
   if(invPStr=="") invPStr = invStr;

   //  two inverse factory objects
   RCP<InverseFactory> invVFact, invPFact;

   // build velocity inverse factory
   invVFact = invLib->getInverseFactory(invVStr);
   invPFact = invVFact; // by default these are the same
   if(invVStr!=invPStr) // if different, build pressure inverse factory
      invPFact = invLib->getInverseFactory(invPStr);

   // based on parameter type build a strategy
   invVelFactory_ = invVFact; 
   invPrsFactory_ = invPFact;
   alpha_ = alpha;
}

} // end namespace NS
} // end namespace PB
