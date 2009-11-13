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
   : invVelFactory_(inverse), invPrsFactory_(inverse), alpha_(alpha), fInverseType_(Diagonal)
{ }

SIMPLEPreconditionerFactory
   ::SIMPLEPreconditionerFactory(const RCP<InverseFactory> & invVFact,
                                 const RCP<InverseFactory> & invPFact,
                                 double alpha)
   : invVelFactory_(invVFact), invPrsFactory_(invPFact), alpha_(alpha), fInverseType_(Diagonal)
{ }

SIMPLEPreconditionerFactory::SIMPLEPreconditionerFactory()
   : alpha_(1.0), fInverseType_(Diagonal)
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

   // get approximation of inv(F) name H
   LinearOp H = getInvDiagonalOp(F);
   if(fInverseType_==Diagonal)
      H = getInvDiagonalOp(F);
   else if(fInverseType_==Lumped)
      H = getInvLumpedMatrix(F);
   else if(fInverseType_==AbsRowSum)
      H = getAbsRowSumInvMatrix(F);

   // build approximate Schur complement: hatS = -C + B*H*Bt
   const LinearOp HBt = explicitMultiply(H,Bt);
   const LinearOp hatS = explicitAdd(C,scale(-1.0,explicitMultiply(B,HBt)));

   // build the inverse for F 
   InverseLinearOp invF = state.getInverse("invF");
   if(invF==Teuchos::null) {
      invF = buildInverse(*invVelFactory_,F);
      state.addInverse("invF",invF); 
   } else {
      rebuildInverse(*invVelFactory_,F,invF);
   }

   // build the approximate Schur complement: This is inefficient! FIXME
   const LinearOp invS = buildInverse(*invPrsFactory_,hatS);

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
   alpha_ = 1.0;

   // "parse" the parameter list
   if(pl.isParameter("Inverse Type"))
      invStr = pl.get<std::string>("Inverse Type");
   if(pl.isParameter("Inverse Velocity Type"))
     invVStr = pl.get<std::string>("Inverse Velocity Type");
   if(pl.isParameter("Inverse Pressure Type")) 
     invPStr = pl.get<std::string>("Inverse Pressure Type");
   if(pl.isParameter("Alpha"))
     alpha_ = pl.get<double>("Alpha");
   if(pl.isParameter("Explicit Velocity Inverse Type")) {
     std::string fInverseStr = pl.get<std::string>("Explicit Velocity Inverse Type");
     if(fInverseStr=="Diagonal")
        fInverseType_ = Diagonal;
     else if(fInverseStr=="Lumped")
        fInverseType_ = Lumped;
     else if(fInverseStr=="AbsRowSum")
        fInverseType_ = AbsRowSum;
     else if(fInverseStr=="Custom")
     { TEST_FOR_EXCEPT(true); }
     else 
     { TEST_FOR_EXCEPT(true); }
   }

   PB_DEBUG_MSG_BEGIN(5)
      DEBUG_STREAM << "SIMPLE Parameters: " << std::endl;
      DEBUG_STREAM << "   inv type   = \"" << invStr  << "\"" << std::endl;
      DEBUG_STREAM << "   inv v type = \"" << invVStr << "\"" << std::endl;
      DEBUG_STREAM << "   inv p type = \"" << invPStr << "\"" << std::endl;
      DEBUG_STREAM << "   alpha    = " << alpha_ << std::endl;
      DEBUG_STREAM << "SIMPLE Parameter list: " << std::endl;
      pl.print(DEBUG_STREAM);
   PB_DEBUG_MSG_END()

   // set defaults as needed
   if(invStr=="") invStr = "Amesos";
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
}

} // end namespace NS
} // end namespace PB
