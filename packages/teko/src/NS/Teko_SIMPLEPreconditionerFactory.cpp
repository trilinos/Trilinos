#include "Teko_SIMPLEPreconditionerFactory.hpp"

#include "PB_Utilities.hpp"
#include "PB_InverseFactory.hpp"
#include "PB_BlockLowerTriInverseOp.hpp"
#include "PB_BlockUpperTriInverseOp.hpp"

#include "Teuchos_Time.hpp"

using Teuchos::RCP;

namespace Teko {
namespace NS {

// Constructor definition
SIMPLEPreconditionerFactory
   ::SIMPLEPreconditionerFactory(const RCP<InverseFactory> & inverse,
                                 double alpha)
   : invVelFactory_(inverse), invPrsFactory_(inverse), alpha_(alpha), fInverseType_(Diagonal), useMass_(false)
{ }

SIMPLEPreconditionerFactory
   ::SIMPLEPreconditionerFactory(const RCP<InverseFactory> & invVFact,
                                 const RCP<InverseFactory> & invPFact,
                                 double alpha)
   : invVelFactory_(invVFact), invPrsFactory_(invPFact), alpha_(alpha), fInverseType_(Diagonal), useMass_(false)
{ }

SIMPLEPreconditionerFactory::SIMPLEPreconditionerFactory()
   : alpha_(1.0), fInverseType_(Diagonal), useMass_(false)
{ }

// Use the factory to build the preconditioner (this is where the work goes)
LinearOp SIMPLEPreconditionerFactory
   ::buildPreconditionerOperator(BlockedLinearOp & blockOp,
                                 BlockPreconditionerState & state) const
{
   Teko_DEBUG_SCOPE("SIMPLEPreconditionerFactory::buildPreconditionerOperator",10);
   Teko_DEBUG_EXPR(Teuchos::Time timer(""));

   int rows = blockRowCount(blockOp);
   int cols = blockColCount(blockOp);
 
   TEUCHOS_ASSERT(rows==2); // sanity checks
   TEUCHOS_ASSERT(cols==2);

   bool buildExplicitSchurComplement = true;

   // extract subblocks
   const LinearOp F  = getBlock(0,0,blockOp);
   const LinearOp Bt = getBlock(0,1,blockOp);
   const LinearOp B  = getBlock(1,0,blockOp);
   const LinearOp C  = getBlock(1,1,blockOp);

   LinearOp matF = F;
   if(useMass_) {
      TEUCHOS_ASSERT(massMatrix_!=Teuchos::null);
      matF = massMatrix_;
   }

   // get approximation of inv(F) name H
   std::string fApproxStr = "<error>";
   LinearOp H;
   if(fInverseType_==NotDiag) {
      H = buildInverse(*customHFactory_,matF);
      fApproxStr = customHFactory_->toString();

      // since H is now implicit, we must build an implicit Schur complement
      buildExplicitSchurComplement = false;
   }
   else {
      // get generic diagonal
      H = getInvDiagonalOp(matF,fInverseType_);
      fApproxStr = getDiagonalName(fInverseType_);
   }

/*
   if(fInverseType_==Diagonal) {
      H = getInvDiagonalOp(matF);
      fApproxStr = "Diagonal";
   }
   else if(fInverseType_==Lumped) {
      H = getInvLumpedMatrix(matF);
      fApproxStr = "Lumped";
   }
   else if(fInverseType_==AbsRowSum) {
      H = getAbsRowSumInvMatrix(matF);
      fApproxStr = "AbsRowSum";
   }
   else if(fInverseType_==Custom) {
      H = buildInverse(*customHFactory_,matF);
      fApproxStr = customHFactory_->toString();

      // since H is now implicit, we must build an implicit Schur complement
      buildExplicitSchurComplement = false;
   }
   else {
      TEUCHOS_ASSERT(false);
   }
*/

   // adjust H for time scaling if it is a mass matrix
   if(useMass_) {
      RCP<const Teuchos::ParameterList> pl = state.getParameterList();

      if(pl->isParameter("stepsize")) {
         // get the step size
         double stepsize = pl->get<double>("stepsize");

         // scale by stepsize only if it is larger than 0
         if(stepsize>0.0)
            H = scale(stepsize,H);
      }
   }

   // build approximate Schur complement: hatS = -C + B*H*Bt
   LinearOp HBt, hatS;

   if(buildExplicitSchurComplement) {
      HBt = explicitMultiply(H,Bt);
      hatS = explicitAdd(C,scale(-1.0,explicitMultiply(B,HBt)));
   }
   else {
      // build an implicit Schur complement
      HBt = multiply(H,Bt);
      hatS = add(C,scale(-1.0,multiply(B,HBt)));
   }

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
   return multiply(invU,invL,"SIMPLE_"+fApproxStr);
}

//! Initialize from a parameter list
void SIMPLEPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & pl)
{
   RCP<const InverseLibrary> invLib = getInverseLibrary();

   // default conditions
   useMass_ = false;
   customHFactory_ = Teuchos::null;
   fInverseType_ = Diagonal;
  
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

      // build inverse types
      fInverseType_ = getDiagonalType(fInverseStr);
      if(fInverseType_==NotDiag)
         customHFactory_ = invLib->getInverseFactory(fInverseStr);

      /*
      if(fInverseStr=="Diagonal")
         fInverseType_ = Diagonal;
      else if(fInverseStr=="Lumped")
         fInverseType_ = Lumped;
      else if(fInverseStr=="AbsRowSum")
         fInverseType_ = AbsRowSum;
      else {
         fInverseType_ = Custom;
         customHFactory_ = invLib->getInverseFactory(fInverseStr);
      } 
      */
   }
   if(pl.isParameter("Use Mass Scaling"))
      useMass_ = pl.get<bool>("Use Mass Scaling");

   Teko_DEBUG_MSG_BEGIN(5)
      DEBUG_STREAM << "SIMPLE Parameters: " << std::endl;
      DEBUG_STREAM << "   inv type   = \"" << invStr  << "\"" << std::endl;
      DEBUG_STREAM << "   inv v type = \"" << invVStr << "\"" << std::endl;
      DEBUG_STREAM << "   inv p type = \"" << invPStr << "\"" << std::endl;
      DEBUG_STREAM << "   alpha      = " << alpha_ << std::endl;
      DEBUG_STREAM << "   use mass   = " << useMass_ << std::endl;
      DEBUG_STREAM << "SIMPLE Parameter list: " << std::endl;
      pl.print(DEBUG_STREAM);
   Teko_DEBUG_MSG_END()

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

//! For assiting in construction of the preconditioner
Teuchos::RCP<Teuchos::ParameterList> SIMPLEPreconditionerFactory::getRequestedParameters() const 
{
   Teuchos::RCP<Teuchos::ParameterList> result;
   Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());

   // grab parameters from F solver
   RCP<Teuchos::ParameterList> vList = invVelFactory_->getRequestedParameters();
   if(vList!=Teuchos::null) {
      Teuchos::ParameterList::ConstIterator itr;
      for(itr=vList->begin();itr!=vList->end();++itr)
         pl->setEntry(itr->first,itr->second);
      result = pl;
   }

   // grab parameters from S solver
   RCP<Teuchos::ParameterList> pList = invPrsFactory_->getRequestedParameters();
   if(pList!=Teuchos::null) {
      Teuchos::ParameterList::ConstIterator itr;
      for(itr=pList->begin();itr!=pList->end();++itr)
         pl->setEntry(itr->first,itr->second);
      result = pl;
   }

   // grab parameters from S solver
   if(customHFactory_!=Teuchos::null) {
      RCP<Teuchos::ParameterList> hList = customHFactory_->getRequestedParameters();
      if(hList!=Teuchos::null) {
         Teuchos::ParameterList::ConstIterator itr;
         for(itr=hList->begin();itr!=hList->end();++itr)
            pl->setEntry(itr->first,itr->second);
         result = pl;
      }
   }

   // use the mass matrix
   if(useMass_) {
      pl->set<Teko::LinearOp>("Velocity Mass Matrix", Teuchos::null,"Velocity mass matrix");
      result = pl;
   }

   return result;
}

//! For assiting in construction of the preconditioner
bool SIMPLEPreconditionerFactory::updateRequestedParameters(const Teuchos::ParameterList & pl) 
{
   Teko_DEBUG_SCOPE("InvLSCStrategy::updateRequestedParameters",10);
   bool result = true;
 
   // update requested parameters in solvers
   result &= invVelFactory_->updateRequestedParameters(pl);
   result &= invPrsFactory_->updateRequestedParameters(pl);
   if(customHFactory_!=Teuchos::null)
      result &= customHFactory_->updateRequestedParameters(pl);

   // set the mass matrix: throw if the strategy is not the right type
   if(useMass_) {
      Teko::LinearOp mass = pl.get<Teko::LinearOp>("Velocity Mass Matrix");

      // we must have a mass matrix
      if(mass==Teuchos::null)
         result &= false;
      else
         setMassMatrix(mass);
   }

   return result;
}

} // end namespace NS
} // end namespace Teko
