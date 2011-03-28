/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
//  
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//  
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//  
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//  
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//  
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission. 
//  
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
// 
// ***********************************************************************
// 
// @HEADER

*/

#include "Teko_SIMPLEPreconditionerFactory.hpp"

#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_BlockLowerTriInverseOp.hpp"
#include "Teko_BlockUpperTriInverseOp.hpp"
#include "Teko_DiagonalPreconditionerFactory.hpp"

#include "Teuchos_Time.hpp"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_PointToBlockDiagPermute.h"
#include "EpetraExt_MatrixMatrix.h"
#include "Thyra_EpetraOperatorWrapper.hpp"
#include "Thyra_EpetraLinearOp.hpp"


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
   else if(fInverseType_==BlkDiag) {
     // Block diagonal approximation for H
     DiagonalPreconditionerFactory Hfact;
     DiagonalPrecondState Hstate;
     Hfact.initializeFromParameterList(BlkDiagList_);           
     H = Hfact.buildPreconditionerOperator(matF,Hstate); 

/*
     // Get a FECrsMarix out of the BDP
     RCP<Epetra_FECrsMatrix> Hcrs=rcp(Hstate.BDP_->CreateFECrsMatrix());
     H=Thyra::epetraLinearOp(Hcrs);
*/

     buildExplicitSchurComplement = true; // NTS: Do I need this? 
                                          // Answer - no, but it is documenting whats going on here.
   }
   else {
      // get generic diagonal
      H = getInvDiagonalOp(matF,fInverseType_);
      fApproxStr = getDiagonalName(fInverseType_);
   }

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
      ModifiableLinearOp & mHBt = state.getModifiableOp("HBt");
      ModifiableLinearOp & mhatS = state.getModifiableOp("hatS");
      ModifiableLinearOp & BHBt = state.getModifiableOp("BHBt");

      // build H*Bt
      mHBt = explicitMultiply(H,Bt,mHBt);
      HBt = mHBt;

      // build B*H*Bt
      BHBt = explicitMultiply(B,HBt,BHBt);

      // build C-B*H*Bt
      mhatS = explicitAdd(C,scale(-1.0,BHBt),mhatS);
      hatS = mhatS;
   }
   else {
      // build an implicit Schur complement
      HBt = multiply(H,Bt);
      
      hatS = add(C,scale(-1.0,multiply(B,HBt)));
   }

   // build the inverse for F 
   ModifiableLinearOp & invF = state.getModifiableOp("invF");
   if(invF==Teuchos::null)
      invF = buildInverse(*invVelFactory_,F);
   else 
      rebuildInverse(*invVelFactory_,F,invF);

   // build the approximate Schur complement: This is inefficient! FIXME
   ModifiableLinearOp & invS = state.getModifiableOp("invS");
   if(invS==Teuchos::null)
      invS = buildInverse(*invPrsFactory_,hatS);
   else
      rebuildInverse(*invPrsFactory_,hatS,invS);

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

      // Grab the sublist if we're using the block diagonal
      if(fInverseType_==BlkDiag)
	BlkDiagList_=pl.sublist("H options");      
   }
   if(pl.isParameter("Use Mass Scaling"))
      useMass_ = pl.get<bool>("Use Mass Scaling");

   Teko_DEBUG_MSG_BEGIN(5)
      DEBUG_STREAM << "SIMPLE Parameters: " << std::endl;
      DEBUG_STREAM << "   inv type    = \"" << invStr  << "\"" << std::endl;
      DEBUG_STREAM << "   inv v type  = \"" << invVStr << "\"" << std::endl;
      DEBUG_STREAM << "   inv p type  = \"" << invPStr << "\"" << std::endl;
      DEBUG_STREAM << "   alpha       = " << alpha_ << std::endl;
      DEBUG_STREAM << "   use mass    = " << useMass_ << std::endl;
      DEBUG_STREAM << "   vel scaling = " << getDiagonalName(fInverseType_) << std::endl;
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

   if(useMass_) {
      Teuchos::RCP<Teko::RequestHandler> rh = getRequestHandler();
      rh->preRequest<Teko::LinearOp>(Teko::RequestMesg("Velocity Mass Matrix"));
      Teko::LinearOp mass 
            = rh->request<Teko::LinearOp>(Teko::RequestMesg("Velocity Mass Matrix"));
      setMassMatrix(mass);
   }
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

   return result;
}

} // end namespace NS
} // end namespace Teko
