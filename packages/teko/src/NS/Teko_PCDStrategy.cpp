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

#include "Teko_PCDStrategy.hpp"

#include "Teuchos_TimeMonitor.hpp"
#include "Teko_Utilities.hpp"

namespace Teko {
namespace NS {

using Teuchos::TimeMonitor;

Teuchos::RCP<Teuchos::Time> PCDStrategy::initTimer_;
Teuchos::RCP<Teuchos::Time> PCDStrategy::invSTimer_;
Teuchos::RCP<Teuchos::Time> PCDStrategy::invFTimer_;
Teuchos::RCP<Teuchos::Time> PCDStrategy::opsTimer_;

void PCDStrategy::buildTimers()
{
   if(initTimer_==Teuchos::null)
      initTimer_ = TimeMonitor::getNewTimer("PCDStrategy::initializePrec");

   if(invSTimer_==Teuchos::null)
      invSTimer_ = TimeMonitor::getNewTimer("PCDStrategy::initializePrec invS");

   if(invFTimer_==Teuchos::null)
      invFTimer_ = TimeMonitor::getNewTimer("PCDStrategy::initializePrec invF");

   if(opsTimer_==Teuchos::null)
      opsTimer_ = TimeMonitor::getNewTimer("PCDStrategy::initializePrec buildOps");
}

PCDStrategy::PCDStrategy() : massInverseType_(Diagonal)
{ 
   buildTimers();
}

//! Constructor to set the inverse factories.
PCDStrategy::PCDStrategy(const Teuchos::RCP<InverseFactory> & invFA,
                                             const Teuchos::RCP<InverseFactory> & invS)
   : invFactoryF_(invFA), invFactoryS_(invS), massInverseType_(Diagonal)
{
   buildTimers();
}

/** returns the first (approximate) inverse of \f$A_{00}\f$ */
const Teko::LinearOp
PCDStrategy::getHatInvA00(const Teko::BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   initializeState(A,state);

   return state.getModifiableOp("invF");
}

/** returns the second (approximate) inverse of \f$A_{00}\f$ */
const Teko::LinearOp
PCDStrategy::getTildeInvA00(const Teko::BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   initializeState(A,state);

   return state.getModifiableOp("invF");
}

/** returns an (approximate) inverse of \f$S = -A_{11} + A_{10} \mbox{diag}(A_{00})^{-1} A_{01}\f$ */
const Teko::LinearOp
PCDStrategy::getInvS(const Teko::BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   initializeState(A,state);

   return state.getLinearOp("invS");
}

void PCDStrategy::initializeState(const Teko::BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   Teko_DEBUG_SCOPE("PCDStrategy::initializeState",10);

   std::string pcdStr      = getPCDString();
   std::string presLapStr  = getPressureLaplaceString();
   std::string presMassStr = getPressureMassString();

   // sanity checks
   TEUCHOS_ASSERT(state.getLinearOp(pcdStr)      !=Teuchos::null);
   TEUCHOS_ASSERT(state.getLinearOp(presLapStr)  !=Teuchos::null);
   TEUCHOS_ASSERT(state.getLinearOp(presMassStr) !=Teuchos::null);

   // no work to be done
   if(state.isInitialized())
      return;

   Teuchos::TimeMonitor timer(*initTimer_,true);

   // extract sub blocks
   LinearOp F  = Teko::getBlock(0,0,A);
   LinearOp Bt = Teko::getBlock(0,1,A);
   LinearOp B  = Teko::getBlock(1,0,A);
   LinearOp C  = Teko::getBlock(1,1,A);

   LinearOp  Qp = state.getLinearOp(presMassStr);

   // build the inverse Laplacian complement
   /////////////////////////////////////////////
   LinearOp iQp;
   if(massInverseType_==NotDiag) {
      ModifiableLinearOp & invMass = state.getModifiableOp("invMass");
      Teko_DEBUG_SCOPE("Building inv(Mass)",10);

      if(invMass==Teuchos::null)
         invMass = buildInverse(*invFactoryS_,Qp);
      else
         rebuildInverse(*invFactoryS_,Qp,invMass);

      iQp = invMass;
   }
   else {
      Teko_DEBUG_MSG("Building inverse mass of type \"" << Teko::getDiagonalName(massInverseType_) << "\"",10);
      iQp = getInvDiagonalOp(Qp,massInverseType_);
   }

   // build the inverse Laplacian complement
   /////////////////////////////////////////////
   ModifiableLinearOp & invLaplace = state.getModifiableOp("invLaplace");
   {
      Teko_DEBUG_SCOPE("Building inv(Laplace)",10);
      Teuchos::TimeMonitor timer(*invSTimer_,true);

      LinearOp laplace = state.getLinearOp(presLapStr);
      if(invLaplace==Teuchos::null)
         invLaplace = buildInverse(*invFactoryS_,laplace);
      else
         rebuildInverse(*invFactoryS_,laplace,invLaplace);
   }

   // build the inverse Schur complement
   /////////////////////////////////////////////
   {
      Teko_DEBUG_SCOPE("Building S",10);
      Teuchos::TimeMonitor timer(*opsTimer_,true);

      // build Schur-complement
      LinearOp pcd = state.getLinearOp(pcdStr);
      LinearOp invL = invLaplace;
      LinearOp invS = multiply(iQp,pcd,invL);

      state.addLinearOp("invS",invS);
   }

   // build inverse F
   /////////////////////////////////////////////
   {
      Teko_DEBUG_SCOPE("Building inv(F)",10);
      Teuchos::TimeMonitor timer(*invFTimer_,true);

      ModifiableLinearOp & invF = state.getModifiableOp("invF"); 
      if(invF==Teuchos::null)
         invF = buildInverse(*invFactoryF_,F);
      else
         rebuildInverse(*invFactoryF_,F,invF);
   }

   // mark state as initialized
   state.setInitialized(true);
}

/** \brief This function builds the internals of the state from a parameter list.
  *        
  * This function builds the internals of the LU 2x2 state
  * from a parameter list. Furthermore, it allows a 
  * developer to easily add a factory to the build system.
  *
  * \param[in] settings Parameter list to use as the internal settings
  * \param[in] invLib Inverse library to use for building inverse factory objects
  *
  * \note The default implementation does nothing.
  */
void PCDStrategy::initializeFromParameterList(const Teuchos::ParameterList & pl,
                                                        const InverseLibrary & invLib)
{
   Teko_DEBUG_SCOPE("PCDStrategy::initializeFromParameterList",10);

   std::string invStr="Amesos", invFStr="", invSStr="";
   massInverseType_ = Diagonal;

   // "parse" the parameter list
   if(pl.isParameter("Inverse Type"))
      invStr = pl.get<std::string>("Inverse Type");
   if(pl.isParameter("Inverse F Type"))
      invFStr = pl.get<std::string>("Inverse F Type");
   if(pl.isParameter("Inverse Laplace Type"))
      invSStr = pl.get<std::string>("Inverse Laplace Type");
   if(pl.isParameter("Inverse Mass Type")) {
      std::string massInverseStr = pl.get<std::string>("Inverse Mass Type");

      // build inverse types
      massInverseType_ = getDiagonalType(massInverseStr);
   }

   // set defaults as needed
   if(invFStr=="") invFStr = invStr;
   if(invSStr=="") invSStr = invStr;

   Teko_DEBUG_MSG_BEGIN(5)
      DEBUG_STREAM << "PCD Strategy Parameters: " << std::endl;
      DEBUG_STREAM << "   inv type   = \"" << invStr  << "\"" << std::endl;
      DEBUG_STREAM << "   inv F type = \"" << invFStr << "\"" << std::endl;
      DEBUG_STREAM << "   inv Laplace type = \"" << invSStr << "\"" << std::endl;
      DEBUG_STREAM << "   inv Mass type = \"" << Teko::getDiagonalName(massInverseType_) << "\"" << std::endl;
      DEBUG_STREAM << "PCD Strategy Parameter list: " << std::endl;
      pl.print(DEBUG_STREAM);
   Teko_DEBUG_MSG_END()

   // build velocity inverse factory
   invFactoryF_ = invLib.getInverseFactory(invFStr);
 
   if(invFStr==invSStr)
      invFactoryS_ = invFactoryF_;
   else
      invFactoryS_ = invLib.getInverseFactory(invSStr);
}

//! For assiting in construction of the preconditioner
Teuchos::RCP<Teuchos::ParameterList> PCDStrategy::getRequestedParameters() const 
{
   Teko_DEBUG_SCOPE("PCDStrategy::getRequestedParameters",10);
   Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());

   // grab parameters from F solver
   RCP<Teuchos::ParameterList> fList = invFactoryF_->getRequestedParameters();
   if(fList!=Teuchos::null) {
      Teuchos::ParameterList::ConstIterator itr;
      for(itr=fList->begin();itr!=fList->end();++itr)
         pl->setEntry(itr->first,itr->second);
   }

   // grab parameters from S solver
   RCP<Teuchos::ParameterList> sList = invFactoryS_->getRequestedParameters();
   if(sList!=Teuchos::null) {
      Teuchos::ParameterList::ConstIterator itr;
      for(itr=sList->begin();itr!=sList->end();++itr)
         pl->setEntry(itr->first,itr->second);
   }

   pl->set("Pressure Mass Operator", false,"Pressure mass matrix");
   pl->set("Pressure Laplace Operator", false,"Pressure Laplacian matrix");
   pl->set("PCD Operator", false,"Pressure convection-diffusion matrix");

   return pl;
}

//! For assiting in construction of the preconditioner
bool PCDStrategy::updateRequestedParameters(const Teuchos::ParameterList & pl) 
{
   Teko_DEBUG_SCOPE("PCDStrategy::updateRequestedParameters",10);
   bool result = true;
 
   // update requested parameters in solvers
   result &= invFactoryF_->updateRequestedParameters(pl);
   result &= invFactoryS_->updateRequestedParameters(pl);

   Teuchos::ParameterList hackList(pl);
   // get required operator acknowledgment...user must set these to true
   bool pmo = hackList.get<bool>("Pressure Mass Operator",false);
   bool plo = hackList.get<bool>("Pressure Laplace Operator",false);
   bool pcdo = hackList.get<bool>("PCD Operator", false);
   
   if(not pmo)  { Teko_DEBUG_MSG("User must acknowledge the use of the \"Pressure Mass Operator\"!",0); }
   if(not plo)  { Teko_DEBUG_MSG("User must acknowledge the use of the \"Pressure Laplace Operator\"!",0); }
   if(not pcdo) { Teko_DEBUG_MSG("User must acknowledge the use of the \"PCD Operator\"!",0); }

   result &= (pmo & plo & pcdo);

   return result;
}

} // end namespace NS
} // end namespace Teko
