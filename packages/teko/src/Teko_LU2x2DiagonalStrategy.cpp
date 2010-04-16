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

#include "Teko_LU2x2DiagonalStrategy.hpp"

#include "Teuchos_TimeMonitor.hpp"

namespace Teko {

using Teuchos::TimeMonitor;

Teuchos::RCP<Teuchos::Time> LU2x2DiagonalStrategy::initTimer_;
Teuchos::RCP<Teuchos::Time> LU2x2DiagonalStrategy::invSTimer_;
Teuchos::RCP<Teuchos::Time> LU2x2DiagonalStrategy::invA00Timer_;
Teuchos::RCP<Teuchos::Time> LU2x2DiagonalStrategy::opsTimer_;

void LU2x2DiagonalStrategy::buildTimers()
{
   if(initTimer_==Teuchos::null)
      initTimer_ = TimeMonitor::getNewTimer("LU2x2DiagonalStrategy::initializePrec");

   if(invSTimer_==Teuchos::null)
      invSTimer_ = TimeMonitor::getNewTimer("LU2x2DiagonalStrategy::initializePrec invS");

   if(invA00Timer_==Teuchos::null)
      invA00Timer_ = TimeMonitor::getNewTimer("LU2x2DiagonalStrategy::initializePrec invA00");

   if(opsTimer_==Teuchos::null)
      opsTimer_ = TimeMonitor::getNewTimer("LU2x2DiagonalStrategy::initializePrec buildOps");
}

LU2x2DiagonalStrategy::LU2x2DiagonalStrategy() 
   : a00InverseType_(Diagonal)
{ 
   buildTimers();
}

//! Constructor to set the inverse factories.
LU2x2DiagonalStrategy::LU2x2DiagonalStrategy(const Teuchos::RCP<InverseFactory> & invFA,
                                             const Teuchos::RCP<InverseFactory> & invS)
   : invFactoryA00_(invFA), invFactoryS_(invS), a00InverseType_(Diagonal)
{
   buildTimers();
}

/** returns the first (approximate) inverse of \f$A_{00}\f$ */
const Teko::LinearOp
LU2x2DiagonalStrategy::getHatInvA00(const Teko::BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   initializeState(A,state);

   return state.getModifiableOp("invA00");
}

/** returns the second (approximate) inverse of \f$A_{00}\f$ */
const Teko::LinearOp
LU2x2DiagonalStrategy::getTildeInvA00(const Teko::BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   initializeState(A,state);

   return state.getModifiableOp("invA00");
}

/** returns an (approximate) inverse of \f$S = -A_{11} + A_{10} \mbox{diag}(A_{00})^{-1} A_{01}\f$ */
const Teko::LinearOp
LU2x2DiagonalStrategy::getInvS(const Teko::BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   initializeState(A,state);

   return state.getModifiableOp("invS");
}

void LU2x2DiagonalStrategy::initializeState(const Teko::BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   Teko_DEBUG_SCOPE("LU2x2DiagonalStrategy::initializeState",10);

   // no work to be done
   if(state.isInitialized())
      return;

   Teuchos::TimeMonitor timer(*initTimer_,true);

   // extract sub blocks
   LinearOp A00 = Teko::getBlock(0,0,A);
   LinearOp A01 = Teko::getBlock(0,1,A);
   LinearOp A10 = Teko::getBlock(1,0,A);
   LinearOp A11 = Teko::getBlock(1,1,A);

   // build the Schur complement
   /////////////////////////////////////////////
   ModifiableLinearOp & S = state.getModifiableOp("S");
   {
      Teko_DEBUG_SCOPE("Building S",5);
      Teuchos::TimeMonitor timer(*opsTimer_,true);

      LinearOp diagA00 = getInvDiagonalOp(A00,a00InverseType_);
   
      // build Schur-complement
      ModifiableLinearOp & triple = state.getModifiableOp("triple");
      triple = explicitMultiply(A10,diagA00,A01,triple);
      S = explicitAdd(scale(-1.0,A11),triple,S);
   }

   // build inverse S
   /////////////////////////////////////////////
   {
      Teko_DEBUG_SCOPE("Building inverse(S)",5);
      Teuchos::TimeMonitor timer(*invSTimer_,true);

      ModifiableLinearOp & invS = state.getModifiableOp("invS"); 
      if(invS==Teuchos::null)
         invS = buildInverse(*invFactoryS_,S);
      else
         rebuildInverse(*invFactoryS_,S,invS);
   }

   // build inverse A00
   /////////////////////////////////////////////
   {
      Teko_DEBUG_SCOPE("Building inverse(A00)",5);
      Teuchos::TimeMonitor timer(*invA00Timer_,true);

      ModifiableLinearOp & invA00 = state.getModifiableOp("invA00"); 
      *getOutputStream() << "(LU2x2) invA00 pointer = " << invA00 << std::endl;
      if(invA00==Teuchos::null)
         invA00 = buildInverse(*invFactoryA00_,A00);
      else
         rebuildInverse(*invFactoryA00_,A00,invA00);
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
void LU2x2DiagonalStrategy::initializeFromParameterList(const Teuchos::ParameterList & pl,
                                                        const InverseLibrary & invLib)
{
   Teko_DEBUG_SCOPE("LU2x2DiagonalStrategy::initializeFromParameterList",10);

   std::string invStr="Amesos", invA00Str="", invSStr="";

   // "parse" the parameter list
   if(pl.isParameter("Inverse Type"))
      invStr = pl.get<std::string>("Inverse Type");
   if(pl.isParameter("Inverse A00 Type"))
      invA00Str = pl.get<std::string>("Inverse A00 Type");
   if(pl.isParameter("Inverse Schur Type"))
      invSStr = pl.get<std::string>("Inverse Schur Type");
   if(pl.isParameter("Diagonal Type")) {
      std::string massInverseStr = pl.get<std::string>("Diagonal Type");

      // build inverse types
      a00InverseType_ = getDiagonalType(massInverseStr);
   }

   // set defaults as needed
   if(invA00Str=="") invA00Str = invStr;
   if(invSStr=="") invSStr = invStr;

   Teko_DEBUG_MSG_BEGIN(5)
      DEBUG_STREAM << "LU2x2 Diagonal Strategy Parameters: " << std::endl;
      DEBUG_STREAM << "   inv type   = \"" << invStr  << "\"" << std::endl;
      DEBUG_STREAM << "   inv A00 type = \"" << invA00Str << "\"" << std::endl;
      DEBUG_STREAM << "   inv S type = \"" << invSStr << "\"" << std::endl;
      DEBUG_STREAM << "LU2x2 Diagonal Strategy Parameter list: " << std::endl;
      pl.print(DEBUG_STREAM);
   Teko_DEBUG_MSG_END()

   // build velocity inverse factory
   invFactoryA00_ = invLib.getInverseFactory(invA00Str);
 
   if(invA00Str==invSStr)
      invFactoryS_ = invFactoryA00_;
   else
      invFactoryS_ = invLib.getInverseFactory(invSStr);
}
} // end namespace Teko
