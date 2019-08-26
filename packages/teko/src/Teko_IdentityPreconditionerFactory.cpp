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

#include "Teko_IdentityPreconditionerFactory.hpp"

#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"

using Teuchos::rcp;

namespace Teko {

/** Build a Identity preconditioner factory from a parameter list 
  */
IdentityPreconditionerFactory::IdentityPreconditionerFactory()
  : scaling_(1.0)
{ }

LinearOp IdentityPreconditionerFactory::buildPreconditionerOperator(LinearOp & lo,PreconditionerState & /* state */) const
{
   return Thyra::scale(scaling_,Thyra::identity(rangeSpace(lo)));
}

//! Initialize from a parameter list
void IdentityPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & pl)
{
   Teko_DEBUG_SCOPE("IdentityPreconditionerFactory::initializeFromParameterList",10);
   Teko_DEBUG_MSG_BEGIN(9);
      DEBUG_STREAM << "Parameter list: " << std::endl;
      pl.print(DEBUG_STREAM);
   Teko_DEBUG_MSG_END();

   // get string specifying default inverse
   std::string scaleStr ="Scaling"; 
   if(pl.isParameter(scaleStr))
      scaling_ = pl.get<double>(scaleStr);
}

} // end namspace Teko
