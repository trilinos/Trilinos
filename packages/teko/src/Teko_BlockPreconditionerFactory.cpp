/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
//  
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// Export of this program may require a license from the United States
// Government.
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

#include "Teko_BlockPreconditionerFactory.hpp"

#include "Teko_Preconditioner.hpp"
#include "Teko_InverseLibrary.hpp"

#include "Thyra_DefaultPreconditioner.hpp"

using namespace Thyra;

namespace Teko {

/////////////////////////////////////////////////////

LinearOp BlockPreconditionerFactory::buildPreconditionerOperator(LinearOp & lo,PreconditionerState & state) const
{
   // get the blocked linear operator
   RCP<LinearOpBase<double> > loA = Teuchos::rcp_const_cast<Thyra::LinearOpBase<double> >(lo);
   BlockedLinearOp A = Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<double> >(loA);

   state.setInitialized(false);

   return buildPreconditionerOperator(A,dynamic_cast<BlockPreconditionerState &>(state));
}

//! is this operator compatiable with the preconditioner factory?
bool BlockPreconditionerFactory::isCompatible(const Thyra::LinearOpSourceBase<double> &fwdOpSrc) const
{
   RCP<const Thyra::PhysicallyBlockedLinearOpBase<double> > A 
         = Teuchos::rcp_dynamic_cast<const Thyra::PhysicallyBlockedLinearOpBase<double> >(fwdOpSrc.getOp());
   return A!=Teuchos::null;
}

} // end namespace Teko
