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

#ifndef __Teko_Preconditioner_hpp__
#define __Teko_Preconditioner_hpp__

// Thyra includes
#include "Thyra_DefaultPreconditioner.hpp"

// Teko includes
#include "Teko_Utilities.hpp"
#include "Teko_PreconditionerState.hpp"

namespace Teko {

using Thyra::LinearOpBase;
using Thyra::DefaultPreconditioner;
using Teuchos::RCP;

/** \brief An extension of the <code>Thyra::DefaultPreconditioner</code>
  *        class with some specializations useful for use within Teko.
  *
  * An extension of the <code>Thyra::DefaultPreconditioner</code>
  * class with some specializations useful for use within Teko. This includes
  * having facilities to store the source vector and the state object.
  */
class Preconditioner : public DefaultPreconditioner<double> {
public:
   //! \name Constructors based from Thyra::DefaultPreconditioner
   //@{
   Preconditioner()
      : DefaultPreconditioner<double>() {}
   Preconditioner(const RCP<LinearOpBase<double> > & leftPrecOp,
                  const RCP<LinearOpBase<double> > & rightPrecOp)
      : DefaultPreconditioner<double>(leftPrecOp,rightPrecOp) {}
   Preconditioner(const RCP<const LinearOpBase<double> > & leftPrecOp,
                  const RCP<const LinearOpBase<double> > & rightPrecOp)
      : DefaultPreconditioner<double>(leftPrecOp,rightPrecOp) {}
   Preconditioner(const RCP<LinearOpBase<double> > & unspecifiedPrecOp)
      : DefaultPreconditioner<double>(unspecifiedPrecOp) {}
   Preconditioner(const RCP<const LinearOpBase<double> > & unspecifiedPrecOp)
      : DefaultPreconditioner<double>(unspecifiedPrecOp) {}
   //@}

   /** Set the vector associated with this operator (think nonlinear system)
     *
     * \param[in] srcVec The source vector associated with this preconditioner.
     */
   virtual void setSourceVector(const RCP<Thyra::MultiVectorBase<double> > & srcVec)
   { if(srcVec!=Teuchos::null) state_->setSourceVector(srcVec);
     else                      state_->setSourceVector(Teuchos::null); }

   /** Set the state object associated with this preconditioner
     *
     * \param[in] state The state object to use
     */
   virtual void setStateObject(const RCP<PreconditionerState> & state)
   { state_ = state; }

   /** Get the state object associated with this preconditioner
     *
     * \returns The state object used by this preconditioner
     */
   virtual const RCP<PreconditionerState> getStateObject()
   { return state_; }

   /** Get the state object associated with this preconditioner
     *
     * \returns The state object used by this preconditioner
     */
   virtual const RCP<const PreconditionerState> getStateObject() const
   { return state_; }

   /** Merge a state object into the one used by this Preconditioner
     * (Merge in this case is the same as that used by PreconditionerState::merge).
     *
     * \param[in] state A state object to merge into the internal state object.
     */
   virtual void mergeStateObject(const PreconditionerState & state)
   { state_->merge(state); }

protected:
   //! User defined state for this preconditioner
   RCP<PreconditionerState> state_;
};

} // end namespace Teko

#endif
