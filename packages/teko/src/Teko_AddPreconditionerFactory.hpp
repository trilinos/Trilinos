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

#ifndef __Teko_AddPreconditionerFactory_hpp__
#define __Teko_AddPreconditionerFactory_hpp__

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_Utilities.hpp"

namespace Teko {

/** Preconditioning factories must supply a 'State' class which
  * is where data specific to the preconditioner construction 
  * is stored. The constructor will be invoked elsewhere.
  */
class AddPrecondState : public Teko::BlockPreconditionerState {
public:
   AddPrecondState() {}

   Teuchos::RCP<BlockPreconditionerState> StateOne_;
   Teuchos::RCP<BlockPreconditionerState> StateTwo_;
};

/** Declaration of preconditioner factory that creates
  * a preconditioner which is the sum (additive) of two
  * other preconditioners.
  */
class AddPreconditionerFactory 
   : public Teko::BlockPreconditionerFactory {
public:
   //! Constructor
   AddPreconditionerFactory(const Teuchos::RCP<const Teko::BlockPreconditionerFactory> & FirstFactory,
                            const Teuchos::RCP<const Teko::BlockPreconditionerFactory> & SecondFactory);

   AddPreconditionerFactory();

   //! Function inherited from Teko::BlockPreconditionerFactory
   Teko::LinearOp buildPreconditionerOperator(Teko::BlockedLinearOp & blo,
                                            Teko::BlockPreconditionerState & state) const;
    
   //! Build the AddPrecondState object
   virtual Teuchos::RCP<Teko::PreconditionerState> buildPreconditionerState() const;

protected:
   // class members
   Teuchos::RCP<const Teko::BlockPreconditionerFactory> FirstFactory_;
   Teuchos::RCP<const Teko::BlockPreconditionerFactory> SecondFactory_;

   //! Initialize from a parameter list
   virtual void initializeFromParameterList(const Teuchos::ParameterList & pl);
};

} // end namespace Teko

#endif
