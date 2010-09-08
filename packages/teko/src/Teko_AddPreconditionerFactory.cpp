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

#include "Teko_AddPreconditionerFactory.hpp"

namespace Teko {

using Teuchos::RCP;

AddPreconditionerFactory::AddPreconditionerFactory(
                            const RCP<const BlockPreconditionerFactory> & FirstFactory,
                            const RCP<const BlockPreconditionerFactory> & SecondFactory)
   : FirstFactory_(FirstFactory), SecondFactory_(SecondFactory)
{}

AddPreconditionerFactory::AddPreconditionerFactory()
{}

//! Build the AddPrecondState object
RCP<PreconditionerState> AddPreconditionerFactory::buildPreconditionerState() const
{ 
   AddPrecondState*   mystate = new AddPrecondState(); 
   mystate->StateOne_ = Teuchos::rcp_dynamic_cast<BlockPreconditionerState>(FirstFactory_->buildPreconditionerState());
   mystate->StateTwo_ = Teuchos::rcp_dynamic_cast<BlockPreconditionerState>(SecondFactory_->buildPreconditionerState());
   return rcp(mystate);
}

// Use the factory to build the preconditioner (this is where the work goes)
LinearOp AddPreconditionerFactory 
   ::buildPreconditionerOperator(BlockedLinearOp & blockOp,
                                 BlockPreconditionerState & state) const
{
   // The main tricky thing here is that we have to take the 'state' object
   // associated with AddPreconditionerFactory(), pull out the states for
   // the individual preconditioners, and pass these on to 
   // buildPreconditionerOperator() for each subpreconditioner.
   
   AddPrecondState *MyState = dynamic_cast<AddPrecondState *> (&state);
   TEUCHOS_ASSERT(MyState != 0);

   LinearOp M1 = FirstFactory_->buildPreconditionerOperator(blockOp, *MyState->StateOne_);
   LinearOp M2 = SecondFactory_->buildPreconditionerOperator(blockOp, *MyState->StateTwo_);

   LinearOp invA = add(M1, M2);

   // return fully constructed preconditioner
   return invA;
}

//! Initialize from a parameter list
void AddPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & pl)
{
   RCP<const InverseLibrary> invLib = getInverseLibrary();

   // get string specifying inverse
   std::string aStr="", bStr="";

   // "parse" the parameter list
   aStr = pl.get<std::string>("Preconditioner A");
   bStr = pl.get<std::string>("Preconditioner B");

   RCP<const Teuchos::ParameterList> aSettings = invLib->getParameterList(aStr);
   RCP<const Teuchos::ParameterList> bSettings = invLib->getParameterList(bStr);

   // build preconditioner from the parameters
   std::string aType = aSettings->get<std::string>("Preconditioner Type");
   RCP<Teko::PreconditionerFactory> precA
         = Teko::PreconditionerFactory::buildPreconditionerFactory(aType,aSettings->sublist("Preconditioner Settings"),invLib);

   // build preconditioner from the parameters
   std::string bType = bSettings->get<std::string>("Preconditioner Type");
   RCP<Teko::PreconditionerFactory> precB
         = Teko::PreconditionerFactory::buildPreconditionerFactory(bType,bSettings->sublist("Preconditioner Settings"),invLib);

   // set precondtioners
   FirstFactory_ = Teuchos::rcp_dynamic_cast<const Teko::BlockPreconditionerFactory>(precA);
   SecondFactory_ = Teuchos::rcp_dynamic_cast<const Teko::BlockPreconditionerFactory>(precB);
}

} // end namespace Teko
