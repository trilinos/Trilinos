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

// Teko includes
#include "Teko_IterativePreconditionerFactory.hpp"

#include "Teko_PreconditionerInverseFactory.hpp"

namespace Teko {

//! Default constructor, for use with the AutoClone class.
IterativePreconditionerFactory::IterativePreconditionerFactory()
   : correctionNum_(0), precFactory_(Teuchos::null)
{ }

/** Construct a preconditioner factory that applies a specified
  * preconditioner, a fixed number of times.
  */
IterativePreconditionerFactory::IterativePreconditionerFactory(unsigned int correctionNum,
                            const Teuchos::RCP<Teko::InverseFactory> & precFactory)
   : correctionNum_(correctionNum), precFactory_(precFactory)
{ }

/** Construct a preconditioner factory that applies a specified
  * preconditioner, a fixed number of times.
  */
IterativePreconditionerFactory::IterativePreconditionerFactory(unsigned int correctionNum,
                            const Teuchos::RCP<Teko::PreconditionerFactory> & precFactory)
   : correctionNum_(correctionNum)
{
   precFactory_ = Teuchos::rcp(new Teko::PreconditionerInverseFactory(precFactory,precFactory->getRequestHandler()));
}


/** \brief Function that is called to build the preconditioner
  *        for the linear operator that is passed in.
  */
LinearOp IterativePreconditionerFactory::buildPreconditionerOperator(LinearOp & lo,PreconditionerState & state) const
{
   TEST_FOR_EXCEPTION(precFactory_==Teuchos::null,std::runtime_error,
                      "ERROR: Teko::IterativePreconditionerFactory::buildPreconditionerOperator requires that a "
                   << "preconditioner factory has been set. Currently it is null!");

   // build user specified preconditioner
   ModifiableLinearOp & invP = state.getModifiableOp("prec");
   if(invP==Teuchos::null)
      invP = Teko::buildInverse(*precFactory_,lo);
   else
      Teko::rebuildInverse(*precFactory_,lo,invP);

   // // no repititions are required
   // if(correctionNum_==0) return invP;

   // now build correction operator
   LinearOp I = Thyra::identity(lo->range(),"I");
   LinearOp AiP = multiply(lo,invP.getConst(),"AiP");
   LinearOp correction = add(I,scale(-1.0,AiP)); // I - A * iPA

   LinearOp resMap = I; // will map r_n to r_{n+1}
   for(unsigned int i=0;i<correctionNum_;i++) 
      resMap = add(I,multiply(resMap,correction)); // resMap = I + resMap*(I-A*iP)
   
   // iP = (I-A*iP)^{correctionNum}
   return multiply(invP.getConst(),resMap);
}

/** \brief This function builds the internals of the preconditioner factory
  *        from a parameter list.
  */
void IterativePreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & settings)
{

   correctionNum_ = 1;
   if(settings.isParameter("Iteration Count"))
      correctionNum_ = settings.get<int>("Iteration Count");
 
   TEST_FOR_EXCEPTION(not settings.isParameter("Preconditioner Type"),std::runtime_error,
                      "Parameter \"Preconditioner Type\" is required by a Teko::IterativePreconditionerFactory");
      
   // grab library and preconditioner name
   Teuchos::RCP<const InverseLibrary> il = getInverseLibrary();
   std::string precName = settings.get<std::string>("Preconditioner Type");

   // build preconditioner factory
   precFactory_ = il->getInverseFactory(precName);
   TEST_FOR_EXCEPTION(precFactory_==Teuchos::null,std::runtime_error,
                      "ERROR: \"Preconditioner Type\" = " << precName 
                   << " could not be found");
}

/** \brief Request the additional parameters this preconditioner factory
  *        needs. 
  */
Teuchos::RCP<Teuchos::ParameterList> IterativePreconditionerFactory::getRequestedParameters() const
{
   TEST_FOR_EXCEPTION(precFactory_==Teuchos::null,std::runtime_error,
                      "ERROR: Teko::IterativePreconditionerFactory::getRequestedParameters requires that a "
                   << "preconditioner factory has been set. Currently it is null!");

   return precFactory_->getRequestedParameters();
}

/** \brief Update this object with the fields from a parameter list.
  */
bool IterativePreconditionerFactory::updateRequestedParameters(const Teuchos::ParameterList & pl)
{
   TEST_FOR_EXCEPTION(precFactory_==Teuchos::null,std::runtime_error,
                      "ERROR: Teko::IterativePreconditionerFactory::updateRequestedParameters requires that a "
                   << "preconditioner factory has been set. Currently it is null!");

   return precFactory_->updateRequestedParameters(pl);
}

} // end namespace Teko
