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

#include "Teko_SolveInverseFactory.hpp"

// Thyra includes
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"

// Stratimikos includes
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

// Teko includes
#include "Teko_Utilities.hpp"
#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_Preconditioner.hpp"
#include "Teko_PreconditionerLinearOp.hpp"

using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;

namespace Teko {

/** \brief Constructor that takes a Thyra solve factory and 
  *        makes it look like an InverseFactory
  *
  * Constructor that takes a Thyra solve factory and 
  * makes it look like an InverseFactory.
  * 
  * \param[in] lowsFactory Thyra LineaerOpWithSolveFactoryBase used for building 
  *                        the inverse.
  */
SolveInverseFactory::SolveInverseFactory(const Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > & lowsFactory)
   : lowsFactory_(lowsFactory)
{ }

//! Copy constructor
SolveInverseFactory::SolveInverseFactory(const SolveInverseFactory & siFactory)
   : lowsFactory_(siFactory.lowsFactory_)
{ }

/** \brief Build an inverse operator
  *
  * Build the inverse operator using this factory.
  *
  * \param[in] linearOp Linear operator needing to be inverted.
  *
  * \returns New linear operator that functions as the inverse
  *          of <code>linearOp</code>.
  */
InverseLinearOp SolveInverseFactory::buildInverse(const LinearOp & linearOp) const
{
   Teko_DEBUG_SCOPE("SolveInverseFactory::buildInverse(linearOp)",10);

   // build and initialize inverse linear op with solve
   Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > invLOWS = lowsFactory_->createOp();
   lowsFactory_->initializeOp(Thyra::defaultLinearOpSource(linearOp),&*invLOWS,Thyra::SUPPORT_SOLVE_FORWARD_ONLY);
   
   return Thyra::nonconstInverse<double>(invLOWS);
}

/** \brief Build a preconditioned inverse operator
  * 
  * Build the inverse operator using this factory and a user specified
  * preconditioning operator. The default behavior is to call buildInverse
  * ignoring the preconditioner. 
  *
  * \param[in] linearOp Linear operator needing to be inverted.
  * \param[in] precOp Preconditioning operator
  *
  * \returns New linear operator that functions as the inverse
  *          of <code>linearOp</code>.
  */
InverseLinearOp SolveInverseFactory::buildInverse(const LinearOp & linearOp,const LinearOp & precOp) const
{ 
   Teko_DEBUG_SCOPE("SolveInverseFactory::buildInverse(linearOp,precOp)",10);

   // build and initialize inverse linear op with solve
   Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > invLOWS = lowsFactory_->createOp();
   lowsFactory_->initializePreconditionedOp(Thyra::defaultLinearOpSource(linearOp),
                                                  Thyra::unspecifiedPrec(precOp),
                                                  &*invLOWS,Thyra::SUPPORT_SOLVE_FORWARD_ONLY);
   
   return Thyra::nonconstInverse<double>(invLOWS);
}

/** \brief Pass in an already constructed inverse operator. Update
  *        the inverse operator based on the new source operator.
  *
  * Pass in an already constructed inverse operator. Update
  * the inverse operator based on the new source operator.
  *
  * \params[in]     source Source operator to be inverted.
  * \params[in,out] dest   Pre constructed inverse operator to be
  *                        rebuilt using the <code>source</code>
  *                        object.
  */
void SolveInverseFactory::rebuildInverse(const LinearOp & source,InverseLinearOp & dest) const
{
   RCP<Thyra::DefaultInverseLinearOp<double> > invDest = rcp_dynamic_cast<Thyra::DefaultInverseLinearOp<double> >(dest);
   RCP<Thyra::LinearOpWithSolveBase<double> > lows = invDest->getNonconstLows();

   // This stems from confusion of if the linear op with solve initializeAndResuseOp actually rebuilds
   // the precondtioner. It seems not to and thus we have a fairly substantial problem.
   // lowsFactory_->initializeAndReuseOp(Thyra::defaultLinearOpSource(source),&*lows);
   lowsFactory_->initializeOp(Thyra::defaultLinearOpSource(source),&*lows);
}

/** \brief Pass in an already constructed inverse operator. Update
  *        the inverse operator based on the new source operator.
  *
  * Pass in an already constructed inverse operator. Update
  * the inverse operator based on the new source operator.
  *
  * \param[in]     source Source operator to be inverted.
  * \param[in]     precOp Preconditioning operator
  * \param[in,out] dest   Pre constructed inverse operator to be
  *                        rebuilt using the <code>source</code>
  *                        object.
  */
void SolveInverseFactory::rebuildInverse(const LinearOp & source,const LinearOp & precOp,InverseLinearOp & dest) const
{
   RCP<Thyra::DefaultInverseLinearOp<double> > invDest = rcp_dynamic_cast<Thyra::DefaultInverseLinearOp<double> >(dest);
   RCP<Thyra::LinearOpWithSolveBase<double> > lows = invDest->getNonconstLows();

   // lowsFactory_->initializeAndReuseOp(Thyra::defaultLinearOpSource(source),&*lows);
   lowsFactory_->initializePreconditionedOp(Thyra::defaultLinearOpSource(source),
                                                  Thyra::unspecifiedPrec(precOp),
                                                  &*lows,Thyra::SUPPORT_SOLVE_FORWARD_ONLY);
}

/** \brief A function that permits inspection of the parameters used to create
  *        this object.
  *
  * A function that permits inspection of the parameters used to create this
  * object. Useful for determining defaults and settings used.
  *
  * \returns A list used to parameterize this object.
  */
Teuchos::RCP<const Teuchos::ParameterList> SolveInverseFactory::getParameterList() const
{ 
   return lowsFactory_->getParameterList(); 
}

} // end namespace Teko
