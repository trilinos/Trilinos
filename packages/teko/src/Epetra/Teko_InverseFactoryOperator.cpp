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

#include "Teko_InverseFactoryOperator.hpp"

// Teko includes
#include "Teko_BasicMappingStrategy.hpp"

// Thyra includes
#include "Thyra_EpetraLinearOp.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::rcp_dynamic_cast;

namespace Teko {
namespace Epetra {

/** \brief Constructor that takes the InverseFactory that will
  *        build the operator.
  *
  * Constructor that takes the InverseFactory that will
  * build the operator.
  */
InverseFactoryOperator::InverseFactoryOperator(const Teuchos::RCP<const InverseFactory> & ifp)
   : inverseFactory_(ifp), firstBuildComplete_(false), setConstFwdOp_(true)
{
}

/** \brief Build the underlying data structure for the inverse operator.
  *        
  * Build the underlying data structure for the inverse operator. This
  * permits the manipulation of the state object for an inverse operator.
  *
  * \param[in] clearOld If true any previously constructed
  *                     operator will be wiped out and
  *                     a new one created. If false, anoperator 
  *                     will be created only if the current one is
  *                     empty (i.e. <code>initPreconditioner</code>
  *                     had not been called).
  */
void InverseFactoryOperator::initInverse(bool clearOld)
{
   if(not clearOld)
      return;
   invOperator_ = Teuchos::null;
}

/** \brief Build this inverse operator from an Epetra_Operator 
  * passed in to this object.
  *
  * Build this inverse opeerator from an Epetra_Operator 
  * passed in to this object. If this Epetra_Operator
  * is an EpetraOperatorWrapper object then the block Thyra components
  * are extracted.
  *
  * \param[in] A The Epetra source operator.
  * \param[in] clear If true, than any previous state saved by the operator 
  *                  is discarded.
  */
void InverseFactoryOperator::buildInverseOperator(const Teuchos::RCP<const Epetra_Operator> & A,bool clear)
{
   Teko_DEBUG_SCOPE("InverseFactoryOperator::buildInverseOperator",10);

   // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
   RCP<const Thyra::LinearOpBase<double> > thyraA = extractLinearOp(A);

   // set the mapping strategy
   SetMapStrategy(rcp(new InverseMappingStrategy(extractMappingStrategy(A))));

   initInverse(clear);

   // actually build the inverse operator
   invOperator_ = Teko::buildInverse(*inverseFactory_,thyraA); 

   SetOperator(invOperator_,false);

   firstBuildComplete_ = true;

   if(setConstFwdOp_)
      fwdOp_ = A;

   setConstFwdOp_ = true;

   TEUCHOS_ASSERT(invOperator_!=Teuchos::null);
   TEUCHOS_ASSERT(getForwardOp()!=Teuchos::null);
   TEUCHOS_ASSERT(getThyraOp()!=Teuchos::null);
   TEUCHOS_ASSERT(getMapStrategy()!=Teuchos::null);
   TEUCHOS_ASSERT(firstBuildComplete_==true);
}

void InverseFactoryOperator::buildInverseOperator(const Teuchos::RCP<Epetra_Operator> & A,bool clear)
{
   setConstFwdOp_ = false;

   fwdOp_.initialize(A);

   buildInverseOperator(A.getConst());

   TEUCHOS_ASSERT(setConstFwdOp_==true);
}

/** \brief Rebuild this inverse from an Epetra_Operator passed
  * in this to object. 
  *
  * Rebuild this inverse from an Epetra_Operator passed
  * in this to object.  If <code>buildInverseOperator</code> has not been called
  * the inverse operator will be built instead. Otherwise efforts are taken
  * to only rebuild what is neccessary. Also, that this Epetra_Operator
  * may be an EpetraOperatorWrapper object, so the block Thyra components
  * can be extracted.
  *
  * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
  */
void InverseFactoryOperator::rebuildInverseOperator(const Teuchos::RCP<const Epetra_Operator> & A)
{
   Teko_DEBUG_SCOPE("InverseFactoryOperator::rebuildPreconditioner",10);

   // if the inverse hasn't been built yet, rebuild from scratch
   if(not firstBuildComplete_) {
      buildInverseOperator(A,false);
      return;
   }

   RCP<const Thyra::LinearOpBase<double> > thyraA = extractLinearOp(A);
   Teko::rebuildInverse(*inverseFactory_,thyraA,invOperator_); 

   if(setConstFwdOp_)
      fwdOp_.initialize(A);

   SetOperator(invOperator_,false);

   setConstFwdOp_ = true;

   TEUCHOS_ASSERT(getForwardOp()!=Teuchos::null);
   TEUCHOS_ASSERT(invOperator_!=Teuchos::null);
   TEUCHOS_ASSERT(getThyraOp()!=Teuchos::null);
   TEUCHOS_ASSERT(firstBuildComplete_==true);
}

void InverseFactoryOperator::rebuildInverseOperator(const Teuchos::RCP<Epetra_Operator> & A)
{
   setConstFwdOp_ = false;

   fwdOp_.initialize(A);
 
   // build from constant epetra operator
   rebuildInverseOperator(A.getConst());

   TEUCHOS_ASSERT(setConstFwdOp_==true);
}

Teuchos::RCP<const Thyra::LinearOpBase<double> > InverseFactoryOperator::extractLinearOp(const Teuchos::RCP<const Epetra_Operator> & A) const
{
   // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
   const RCP<const EpetraOperatorWrapper> & eow = rcp_dynamic_cast<const EpetraOperatorWrapper>(A);
  
   // if it is an EpetraOperatorWrapper, then get the Thyra operator
   if(eow!=Teuchos::null)
      return eow->getThyraOp(); 

   // otherwise wrap it up as a thyra operator 
   return Thyra::epetraLinearOp(A);
}

Teuchos::RCP<const MappingStrategy> InverseFactoryOperator::extractMappingStrategy(const Teuchos::RCP<const Epetra_Operator> & A) const
{
   // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
   const RCP<const EpetraOperatorWrapper> & eow = rcp_dynamic_cast<const EpetraOperatorWrapper>(A);
  
   // if it is an EpetraOperatorWrapper, then get the Thyra operator
   if(eow!=Teuchos::null)
      return eow->getMapStrategy(); 

   // otherwise wrap it up as a thyra operator 
   RCP<const Epetra_Map> range = rcpFromRef(A->OperatorRangeMap());
   RCP<const Epetra_Map> domain = rcpFromRef(A->OperatorDomainMap());
   return rcp(new BasicMappingStrategy(range,domain,A->Comm()));
}

} // end namespace Epetra
} // end namespace Teko
