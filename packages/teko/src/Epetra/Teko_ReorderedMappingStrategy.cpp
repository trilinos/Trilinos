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

#include "Epetra/Teko_ReorderedMappingStrategy.hpp"

#include "Teko_BlockedReordering.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcpFromRef;

namespace Teko {
namespace Epetra {

ReorderedMappingStrategy::ReorderedMappingStrategy(const BlockReorderManager & brm,const Teuchos::RCP<const MappingStrategy> & map)
   : reorderManager_(brm), mapStrategy_(map)
{
   rangeMap_ = mapStrategy_->rangeMap();
   domainMap_ = mapStrategy_->domainMap();
}

void ReorderedMappingStrategy::copyEpetraIntoThyra(const Epetra_MultiVector& X,
                                                 const Teuchos::Ptr<Thyra::MultiVectorBase<double> > & thyra_X,
                                                 const Teko::Epetra::EpetraOperatorWrapper & eow) const
{
   using Teuchos::ptr_const_cast;
   using Teuchos::rcp_const_cast;

   // first flatten the vector: notice this just works on the block structure
   RCP<Thyra::ProductMultiVectorBase<double> > prod_X = rcp_dynamic_cast<Thyra::ProductMultiVectorBase<double> >(rcpFromRef(*thyra_X));
   RCP<Thyra::MultiVectorBase<double> > flat_X = buildFlatMultiVector(reorderManager_,prod_X);

   // now use the underlying mapping strategy to copy the flat vector 
   mapStrategy_->copyEpetraIntoThyra(X,flat_X.ptr(),eow);
}

void ReorderedMappingStrategy::copyThyraIntoEpetra(const RCP<const Thyra::MultiVectorBase<double> > & thyra_Y,
                                                 Epetra_MultiVector& Y,
                                                 const Teko::Epetra::EpetraOperatorWrapper & eow) const
{
   // first flatten the vector: notice this just works on the block structure
   RCP<const Thyra::ProductMultiVectorBase<double> > prod_Y = rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<double> >(rcpFromRef(*thyra_Y));
   RCP<const Thyra::MultiVectorBase<double> > flat_Y = buildFlatMultiVector(reorderManager_,prod_Y);

   // now use the underlying mapping strategy to copy the flat vector 
   mapStrategy_->copyThyraIntoEpetra(flat_Y,Y,eow);
}

} // end namespace Epetra
} // end namespace Teko
