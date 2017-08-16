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

#include "Teko_BlockImplicitLinearOp.hpp"

namespace Teko {

using Teuchos::rcpFromRef;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcp_const_cast;
using Teuchos::RCP;

using Thyra::ProductMultiVectorBase;

void 
BlockImplicitLinearOp::
implicitApply(const Thyra::EOpTransp M_trans,
              const BlockedMultiVector & x, BlockedMultiVector & y,
              const double alpha, const double beta) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(M_trans!=Thyra::NOTRANS, std::runtime_error,
     "Linear operators of inherited type BlockImplicitLinearOp "
     "cannot handle conjugation (yet!)");

   // call apply
   implicitApply(x,y,alpha,beta);
}

bool BlockImplicitLinearOp::opSupportedImpl(const Thyra::EOpTransp M_trans) const
{
  return (M_trans == Thyra::NOTRANS);
}

void BlockImplicitLinearOp::applyImpl(
  const Thyra::EOpTransp M_trans,
  const Thyra::MultiVectorBase<double> & x,
  const Teuchos::Ptr<Thyra::MultiVectorBase<double> > & y,
  const double alpha,
  const double beta
  ) const
{
   // cast source vector
   RCP<const ProductMultiVectorBase<double> > src =
     rcp_dynamic_cast<const ProductMultiVectorBase<double> >(rcpFromRef(x));
   BlockedMultiVector srcX = rcp_const_cast<ProductMultiVectorBase<double> >(src);

   // cast destination vector
   BlockedMultiVector destY =
     rcp_dynamic_cast<ProductMultiVectorBase<double> >(rcpFromPtr(y));

   // call apply
   implicitApply(M_trans,srcX,destY,alpha,beta);
}
 
} // end namespace Teko
