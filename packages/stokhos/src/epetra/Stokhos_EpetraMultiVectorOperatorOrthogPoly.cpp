// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Stokhos_EpetraMultiVectorOperatorOrthogPoly.hpp"

Stokhos::EpetraMultiVectorOperatorOrthogPoly::
EpetraMultiVectorOperatorOrthogPoly(
  const Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly>& sg_mv_,
  bool is_multi_vec_transposed) : 
  ProductContainer<Epetra_Operator>(sg_mv_->map()),
  VectorOrthogPoly<Epetra_Operator>(sg_mv_->basis(), sg_mv_->map()),
  ProductEpetraOperator(sg_mv_->map(), sg_mv_->productComm()),
  EpetraOperatorOrthogPoly(sg_mv_->basis(), sg_mv_->map(), 
			   sg_mv_->productComm()),
  ProductEpetraMultiVectorOperator(sg_mv_, is_multi_vec_transposed),
  sg_mv(sg_mv_)
{
  // Now need to call setup() as it is called in 
  // ProductEpetraMultiVectorOperator constructor
}
    
Stokhos::EpetraMultiVectorOperatorOrthogPoly::
EpetraMultiVectorOperatorOrthogPoly(
  const Stokhos::EpetraMultiVectorOperatorOrthogPoly& v) :
  ProductContainer<Epetra_Operator>(v),
  VectorOrthogPoly<Epetra_Operator>(v),
  ProductEpetraOperator(v),
  EpetraOperatorOrthogPoly(v),
  ProductEpetraMultiVectorOperator(v)
{
  sg_mv = v.sg_mv;
}

Stokhos::EpetraMultiVectorOperatorOrthogPoly::
~EpetraMultiVectorOperatorOrthogPoly() {}

Stokhos::EpetraMultiVectorOperatorOrthogPoly& 
Stokhos::EpetraMultiVectorOperatorOrthogPoly::
operator=(const Stokhos::EpetraMultiVectorOperatorOrthogPoly& v) {
  ProductEpetraMultiVectorOperator::operator=(v);
  this->basis_ = v.basis_;
  sg_mv = v.sg_mv;
  return *this;
}

Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> 
Stokhos::EpetraMultiVectorOperatorOrthogPoly::
multiVectorOrthogPoly() const
{
  return sg_mv;
}
