// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
