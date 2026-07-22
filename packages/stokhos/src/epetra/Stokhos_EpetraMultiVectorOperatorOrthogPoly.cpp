// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
