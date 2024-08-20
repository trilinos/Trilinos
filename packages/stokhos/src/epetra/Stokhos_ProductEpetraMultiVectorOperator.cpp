// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_ProductEpetraMultiVectorOperator.hpp"
#include "Stokhos_EpetraMultiVectorOperator.hpp"
#include "Epetra_LocalMap.h"

Stokhos::ProductEpetraMultiVectorOperator::
ProductEpetraMultiVectorOperator(
  const Teuchos::RCP<Stokhos::ProductEpetraMultiVector>& product_mv_,
  bool is_multi_vec_transposed) : 
  ProductContainer<Epetra_Operator>(product_mv_->map()),
  ProductEpetraOperator(product_mv_->map(),
			product_mv_->productComm()),
  product_mv(product_mv_)
{
  // Create domain, range maps
  Teuchos::RCP<const Epetra_Map> domain_map, range_map;
  int nv = product_mv->numVectors();
  const Epetra_Comm& coeff_comm = product_mv->productComm()->SubDomainComm();
  Teuchos::RCP<Epetra_LocalMap> local_map = 
    Teuchos::rcp(new Epetra_LocalMap(nv, 0, coeff_comm));
  if (is_multi_vec_transposed) {
    domain_map = 
      Teuchos::rcp_dynamic_cast<const Epetra_Map>(product_mv->coefficientMap());
    range_map = local_map;
  }
  else {
    domain_map = local_map;
    range_map = 
      Teuchos::rcp_dynamic_cast<const Epetra_Map>(product_mv->coefficientMap());
  }
  ProductEpetraOperator::setup(domain_map, range_map);

  // Set coefficients as operators
  for (int i=0; i<this->size(); i++) {
    Teuchos::RCP<Stokhos::EpetraMultiVectorOperator> op = 
      Teuchos::rcp(new Stokhos::EpetraMultiVectorOperator(
		     product_mv->getCoeffPtr(i), is_multi_vec_transposed));
    this->setCoeffPtr(i, op);
  }
}
    
Stokhos::ProductEpetraMultiVectorOperator::
ProductEpetraMultiVectorOperator(
  const Stokhos::ProductEpetraMultiVectorOperator& v) :
  ProductContainer<Epetra_Operator>(v),
  ProductEpetraOperator(v),
  product_mv(v.product_mv)
{
}

Stokhos::ProductEpetraMultiVectorOperator::
~ProductEpetraMultiVectorOperator() {}

Stokhos::ProductEpetraMultiVectorOperator& 
Stokhos::ProductEpetraMultiVectorOperator::
operator=(const Stokhos::ProductEpetraMultiVectorOperator& v) {
  ProductEpetraOperator::operator=(v);
  product_mv = v.product_mv;
  return *this;
}

Teuchos::RCP<Stokhos::ProductEpetraMultiVector> 
Stokhos::ProductEpetraMultiVectorOperator::
productMultiVector() const {
  return product_mv;
}
