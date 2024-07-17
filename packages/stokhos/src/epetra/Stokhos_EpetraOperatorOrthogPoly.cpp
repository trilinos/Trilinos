// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_EpetraOperatorOrthogPoly.hpp"

Stokhos::EpetraOperatorOrthogPoly::
EpetraOperatorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_Map>& domain_base_map,
  const Teuchos::RCP<const Epetra_Map>& range_base_map,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm) : 
  ProductContainer<Epetra_Operator>(block_map),
  VectorOrthogPoly<Epetra_Operator>(basis, block_map),
  ProductEpetraOperator(block_map, domain_base_map, range_base_map, 
			product_comm) 
{
}

Stokhos::EpetraOperatorOrthogPoly::
EpetraOperatorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_Map>& domain_base_map,
  const Teuchos::RCP<const Epetra_Map>& range_base_map,
  const Teuchos::RCP<const Epetra_Map>& product_range_map,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm) :
  ProductContainer<Epetra_Operator>(block_map),
  VectorOrthogPoly<Epetra_Operator>(basis, block_map),
  ProductEpetraOperator(block_map, domain_base_map, range_base_map, 
			product_range_map, product_comm) 
{
}
    
Stokhos::EpetraOperatorOrthogPoly::
EpetraOperatorOrthogPoly(const Stokhos::EpetraOperatorOrthogPoly& v) :
  ProductContainer<Epetra_Operator>(v),
  VectorOrthogPoly<Epetra_Operator>(v),
  ProductEpetraOperator(v)
{
}

Stokhos::EpetraOperatorOrthogPoly::
~EpetraOperatorOrthogPoly() {}

Stokhos::EpetraOperatorOrthogPoly& 
Stokhos::EpetraOperatorOrthogPoly::
operator=(const Stokhos::EpetraOperatorOrthogPoly& v) {
  ProductEpetraOperator::operator=(v);
  this->basis_ = v.basis_;
  return *this;
}

Stokhos::EpetraOperatorOrthogPoly::
EpetraOperatorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_) : 
  ProductContainer<Epetra_Operator>(block_map),
  VectorOrthogPoly<Epetra_Operator>(basis, block_map),
  ProductEpetraOperator(block_map, product_comm_) 
{
}

void
Stokhos::EpetraOperatorOrthogPoly::
setup(const Teuchos::RCP<const Epetra_Map>& domain_base_map,
      const Teuchos::RCP<const Epetra_Map>& range_base_map)
{
  ProductEpetraOperator::setup(domain_base_map, range_base_map);
}
