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
  const Teuchos::RCP<const Epetra_Map>& _domain_base_map,
  const Teuchos::RCP<const Epetra_Map>& _range_base_map,
  const Teuchos::RCP<const EpetraExt::MultiComm>& _product_comm) :
  ProductContainer<Epetra_Operator>(block_map),
  VectorOrthogPoly<Epetra_Operator>(basis, block_map),
  ProductEpetraOperator(block_map, _domain_base_map, _range_base_map,
			_product_comm)
{
}

Stokhos::EpetraOperatorOrthogPoly::
EpetraOperatorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_Map>& _domain_base_map,
  const Teuchos::RCP<const Epetra_Map>& _range_base_map,
  const Teuchos::RCP<const Epetra_Map>& _product_range_map,
  const Teuchos::RCP<const EpetraExt::MultiComm>& _product_comm) :
  ProductContainer<Epetra_Operator>(block_map),
  VectorOrthogPoly<Epetra_Operator>(basis, block_map),
  ProductEpetraOperator(block_map, _domain_base_map, _range_base_map,
			_product_range_map, _product_comm)
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
  const Teuchos::RCP<const EpetraExt::MultiComm>& _product_comm_) :
  ProductContainer<Epetra_Operator>(block_map),
  VectorOrthogPoly<Epetra_Operator>(basis, block_map),
  ProductEpetraOperator(block_map, _product_comm_)
{
}

void
Stokhos::EpetraOperatorOrthogPoly::
setup(const Teuchos::RCP<const Epetra_Map>& _domain_base_map,
      const Teuchos::RCP<const Epetra_Map>& _range_base_map)
{
  ProductEpetraOperator::setup(_domain_base_map, _range_base_map);
}
