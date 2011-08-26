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
