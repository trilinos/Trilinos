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

#include "Stokhos_ProductEpetraOperator.hpp"
#include "EpetraExt_BlockUtility.h"
#include "EpetraExt_BlockMultiVector.h"

Stokhos::ProductEpetraOperator::
ProductEpetraOperator(
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_Map>& domain_base_map_,
  const Teuchos::RCP<const Epetra_Map>& range_base_map_,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_) : 
  ProductContainer<Epetra_Operator>(block_map),
  domain_base_map(domain_base_map_),
  range_base_map(range_base_map_),
  product_range_map(Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*range_base_map, *block_map, *product_comm_))),
  product_comm(product_comm_),
  useTranspose(false)
{
}

Stokhos::ProductEpetraOperator::
ProductEpetraOperator(
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_Map>& domain_base_map_,
  const Teuchos::RCP<const Epetra_Map>& range_base_map_,
  const Teuchos::RCP<const Epetra_Map>& product_range_map_,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_) :
  ProductContainer<Epetra_Operator>(block_map),
  domain_base_map(domain_base_map_),
  range_base_map(range_base_map_),
  product_range_map(product_range_map_),
  product_comm(product_comm_),
  useTranspose(false)
{
}
    
Stokhos::ProductEpetraOperator::
ProductEpetraOperator(const Stokhos::ProductEpetraOperator& v) :
  ProductContainer<Epetra_Operator>(v),
  domain_base_map(v.domain_base_map),
  range_base_map(v.range_base_map),
  product_range_map(v.product_range_map),
  product_comm(v.product_comm),
  useTranspose(v.useTranspose)
{
}

Stokhos::ProductEpetraOperator::
~ProductEpetraOperator() {}

Stokhos::ProductEpetraOperator& 
Stokhos::ProductEpetraOperator::
operator=(const Stokhos::ProductEpetraOperator& v) {
  ProductContainer<Epetra_Operator>::operator=(v);
  domain_base_map = v.domain_base_map;
  range_base_map = v.range_base_map;
  product_range_map = v.product_range_map;
  product_comm = v.product_comm;
  useTranspose = v.useTranspose;
  return *this;
}

Teuchos::RCP<const EpetraExt::MultiComm> 
Stokhos::ProductEpetraOperator::
productComm() const {
  return product_comm;
}

int 
Stokhos::ProductEpetraOperator::
SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  for (int i=0; i<coeff_.size(); i++)
    coeff_[i]->SetUseTranspose(useTranspose);

  return 0;
}

int 
Stokhos::ProductEpetraOperator::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  if (useTranspose) {
    EpetraExt::BlockMultiVector sg_input(View, *range_base_map, Input);
    Epetra_MultiVector tmp(Result.Map(), Result.NumVectors());
    Result.PutScalar(0.0);
    for (int i=0; i<coeff_.size(); i++) {
      coeff_[i]->Apply(*(sg_input.GetBlock(i)), tmp);
      Result.Update(1.0, tmp, 1.0);
    }
  }
  else {
    EpetraExt::BlockMultiVector sg_result(View, *range_base_map, Result);
    for (int i=0; i<coeff_.size(); i++)
      coeff_[i]->Apply(Input, *(sg_result.GetBlock(i)));
  }

  return 0;
}

int 
Stokhos::ProductEpetraOperator::
ApplyInverse(const Epetra_MultiVector& Input, 
	     Epetra_MultiVector& Result) const
{
  throw "ProductEpetraOperator::ApplyInverse not defined!";
  return -1;
}

double 
Stokhos::ProductEpetraOperator::
NormInf() const
{
  return 1.0;
}


const char* 
Stokhos::ProductEpetraOperator::
Label () const
{
  return "Stokhos::ProductEpetraOperator";
}
  
bool 
Stokhos::ProductEpetraOperator::
UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::ProductEpetraOperator::
HasNormInf() const
{
  return false;
}

const Epetra_Comm & 
Stokhos::ProductEpetraOperator::
Comm() const
{
  return *product_comm;
}
const Epetra_Map& 
Stokhos::ProductEpetraOperator::
OperatorDomainMap() const
{
  if (useTranspose)
    return *product_range_map;
  return *domain_base_map;
}

const Epetra_Map& 
Stokhos::ProductEpetraOperator::
OperatorRangeMap() const
{
  if (useTranspose)
    return *domain_base_map;
  return *product_range_map;
}
