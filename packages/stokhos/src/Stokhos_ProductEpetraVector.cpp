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

#include "Stokhos_ProductEpetraVector.hpp"
#include "EpetraExt_BlockUtility.h"
#include "Epetra_Map.h"

Stokhos::ProductEpetraVector::
ProductEpetraVector() :
  ProductContainer<Epetra_Vector>()
{
}

Stokhos::ProductEpetraVector::
ProductEpetraVector(const Teuchos::RCP<const Epetra_BlockMap>& block_map) :
  ProductContainer<Epetra_Vector>(block_map)
{
}

Stokhos::ProductEpetraVector::
ProductEpetraVector(
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map_,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_) : 
  ProductContainer<Epetra_Vector>(block_map),
  coeff_map(coeff_map_),
  product_comm(product_comm_),
  product_map(Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*coeff_map,*block_map, *product_comm))),
  bv(Teuchos::rcp(new EpetraExt::BlockVector(*coeff_map, *product_map)))
{
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

Stokhos::ProductEpetraVector::
ProductEpetraVector(
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map_,
  const Teuchos::RCP<const Epetra_BlockMap>& product_map_,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_) :
  ProductContainer<Epetra_Vector>(block_map),
  coeff_map(coeff_map_),
  product_comm(product_comm_),
  product_map(product_map_),
  bv(Teuchos::rcp(new EpetraExt::BlockVector(*coeff_map, *product_map)))
{
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

Stokhos::ProductEpetraVector::
ProductEpetraVector(
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map_,
  const Teuchos::RCP<const Epetra_BlockMap>& product_map_,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_,
  Epetra_DataAccess CV,
  const Epetra_Vector& block_vector) :
  ProductContainer<Epetra_Vector>(block_map),
  coeff_map(coeff_map_),
  product_comm(product_comm_),
  product_map(product_map_),
  bv(Teuchos::rcp(new EpetraExt::BlockVector(CV, *coeff_map, block_vector)))
{
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}
    
Stokhos::ProductEpetraVector::
ProductEpetraVector(const Stokhos::ProductEpetraVector& v) :
  ProductContainer<Epetra_Vector>(v),
  coeff_map(v.coeff_map),
  product_comm(v.product_comm),
  product_map(v.product_map),
  bv(v.bv)
{
}

Stokhos::ProductEpetraVector::
~ProductEpetraVector() {}

Stokhos::ProductEpetraVector& 
Stokhos::ProductEpetraVector::
operator=(const Stokhos::ProductEpetraVector& v) {
  ProductContainer<Epetra_Vector>::operator=(v);
  coeff_map = v.coeff_map;
  product_comm = v.product_comm;
  product_map = v.product_map;
  bv = v.bv;  // Note this is a shallow copy, which is consistent with above
  return *this;
}

Stokhos::ProductEpetraVector& 
Stokhos::ProductEpetraVector::
operator=(const Epetra_Vector& v) {
  if (this->size() > 0) {
    if (bv != Teuchos::null)
      bv->Update(1.0, v, 0.0);
    else {
      EpetraExt::BlockVector block_v(View, *coeff_map, v);
      for (int i=0; i<this->size(); i++)
	*(coeff_[i]) = *(block_v.GetBlock(i));
    }
  }
  return *this;
}

void 
Stokhos::ProductEpetraVector::
assignToBlockVector(Epetra_Vector& v) const 
{
  if (this->size() > 0) {
    if (bv != Teuchos::null)
      v.Update(1.0, *bv, 0.0);
    else {
      EpetraExt::BlockVector block_v(View, *coeff_map, v);
      for (int i=0; i<this->size(); i++)
	*(block_v.GetBlock(i)) = *(coeff_[i]);
    }
  }
}

void 
Stokhos::ProductEpetraVector::
assignFromBlockVector(const Epetra_Vector& v) 
{
  if (this->size() > 0) {
    if (bv != Teuchos::null)
      bv->Update(1.0, v, 0.0);
    else {
      EpetraExt::BlockVector block_v(View, *coeff_map, v);
      for (int i=0; i<this->size(); i++)
	*(coeff_[i]) = *(block_v.GetBlock(i));
    }
  }
}

Teuchos::RCP<const Epetra_BlockMap> 
Stokhos::ProductEpetraVector::
coefficientMap() const {
  return coeff_map;
}

Teuchos::RCP<const Epetra_BlockMap> 
Stokhos::ProductEpetraVector::
productMap() const {
  return product_map;
}

Teuchos::RCP<const EpetraExt::MultiComm> 
Stokhos::ProductEpetraVector::
productComm() const {
  return product_comm;
}

void 
Stokhos::ProductEpetraVector::
reset(const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map_,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_) 
{
  ProductContainer<Epetra_Vector>::reset(block_map);
  coeff_map = coeff_map_;
  product_comm = product_comm_;
  product_map = 
    Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*coeff_map,
							   *block_map, 
							   *product_comm));
  bv = Teuchos::rcp(new EpetraExt::BlockVector(*coeff_map, *product_map));
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

void 
Stokhos::ProductEpetraVector::
reset(const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map_,
      const Teuchos::RCP<const Epetra_BlockMap>& product_map_,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_) 
{
  ProductContainer<Epetra_Vector>::reset(block_map);
  coeff_map = coeff_map_;
  product_comm = product_comm_;
  product_map = product_map_;
  bv = Teuchos::rcp(new EpetraExt::BlockVector(*coeff_map, *product_map));
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

void 
Stokhos::ProductEpetraVector::
resetCoefficients(Epetra_DataAccess CV,
		  const Epetra_Vector& block_vector) 
{
  bv = Teuchos::rcp(new EpetraExt::BlockVector(CV, *coeff_map, block_vector));
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

Teuchos::RCP<EpetraExt::BlockVector> 
Stokhos::ProductEpetraVector::
getBlockVector() 
{
  return bv;
}

Teuchos::RCP<const EpetraExt::BlockVector> 
Stokhos::ProductEpetraVector::
getBlockVector() const 
{
  return bv;
}

void 
Stokhos::ProductEpetraVector::
setBlockVector(const Teuchos::RCP<EpetraExt::BlockVector>& block_vec) 
{
  bv = block_vec;
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

void
Stokhos::ProductEpetraVector::
sumAll()
{
  Epetra_Vector v(*coeff_map);
  int sz = this->size();
  for (int i=0; i<sz; i++) {
    v.Scale(1.0, *(this->coeff_[i]));
    this->map_->Comm().SumAll(v.Values(), 
			      this->coeff_[i]->Values(), 
			      this->coeff_[i]->MyLength());
  }
}
