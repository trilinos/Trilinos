// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_ProductEpetraMultiVector.hpp"
#include "EpetraExt_BlockUtility.h"
#include "Epetra_Map.h"

Stokhos::ProductEpetraMultiVector::
ProductEpetraMultiVector() :
  ProductContainer<Epetra_MultiVector>()
{
}

Stokhos::ProductEpetraMultiVector::
ProductEpetraMultiVector(const Teuchos::RCP<const Epetra_BlockMap>& block_map) :
  ProductContainer<Epetra_MultiVector>(block_map)
{
}

Stokhos::ProductEpetraMultiVector::
ProductEpetraMultiVector(
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map_,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_,
  int num_vectors) : 
  ProductContainer<Epetra_MultiVector>(block_map),
  coeff_map(coeff_map_),
  product_comm(product_comm_),
  product_map(Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*coeff_map,*block_map, *product_comm))),
  bv(Teuchos::rcp(new EpetraExt::BlockMultiVector(*coeff_map, *product_map,
						  num_vectors)))
{
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

Stokhos::ProductEpetraMultiVector::
ProductEpetraMultiVector(
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map_,
  const Teuchos::RCP<const Epetra_BlockMap>& product_map_,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_,
  int num_vectors) :
  ProductContainer<Epetra_MultiVector>(block_map),
  coeff_map(coeff_map_),
  product_comm(product_comm_),
  product_map(product_map_),
  bv(Teuchos::rcp(new EpetraExt::BlockMultiVector(*coeff_map, *product_map,
						  num_vectors)))
{
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

Stokhos::ProductEpetraMultiVector::
ProductEpetraMultiVector(
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map_,
  const Teuchos::RCP<const Epetra_BlockMap>& product_map_,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_,
  Epetra_DataAccess CV,
  const Epetra_MultiVector& block_vector) :
  ProductContainer<Epetra_MultiVector>(block_map),
  coeff_map(coeff_map_),
  product_comm(product_comm_),
  product_map(product_map_),
  bv(Teuchos::rcp(new EpetraExt::BlockMultiVector(CV, *coeff_map, 
						  block_vector)))
{
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}
    
Stokhos::ProductEpetraMultiVector::
ProductEpetraMultiVector(const Stokhos::ProductEpetraMultiVector& v) :
  ProductContainer<Epetra_MultiVector>(v),
  coeff_map(v.coeff_map),
  product_comm(v.product_comm),
  product_map(v.product_map),
  bv(v.bv)
{
}

Stokhos::ProductEpetraMultiVector::
~ProductEpetraMultiVector() {}

Stokhos::ProductEpetraMultiVector& 
Stokhos::ProductEpetraMultiVector::
operator=(const Stokhos::ProductEpetraMultiVector& v) {
  ProductContainer<Epetra_MultiVector>::operator=(v);
  coeff_map = v.coeff_map;
  product_comm = v.product_comm;
  product_map = v.product_map;
  bv = v.bv;  // Note this is a shallow copy, which is consistent with above
  return *this;
}

Stokhos::ProductEpetraMultiVector& 
Stokhos::ProductEpetraMultiVector::
operator=(const Epetra_MultiVector& v) {
  if (this->size() > 0) {
    if (bv != Teuchos::null)
      bv->Update(1.0, v, 0.0);
    else {
      EpetraExt::BlockMultiVector block_v(View, *coeff_map, v);
      for (int i=0; i<this->size(); i++)
	*(coeff_[i]) = *(block_v.GetBlock(i));
    }
  }
  return *this;
}

void 
Stokhos::ProductEpetraMultiVector::
assignToBlockMultiVector(Epetra_MultiVector& v) const 
{
  if (this->size() > 0) {
    if (bv != Teuchos::null)
      v.Update(1.0, *bv, 0.0);
    else {
      EpetraExt::BlockMultiVector block_v(View, *coeff_map, v);
      for (int i=0; i<this->size(); i++)
	*(block_v.GetBlock(i)) = *(coeff_[i]);
    }
  }
}

void 
Stokhos::ProductEpetraMultiVector::
assignFromBlockMultiVector(const Epetra_MultiVector& v) 
{
  if (this->size() > 0) {
    if (bv != Teuchos::null)
      bv->Update(1.0, v, 0.0);
    else {
      EpetraExt::BlockMultiVector block_v(View, *coeff_map, v);
      for (int i=0; i<this->size(); i++)
	*(coeff_[i]) = *(block_v.GetBlock(i));
    }
  }
}

Teuchos::RCP<const Epetra_BlockMap> 
Stokhos::ProductEpetraMultiVector::
coefficientMap() const {
  return coeff_map;
}

Teuchos::RCP<const Epetra_BlockMap> 
Stokhos::ProductEpetraMultiVector::
productMap() const {
  return product_map;
}

Teuchos::RCP<const EpetraExt::MultiComm> 
Stokhos::ProductEpetraMultiVector::
productComm() const {
  return product_comm;
}

int
Stokhos::ProductEpetraMultiVector::
numVectors() const {
  if (bv != Teuchos::null)
    return bv->NumVectors();
  else if (this->size() > 0 && this->coeff_[0] != Teuchos::null)
    return this->coeff_[0]->NumVectors();
  return -1;
}
      
void 
Stokhos::ProductEpetraMultiVector::
reset(const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map_,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_,
      int num_vectors) 
{
  ProductContainer<Epetra_MultiVector>::reset(block_map);
  coeff_map = coeff_map_;
  product_comm = product_comm_;
  product_map = 
    Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*coeff_map,
							   *block_map, 
							   *product_comm));
  bv = Teuchos::rcp(new EpetraExt::BlockMultiVector(*coeff_map, *product_map, 
						    num_vectors));
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

void 
Stokhos::ProductEpetraMultiVector::
reset(const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map_,
      const Teuchos::RCP<const Epetra_BlockMap>& product_map_,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_,
      int num_vectors) 
{
  ProductContainer<Epetra_MultiVector>::reset(block_map);
  coeff_map = coeff_map_;
  product_comm = product_comm_;
  product_map = product_map_;
  bv = Teuchos::rcp(new EpetraExt::BlockMultiVector(*coeff_map, *product_map,
						    num_vectors));
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

void 
Stokhos::ProductEpetraMultiVector::
resetCoefficients(Epetra_DataAccess CV,
		  const Epetra_MultiVector& block_vector) 
{
  bv = 
    Teuchos::rcp(new EpetraExt::BlockMultiVector(CV, *coeff_map, block_vector));
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

Teuchos::RCP<EpetraExt::BlockMultiVector> 
Stokhos::ProductEpetraMultiVector::
getBlockMultiVector() 
{
  return bv;
}

Teuchos::RCP<const EpetraExt::BlockMultiVector> 
Stokhos::ProductEpetraMultiVector::
getBlockMultiVector() const 
{
  return bv;
}

void 
Stokhos::ProductEpetraMultiVector::
setBlockMultiVector(const Teuchos::RCP<EpetraExt::BlockMultiVector>& block_vec) 
{
  bv = block_vec;
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}
