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
