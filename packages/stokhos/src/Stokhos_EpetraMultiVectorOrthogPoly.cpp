// $Id$ 
// $Source$ 
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
#include "Stokhos_EpetraMultiVectorOrthogPoly.hpp"

Stokhos::EpetraMultiVectorOrthogPoly::
EpetraMultiVectorOrthogPoly() :
  VectorOrthogPoly<Epetra_MultiVector>(),
  bv() 
{
}

Stokhos::EpetraMultiVectorOrthogPoly::
EpetraMultiVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis) :
  VectorOrthogPoly<Epetra_MultiVector>(basis),
  bv() 
{
}

Stokhos::EpetraMultiVectorOrthogPoly::
EpetraMultiVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  int sz) : 
  VectorOrthogPoly<Epetra_MultiVector>(basis, sz),
  bv() 
{
}

Stokhos::EpetraMultiVectorOrthogPoly::
EpetraMultiVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  const Epetra_BlockMap& coeff_map,
  int num_vectors) : 
  VectorOrthogPoly<Epetra_MultiVector>(basis, 
				       EpetraMultiVectorCloner(coeff_map, 
							       num_vectors)),
  bv() 
{
}

Stokhos::EpetraMultiVectorOrthogPoly::
EpetraMultiVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  const Epetra_BlockMap& coeff_map,
  int num_vectors,
  int sz) : 
  VectorOrthogPoly<Epetra_MultiVector>(basis, 
				       EpetraMultiVectorCloner(coeff_map,
							       num_vectors), 
				       sz),
  bv() 
{
}

Stokhos::EpetraMultiVectorOrthogPoly::
EpetraMultiVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  Epetra_DataAccess CV,
  const Epetra_BlockMap& coeff_map,
  const Epetra_BlockMap& block_map,
  int num_vectors) :
  VectorOrthogPoly<Epetra_MultiVector>(basis),
  bv(Teuchos::rcp(new EpetraExt::BlockMultiVector(coeff_map, block_map, 
						  num_vectors)))
{
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

Stokhos::EpetraMultiVectorOrthogPoly::
EpetraMultiVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  Epetra_DataAccess CV,
  const Epetra_BlockMap& coeff_map,
  const Epetra_MultiVector& block_vector) :
  VectorOrthogPoly<Epetra_MultiVector>(basis),
  bv(Teuchos::rcp(new EpetraExt::BlockMultiVector(CV, coeff_map, block_vector)))
{
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}
    
Stokhos::EpetraMultiVectorOrthogPoly::
EpetraMultiVectorOrthogPoly(const Stokhos::EpetraMultiVectorOrthogPoly& v) :
  VectorOrthogPoly<Epetra_MultiVector>(v),
  bv(v.bv)
{
}

Stokhos::EpetraMultiVectorOrthogPoly::
~EpetraMultiVectorOrthogPoly() {}

Stokhos::EpetraMultiVectorOrthogPoly& 
Stokhos::EpetraMultiVectorOrthogPoly::
operator=(const Stokhos::EpetraMultiVectorOrthogPoly& v) {
  VectorOrthogPoly<Epetra_MultiVector>::operator=(v);
  bv = v.bv;
  return *this;
}

Stokhos::EpetraMultiVectorOrthogPoly& 
Stokhos::EpetraMultiVectorOrthogPoly::
operator=(const Epetra_MultiVector& v) {
  if (this->size() > 0) {
    if (bv != Teuchos::null)
      bv->Update(1.0, v, 0.0);
    else {
      EpetraExt::BlockMultiVector block_v(View, this->coeff_[0]->Map(), v);
      for (int i=0; i<this->size(); i++)
	*(coeff_[i]) = *(block_v.GetBlock(i));
    }
  }
  return *this;
}

void 
Stokhos::EpetraMultiVectorOrthogPoly::
assignToBlockMultiVector(Epetra_MultiVector& v) const 
{
  if (this->size() > 0) {
    if (bv != Teuchos::null)
      v.Update(1.0, *bv, 0.0);
    else {
      EpetraExt::BlockMultiVector block_v(View, this->coeff_[0]->Map(), v);
      for (int i=0; i<this->size(); i++)
	*(block_v.GetBlock(i)) = *(coeff_[i]);
    }
  }
}

void 
Stokhos::EpetraMultiVectorOrthogPoly::
assignFromBlockMultiVector(const Epetra_MultiVector& v) 
{
  if (this->size() > 0) {
    if (bv != Teuchos::null)
      bv->Update(1.0, v, 0.0);
    else {
      EpetraExt::BlockMultiVector block_v(View, this->coeff_[0]->Map(), v);
      for (int i=0; i<this->size(); i++)
	*(coeff_[i]) = *(block_v.GetBlock(i));
    }
  }
}
      
void 
Stokhos::EpetraMultiVectorOrthogPoly::
reset(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& new_basis,
  const Epetra_BlockMap& coeff_map,
  int num_vectors) 
{
  VectorOrthogPoly<Epetra_MultiVector>::reset(
    new_basis, EpetraMultiVectorCloner(coeff_map, num_vectors));
  bv = Teuchos::null;
}

void 
Stokhos::EpetraMultiVectorOrthogPoly::
resetCoefficients(Epetra_DataAccess CV,
		  const Epetra_BlockMap& coeff_map,
		  const Epetra_MultiVector& block_vector) 
{
  bv = 
    Teuchos::rcp(new EpetraExt::BlockMultiVector(CV, coeff_map, block_vector));
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

Teuchos::RCP<EpetraExt::BlockMultiVector> 
Stokhos::EpetraMultiVectorOrthogPoly::
getBlockMultiVector() 
{
  return bv;
}

Teuchos::RCP<const EpetraExt::BlockMultiVector> 
Stokhos::EpetraMultiVectorOrthogPoly::
getBlockMultiVector() const 
{
  return bv;
}

void 
Stokhos::EpetraMultiVectorOrthogPoly::
setBlockMultiVector(const Teuchos::RCP<EpetraExt::BlockMultiVector>& block_vec) 
{
  bv = block_vec;
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

