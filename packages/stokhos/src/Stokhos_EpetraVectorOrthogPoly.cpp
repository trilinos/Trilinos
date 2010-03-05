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
#include "Stokhos_EpetraVectorOrthogPoly.hpp"
Stokhos::EpetraVectorOrthogPoly::
EpetraVectorOrthogPoly() :
  VectorOrthogPoly<Epetra_Vector>(),
  bv() 
{
}

Stokhos::EpetraVectorOrthogPoly::
EpetraVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis) :
  VectorOrthogPoly<Epetra_Vector>(basis),
  bv() 
{
}

Stokhos::EpetraVectorOrthogPoly::
EpetraVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  int sz) : 
  VectorOrthogPoly<Epetra_Vector>(basis, sz),
  bv() 
{
}

Stokhos::EpetraVectorOrthogPoly::
EpetraVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  const Epetra_BlockMap& coeff_map) : 
  VectorOrthogPoly<Epetra_Vector>(basis, EpetraVectorCloner(coeff_map)),
  bv() 
{
}

Stokhos::EpetraVectorOrthogPoly::
EpetraVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  const Epetra_BlockMap& coeff_map,
  int sz) : 
  VectorOrthogPoly<Epetra_Vector>(basis, EpetraVectorCloner(coeff_map), sz),
  bv() 
{
}

Stokhos::EpetraVectorOrthogPoly::
EpetraVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  Epetra_DataAccess CV,
  const Epetra_BlockMap& coeff_map,
  const Epetra_BlockMap& block_map) :
  VectorOrthogPoly<Epetra_Vector>(basis),
  bv(Teuchos::rcp(new EpetraExt::BlockVector(coeff_map, block_map)))
{
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

Stokhos::EpetraVectorOrthogPoly::
EpetraVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  Epetra_DataAccess CV,
  const Epetra_BlockMap& coeff_map,
  const Epetra_Vector& block_vector) :
  VectorOrthogPoly<Epetra_Vector>(basis),
  bv(Teuchos::rcp(new EpetraExt::BlockVector(CV, coeff_map, block_vector)))
{
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}
    
Stokhos::EpetraVectorOrthogPoly::
EpetraVectorOrthogPoly(const Stokhos::EpetraVectorOrthogPoly& v) :
  VectorOrthogPoly<Epetra_Vector>(v),
  bv(v.bv)
{
}

Stokhos::EpetraVectorOrthogPoly::
~EpetraVectorOrthogPoly() {}

Stokhos::EpetraVectorOrthogPoly& 
Stokhos::EpetraVectorOrthogPoly::
operator=(const Stokhos::EpetraVectorOrthogPoly& v) {
  VectorOrthogPoly<Epetra_Vector>::operator=(v);
  bv = v.bv;
  return *this;
}

Stokhos::EpetraVectorOrthogPoly& 
Stokhos::EpetraVectorOrthogPoly::
operator=(const Epetra_Vector& v) {
  if (this->size() > 0) {
    if (bv != Teuchos::null)
      bv->Update(1.0, v, 0.0);
    else {
      EpetraExt::BlockVector block_v(View, this->coeff_[0]->Map(), v);
      for (int i=0; i<this->size(); i++)
	*(coeff_[i]) = *(block_v.GetBlock(i));
    }
  }
  return *this;
}

void 
Stokhos::EpetraVectorOrthogPoly::
assignToBlockVector(Epetra_Vector& v) const 
{
  if (this->size() > 0) {
    if (bv != Teuchos::null)
      v.Update(1.0, *bv, 0.0);
    else {
      EpetraExt::BlockVector block_v(View, this->coeff_[0]->Map(), v);
      for (int i=0; i<this->size(); i++)
	*(block_v.GetBlock(i)) = *(coeff_[i]);
    }
  }
}

void 
Stokhos::EpetraVectorOrthogPoly::
assignFromBlockVector(const Epetra_Vector& v) 
{
  if (this->size() > 0) {
    if (bv != Teuchos::null)
      bv->Update(1.0, v, 0.0);
    else {
      EpetraExt::BlockVector block_v(View, this->coeff_[0]->Map(), v);
      for (int i=0; i<this->size(); i++)
	*(coeff_[i]) = *(block_v.GetBlock(i));
    }
  }
}
      
void 
Stokhos::EpetraVectorOrthogPoly::
reset(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& new_basis,
  const Epetra_BlockMap& coeff_map) 
{
  VectorOrthogPoly<Epetra_Vector>::reset(new_basis, 
					 EpetraVectorCloner(coeff_map));
  bv = Teuchos::null;
}

void 
Stokhos::EpetraVectorOrthogPoly::
resetCoefficients(Epetra_DataAccess CV,
		  const Epetra_BlockMap& coeff_map,
		  const Epetra_Vector& block_vector) 
{
  bv = 
    Teuchos::rcp(new EpetraExt::BlockVector(CV, coeff_map, block_vector));
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

Teuchos::RCP<EpetraExt::BlockVector> 
Stokhos::EpetraVectorOrthogPoly::
getBlockVector() 
{
  return bv;
}

Teuchos::RCP<const EpetraExt::BlockVector> 
Stokhos::EpetraVectorOrthogPoly::
getBlockVector() const 
{
  return bv;
}

void 
Stokhos::EpetraVectorOrthogPoly::
setBlockVector(const Teuchos::RCP<EpetraExt::BlockVector>& block_vec) 
{
  bv = block_vec;
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

