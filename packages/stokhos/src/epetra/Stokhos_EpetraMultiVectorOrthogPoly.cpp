// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_EpetraMultiVectorOrthogPoly.hpp"

Stokhos::EpetraMultiVectorOrthogPoly::
EpetraMultiVectorOrthogPoly() :
  ProductContainer<Epetra_MultiVector>(),
  VectorOrthogPoly<Epetra_MultiVector>(),
  ProductEpetraMultiVector() 
{
}

Stokhos::EpetraMultiVectorOrthogPoly::
EpetraMultiVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  const Teuchos::RCP<const Epetra_BlockMap>& block_map) :
  ProductContainer<Epetra_MultiVector>(block_map),
  VectorOrthogPoly<Epetra_MultiVector>(basis, block_map),
  ProductEpetraMultiVector(block_map) 
{
}

Stokhos::EpetraMultiVectorOrthogPoly::
EpetraMultiVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm,
  int num_vectors) : 
  ProductContainer<Epetra_MultiVector>(block_map),
  VectorOrthogPoly<Epetra_MultiVector>(basis, block_map),
  ProductEpetraMultiVector(block_map, coeff_map, product_comm, num_vectors) 
{
}

Stokhos::EpetraMultiVectorOrthogPoly::
EpetraMultiVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
  const Teuchos::RCP<const Epetra_BlockMap>& product_map,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm,
  int num_vectors) :
  ProductContainer<Epetra_MultiVector>(block_map),
  VectorOrthogPoly<Epetra_MultiVector>(basis, block_map),
  ProductEpetraMultiVector(block_map, coeff_map, product_map, product_comm, 
			   num_vectors)
{
}

Stokhos::EpetraMultiVectorOrthogPoly::
EpetraMultiVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
  const Teuchos::RCP<const Epetra_BlockMap>& product_map,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm,
  Epetra_DataAccess CV,
  const Epetra_MultiVector& block_vector) :
  ProductContainer<Epetra_MultiVector>(block_map),
  VectorOrthogPoly<Epetra_MultiVector>(basis, block_map),
  ProductEpetraMultiVector(block_map, coeff_map, product_map, product_comm, CV,
			   block_vector)
{
}
    
Stokhos::EpetraMultiVectorOrthogPoly::
EpetraMultiVectorOrthogPoly(const Stokhos::EpetraMultiVectorOrthogPoly& v) :
  ProductContainer<Epetra_MultiVector>(v),
  VectorOrthogPoly<Epetra_MultiVector>(v),
  ProductEpetraMultiVector(v)
{
}

Stokhos::EpetraMultiVectorOrthogPoly::
~EpetraMultiVectorOrthogPoly() {}

Stokhos::EpetraMultiVectorOrthogPoly& 
Stokhos::EpetraMultiVectorOrthogPoly::
operator=(const Stokhos::EpetraMultiVectorOrthogPoly& v) {
  ProductEpetraMultiVector::operator=(v);
  this->basis_ = v.basis_;
  return *this;
}
      
void 
Stokhos::EpetraMultiVectorOrthogPoly::
reset(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& new_basis,
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm,
  int num_vectors) 
{
  ProductEpetraMultiVector::reset(block_map, coeff_map, product_comm, 
				  num_vectors);
  this->basis_ = new_basis;
}

void 
Stokhos::EpetraMultiVectorOrthogPoly::
reset(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& new_basis,
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
  const Teuchos::RCP<const Epetra_BlockMap>& product_map,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm,
  int num_vectors) 
{
  ProductEpetraMultiVector::reset(block_map, coeff_map, product_map, 
				  product_comm, num_vectors);
  this->basis_ = new_basis;
}

void
Stokhos::EpetraMultiVectorOrthogPoly::
computeMean(Epetra_MultiVector& v) const
{
  int lid = this->map_->LID(0);
  v.Scale(1.0, *(this->coeff_[lid]));
}

void
Stokhos::EpetraMultiVectorOrthogPoly::
computeVariance(Epetra_MultiVector& v) const
{
  const Teuchos::Array<double>& nrm2 = this->basis_->norm_squared();
  v.PutScalar(0.0);
  int i_gid;
  for (int i=1; i<this->size(); i++) {
    i_gid = this->map_->GID(i);
    v.Multiply(nrm2[i_gid], *(this->coeff_[i]), *(this->coeff_[i]), 1.0);
  }
}

void
Stokhos::EpetraMultiVectorOrthogPoly::
computeStandardDeviation(Epetra_MultiVector& v) const
{
  computeVariance(v);

  for (int j=0; j<v.NumVectors(); j++)
    for (int i=0; i<v.MyLength(); i++)
      v[j][i] = std::sqrt(v[j][i]);
}
