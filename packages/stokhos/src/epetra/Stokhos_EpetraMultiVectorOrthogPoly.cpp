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
