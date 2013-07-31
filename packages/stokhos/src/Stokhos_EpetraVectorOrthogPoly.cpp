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

#include "Stokhos_EpetraVectorOrthogPoly.hpp"

Stokhos::EpetraVectorOrthogPoly::
EpetraVectorOrthogPoly() :
  ProductContainer<Epetra_Vector>(),
  VectorOrthogPoly<Epetra_Vector>(),
  ProductEpetraVector()
{
}

Stokhos::EpetraVectorOrthogPoly::
EpetraVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  const Teuchos::RCP<const Epetra_BlockMap>& block_map) :
  ProductContainer<Epetra_Vector>(block_map),
  VectorOrthogPoly<Epetra_Vector>(basis, block_map),
  ProductEpetraVector(block_map) 
{
}

Stokhos::EpetraVectorOrthogPoly::
EpetraVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm) : 
  ProductContainer<Epetra_Vector>(block_map),
  VectorOrthogPoly<Epetra_Vector>(basis, block_map),
  ProductEpetraVector(block_map, coeff_map, product_comm) 
{
}

Stokhos::EpetraVectorOrthogPoly::
EpetraVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
  const Teuchos::RCP<const Epetra_BlockMap>& product_map,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm) :
  ProductContainer<Epetra_Vector>(block_map),
  VectorOrthogPoly<Epetra_Vector>(basis, block_map),
  ProductEpetraVector(block_map, coeff_map, product_map, product_comm) 
{
}

Stokhos::EpetraVectorOrthogPoly::
EpetraVectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
  const Teuchos::RCP<const Epetra_BlockMap>& product_map,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm,
  Epetra_DataAccess CV,
  const Epetra_Vector& block_vector) :
  ProductContainer<Epetra_Vector>(block_map),
  VectorOrthogPoly<Epetra_Vector>(basis, block_map),
  ProductEpetraVector(block_map, coeff_map, product_map, product_comm, CV,
		      block_vector) 
{
}
    
Stokhos::EpetraVectorOrthogPoly::
EpetraVectorOrthogPoly(const Stokhos::EpetraVectorOrthogPoly& v) :
  ProductContainer<Epetra_Vector>(v),
  VectorOrthogPoly<Epetra_Vector>(v),
  ProductEpetraVector(v)
{
}

Stokhos::EpetraVectorOrthogPoly::
~EpetraVectorOrthogPoly() {}

Stokhos::EpetraVectorOrthogPoly& 
Stokhos::EpetraVectorOrthogPoly::
operator=(const Stokhos::EpetraVectorOrthogPoly& v) {
  ProductEpetraVector::operator=(v);
  this->basis_ = v.basis_;
  return *this;
}
      
void 
Stokhos::EpetraVectorOrthogPoly::
reset(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& new_basis,
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm) 
{
  ProductEpetraVector::reset(block_map, coeff_map, product_comm);
  this->basis_ = new_basis;
}

void 
Stokhos::EpetraVectorOrthogPoly::
reset(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& new_basis,
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
  const Teuchos::RCP<const Epetra_BlockMap>& product_map,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm) 
{
  ProductEpetraVector::reset(block_map, coeff_map, product_map, product_comm);
  this->basis_ = new_basis;
}

void
Stokhos::EpetraVectorOrthogPoly::
computeMean(Epetra_Vector& v) const
{
  if (this->map_->Comm().NumProc() == 1 || !this->map_->DistributedGlobal()) {
    v.Scale(1.0, *(this->coeff_[0]));
    return;
  }

  int root = -1;
  int gid = 0;
  int lid = -1;
  this->map_->RemoteIDList(1, &gid, &root, &lid);
  if (this->map_->Comm().MyPID() == root) {
    v.Scale(1.0, *(this->coeff_[lid]));
  }
  this->map_->Comm().Broadcast(v.Values(), v.MyLength(), root);
}

void
Stokhos::EpetraVectorOrthogPoly::
computeVariance(Epetra_Vector& v) const
{
  Epetra_Vector *v_local = NULL;
  bool is_parallel = (this->map_->Comm().NumProc() > 1) && 
    this->map_->DistributedGlobal();
  if (is_parallel)
    v_local = new Epetra_Vector(v.Map());
  else
    v_local = &v;

  // Compute partial variance with terms on this processor
  const Teuchos::Array<double>& nrm2 = this->basis_->norm_squared();
  v_local->PutScalar(0.0);
  int i_gid;
  for (int i=0; i<this->size(); i++) {
    i_gid = this->map_->GID(i);
    if (i_gid != 0) 
      v_local->Multiply(nrm2[i_gid], *(this->coeff_[i]), *(this->coeff_[i]), 
			1.0);
  }
  
  // Compute total variance by summing across all processors
  if (is_parallel)
    this->map_->Comm().SumAll(v_local->Values(), v.Values(), v.MyLength());
}

void
Stokhos::EpetraVectorOrthogPoly::
computeStandardDeviation(Epetra_Vector& v) const
{
  // Compute variance
  computeVariance(v);

  // Compute standard deviation
  for (int i=0; i<v.MyLength(); i++)
    v[i] = std::sqrt(v[i]);
}
