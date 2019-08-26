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

Stokhos::EpetraOperatorOrthogPoly::
EpetraOperatorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_) : 
  ProductContainer<Epetra_Operator>(block_map),
  VectorOrthogPoly<Epetra_Operator>(basis, block_map),
  ProductEpetraOperator(block_map, product_comm_) 
{
}

void
Stokhos::EpetraOperatorOrthogPoly::
setup(const Teuchos::RCP<const Epetra_Map>& domain_base_map,
      const Teuchos::RCP<const Epetra_Map>& range_base_map)
{
  ProductEpetraOperator::setup(domain_base_map, range_base_map);
}
