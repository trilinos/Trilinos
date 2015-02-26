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
SetUseTranspose(bool UseTheTranspose)
{
  useTranspose = UseTheTranspose;
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

Stokhos::ProductEpetraOperator::
ProductEpetraOperator(
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_) :
  ProductContainer<Epetra_Operator>(block_map),
  product_comm(product_comm_),
  useTranspose(false)
{
}

void
Stokhos::ProductEpetraOperator::
setup(const Teuchos::RCP<const Epetra_Map>& domain_base_map_,
      const Teuchos::RCP<const Epetra_Map>& range_base_map_)
{
  domain_base_map = domain_base_map_;
  range_base_map = range_base_map_;
  product_range_map =
    Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*range_base_map,
                                                           *(this->map_),
                                                           *product_comm));
}
