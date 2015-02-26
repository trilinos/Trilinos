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

#include "Stokhos_BlockDiagonalOperator.hpp"
#include "EpetraExt_BlockMultiVector.h"
#include "Epetra_Map.h"

Stokhos::BlockDiagonalOperator::
BlockDiagonalOperator(
  const Teuchos::RCP<const EpetraExt::MultiComm>& mp_comm_,
  int num_mp_blocks_,
  const Teuchos::RCP<const Epetra_Map>& domain_base_map_,
  const Teuchos::RCP<const Epetra_Map>& range_base_map_,
  const Teuchos::RCP<const Epetra_Map>& domain_mp_map_,
  const Teuchos::RCP<const Epetra_Map>& range_mp_map_) :
  label("Stokhos Block Diagonal Operator"),
  mp_comm(mp_comm_),
  num_mp_blocks(num_mp_blocks_),
  domain_base_map(domain_base_map_),
  range_base_map(range_base_map_),
  domain_mp_map(domain_mp_map_),
  range_mp_map(range_mp_map_),
  block_ops(),
  useTranspose(false)
{
}

Stokhos::BlockDiagonalOperator::
~BlockDiagonalOperator()
{
}

void
Stokhos::BlockDiagonalOperator::
setupOperator(
   const Teuchos::RCP<Stokhos::ProductEpetraOperator >& ops)
{
  block_ops = ops;
}

Teuchos::RCP< Stokhos::ProductEpetraOperator >
Stokhos::BlockDiagonalOperator::
getMPOps()
{
  return block_ops;
}

Teuchos::RCP<const Stokhos::ProductEpetraOperator >
Stokhos::BlockDiagonalOperator::
getMPOps() const
{
  return block_ops;
}

int
Stokhos::BlockDiagonalOperator::
SetUseTranspose(bool UseTheTranspose)
{
  useTranspose = UseTheTranspose;
  for (int i=0; i<num_mp_blocks; i++)
    (*block_ops)[i].SetUseTranspose(useTranspose);

  return 0;
}

int
Stokhos::BlockDiagonalOperator::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  Teuchos::RCP<const Epetra_BlockMap> input_base_map, result_base_map;
  if (useTranspose) {
    input_base_map = range_base_map;
    result_base_map = domain_base_map;
  }
  else {
    input_base_map = domain_base_map;
    result_base_map = range_base_map;
  }
  EpetraExt::BlockMultiVector mp_input(View, *input_base_map, Input);
  EpetraExt::BlockMultiVector mp_result(View, *result_base_map, Result);
  for (int i=0; i<num_mp_blocks; i++) {
    (*block_ops)[i].Apply(*(mp_input.GetBlock(i)), *(mp_result.GetBlock(i)));
  }

  return 0;
}

int
Stokhos::BlockDiagonalOperator::ApplyInverse(const Epetra_MultiVector& Input,
                                             Epetra_MultiVector& Result) const
{
  Teuchos::RCP<const Epetra_BlockMap> input_base_map, result_base_map;
  if (useTranspose) {
    input_base_map = domain_base_map;
    result_base_map = range_base_map;
  }
  else {
    input_base_map = range_base_map;
    result_base_map = domain_base_map;
  }
  EpetraExt::BlockMultiVector mp_input(View, *input_base_map, Input);
  EpetraExt::BlockMultiVector mp_result(View, *range_base_map, Result);
  for (int i=0; i<num_mp_blocks; i++) {
    (*block_ops)[i].ApplyInverse(*(mp_input.GetBlock(i)),
                                 *(mp_result.GetBlock(i)));
  }

  return 0;
}

double
Stokhos::BlockDiagonalOperator::NormInf() const
{
  double product_nrm = 0.0;
  for (int i=0; i<num_mp_blocks; i++) {
    double nrm = (*block_ops)[i].NormInf();
    if (nrm > product_nrm)
      product_nrm = nrm;
  }

  return product_nrm;
}


const char*
Stokhos::BlockDiagonalOperator::Label () const
{
  return const_cast<char*>(label.c_str());
}

bool
Stokhos::BlockDiagonalOperator::UseTranspose() const
{
  return useTranspose;
}

bool
Stokhos::BlockDiagonalOperator::HasNormInf() const
{
  if (num_mp_blocks == 0)
    return false;
  return (*block_ops)[0].HasNormInf();
}

const Epetra_Comm &
Stokhos::BlockDiagonalOperator::Comm() const
{
  return *mp_comm;
}
const Epetra_Map&
Stokhos::BlockDiagonalOperator::OperatorDomainMap() const
{
  if (useTranspose)
    return *range_mp_map;
  return *domain_mp_map;
}

const Epetra_Map&
Stokhos::BlockDiagonalOperator::OperatorRangeMap() const
{
  if (useTranspose)
    return *domain_mp_map;
  return *range_mp_map;
}
