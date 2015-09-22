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

#include "Stokhos_MPBlockDiagonalPreconditioner.hpp"
#include "EpetraExt_BlockMultiVector.h"

Stokhos::MPBlockDiagonalPreconditioner::
MPBlockDiagonalPreconditioner(
  const Teuchos::RCP<const EpetraExt::MultiComm>& mp_comm_,
  int num_mp_blocks_,
  const Teuchos::RCP<const Epetra_Map>& base_map_,
  const Teuchos::RCP<const Epetra_Map>& mp_map_,
  const Teuchos::RCP<Stokhos::AbstractPreconditionerFactory>& prec_factory_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_) :
  label("Stokhos MP Block Diagonal Preconditioner"),
  mp_comm(mp_comm_),
  num_mp_blocks(num_mp_blocks_),
  base_map(base_map_),
  mp_map(mp_map_),
  prec_factory(prec_factory_),
  block_precs(num_mp_blocks),
  useTranspose(false)
{
}

Stokhos::MPBlockDiagonalPreconditioner::
~MPBlockDiagonalPreconditioner()
{
}

void
Stokhos::MPBlockDiagonalPreconditioner::
setupPreconditioner(const Teuchos::RCP<Stokhos::BlockDiagonalOperator>& mp_op,
                    const Epetra_Vector& x)
{
   TEUCHOS_TEST_FOR_EXCEPTION(prec_factory == Teuchos::null, std::logic_error,
                      "Error!  setupPreconditioner() cannot be called when " <<
                      "prec_factory is null!" << std::endl);

   Teuchos::RCP<Stokhos::ProductContainer<Epetra_Operator> > mp_ops =
     mp_op->getMPOps();
   for (int i=0; i<num_mp_blocks; i++)
     block_precs[i] = prec_factory->compute(mp_ops->getCoeffPtr(i));
   if (num_mp_blocks > 0) {
     label = std::string("Stokhos MP Block Diagonal Preconditioner:\n") +
       std::string("            ***** ") +
       std::string(block_precs[0]->Label());
   }
}

int
Stokhos::MPBlockDiagonalPreconditioner::
SetUseTranspose(bool UseTheTranspose)
{
  useTranspose = UseTheTranspose;
  for (int i=0; i<num_mp_blocks; i++)
    block_precs[i]->SetUseTranspose(useTranspose);

  return 0;
}

int
Stokhos::MPBlockDiagonalPreconditioner::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  EpetraExt::BlockMultiVector mp_input(View, *base_map, Input);
  EpetraExt::BlockMultiVector mp_result(View, *base_map, Result);
  for (int i=0; i<num_mp_blocks; i++) {
    block_precs[i]->Apply(*(mp_input.GetBlock(i)), *(mp_result.GetBlock(i)));
  }

  return 0;
}

int
Stokhos::MPBlockDiagonalPreconditioner::
ApplyInverse(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  EpetraExt::BlockMultiVector mp_input(View, *base_map, Input);
  EpetraExt::BlockMultiVector mp_result(View, *base_map, Result);
  for (int i=0; i<num_mp_blocks; i++) {
    block_precs[i]->ApplyInverse(*(mp_input.GetBlock(i)),
                                 *(mp_result.GetBlock(i)));
  }

  return 0;
}

double
Stokhos::MPBlockDiagonalPreconditioner::
NormInf() const
{
  double product_nrm = 0.0;
  for (int i=0; i<num_mp_blocks; i++) {
    double nrm = block_precs[i]->NormInf();
    if (nrm > product_nrm)
      product_nrm = nrm;
  }

  return product_nrm;
}


const char*
Stokhos::MPBlockDiagonalPreconditioner::
Label () const
{
  return const_cast<char*>(label.c_str());
}

bool
Stokhos::MPBlockDiagonalPreconditioner::
UseTranspose() const
{
  return useTranspose;
}

bool
Stokhos::MPBlockDiagonalPreconditioner::
HasNormInf() const
{
  if (num_mp_blocks == 0)
    return false;
  return block_precs[0]->HasNormInf();
}

const Epetra_Comm &
Stokhos::MPBlockDiagonalPreconditioner::
Comm() const
{
  return *mp_comm;
}
const Epetra_Map&
Stokhos::MPBlockDiagonalPreconditioner::
OperatorDomainMap() const
{
  return *mp_map;
}

const Epetra_Map&
Stokhos::MPBlockDiagonalPreconditioner::
OperatorRangeMap() const
{
  return *mp_map;
}
