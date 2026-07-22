// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_MPMeanBasedPreconditioner.hpp"
#include "EpetraExt_BlockMultiVector.h"

Stokhos::MPMeanBasedPreconditioner::
MPMeanBasedPreconditioner(
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
  mean_prec(),
  useTranspose(false)
{
}

Stokhos::MPMeanBasedPreconditioner::
~MPMeanBasedPreconditioner()
{
}

void
Stokhos::MPMeanBasedPreconditioner::
setupPreconditioner(const Teuchos::RCP<Stokhos::BlockDiagonalOperator>& mp_op,
                    const Epetra_Vector& x)
{
   TEUCHOS_TEST_FOR_EXCEPTION(prec_factory == Teuchos::null, std::logic_error,
                      "Error!  setupPreconditioner() cannot be called when " <<
                      "prec_factory is null!" << std::endl);

   Teuchos::RCP<Stokhos::ProductContainer<Epetra_Operator> > mp_ops =
     mp_op->getMPOps();
   mean_prec = prec_factory->compute(mp_ops->getCoeffPtr(0));
   if (num_mp_blocks > 0) {
     label = std::string("Stokhos MP Mean Preconditioner:\n") +
       std::string("            ***** ") +
       std::string(mean_prec->Label());
   }
}

int
Stokhos::MPMeanBasedPreconditioner::
SetUseTranspose(bool UseTheTranspose)
{
  useTranspose = UseTheTranspose;
  mean_prec->SetUseTranspose(useTranspose);

  return 0;
}

int
Stokhos::MPMeanBasedPreconditioner::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  EpetraExt::BlockMultiVector mp_input(View, *base_map, Input);
  EpetraExt::BlockMultiVector mp_result(View, *base_map, Result);
  for (int i=0; i<num_mp_blocks; i++) {
    mean_prec->Apply(*(mp_input.GetBlock(i)), *(mp_result.GetBlock(i)));
  }

  return 0;
}

int
Stokhos::MPMeanBasedPreconditioner::
ApplyInverse(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  EpetraExt::BlockMultiVector mp_input(View, *base_map, Input);
  EpetraExt::BlockMultiVector mp_result(View, *base_map, Result);
  for (int i=0; i<num_mp_blocks; i++) {
    mean_prec->ApplyInverse(*(mp_input.GetBlock(i)), *(mp_result.GetBlock(i)));
  }

  return 0;
}

double
Stokhos::MPMeanBasedPreconditioner::
NormInf() const
{
  return mean_prec->NormInf();
}


const char*
Stokhos::MPMeanBasedPreconditioner::
Label () const
{
  return const_cast<char*>(label.c_str());
}

bool
Stokhos::MPMeanBasedPreconditioner::
UseTranspose() const
{
  return useTranspose;
}

bool
Stokhos::MPMeanBasedPreconditioner::
HasNormInf() const
{
  return mean_prec->HasNormInf();
}

const Epetra_Comm &
Stokhos::MPMeanBasedPreconditioner::
Comm() const
{
  return *mp_comm;
}
const Epetra_Map&
Stokhos::MPMeanBasedPreconditioner::
OperatorDomainMap() const
{
  return *mp_map;
}

const Epetra_Map&
Stokhos::MPMeanBasedPreconditioner::
OperatorRangeMap() const
{
  return *mp_map;
}
