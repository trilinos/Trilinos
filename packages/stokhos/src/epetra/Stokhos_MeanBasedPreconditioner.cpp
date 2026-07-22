// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_MeanBasedPreconditioner.hpp"
#include "EpetraExt_BlockMultiVector.h"

Stokhos::MeanBasedPreconditioner::
MeanBasedPreconditioner(
  const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm_,
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
  const Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>& epetraCijk_,
  const Teuchos::RCP<const Epetra_Map>& base_map_,
  const Teuchos::RCP<const Epetra_Map>& sg_map_,
  const Teuchos::RCP<Stokhos::AbstractPreconditionerFactory>& prec_factory_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_) :
  label("Stokhos Mean-Based Preconditioner"),
  sg_comm(sg_comm_),
  sg_basis(sg_basis_),
  epetraCijk(epetraCijk_),
  base_map(base_map_),
  sg_map(sg_map_),
  num_blocks(0),
  prec_factory(prec_factory_),
  mean_prec(),
  useTranspose(false),
  use_block_apply(true)
{
  use_block_apply = params_->get("Use Block Apply", true);
}

Stokhos::MeanBasedPreconditioner::
~MeanBasedPreconditioner()
{
}

void
Stokhos::MeanBasedPreconditioner::
setupPreconditioner(const Teuchos::RCP<Stokhos::SGOperator>& sg_op,
                    const Epetra_Vector& x)
{
   TEUCHOS_TEST_FOR_EXCEPTION(prec_factory == Teuchos::null, std::logic_error,
                      "Error!  setupPreconditioner() cannot be called when " <<
                      "prec_factory is null!" << std::endl);

   Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly > sg_poly =
     sg_op->getSGPolynomial();
   mean_prec = prec_factory->compute(sg_poly->getCoeffPtr(0));
   label = std::string("Stokhos Mean-Based Preconditioner:\n") +
     std::string("              ***** ") +
     std::string(mean_prec->Label());
   num_blocks = sg_basis()->size();
}

int
Stokhos::MeanBasedPreconditioner::
SetUseTranspose(bool UseTheTranspose)
{
  useTranspose = UseTheTranspose;
  mean_prec->SetUseTranspose(useTranspose);

  return 0;
}

int
Stokhos::MeanBasedPreconditioner::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  int myBlockRows = epetraCijk->numMyRows();

  if (!use_block_apply) {
    EpetraExt::BlockMultiVector sg_input(View, *base_map, Input);
    EpetraExt::BlockMultiVector sg_result(View, *base_map, Result);
    for (int i=0; i<myBlockRows; i++) {
      mean_prec->Apply(*(sg_input.GetBlock(i)), *(sg_result.GetBlock(i)));
    }
  }

  else {
    int m = Input.NumVectors();
    Epetra_MultiVector input_block(
      View, *base_map, Input.Values(), base_map->NumMyElements(),
      m*myBlockRows);
    Epetra_MultiVector result_block(
      View, *base_map, Result.Values(), base_map->NumMyElements(),
      m*myBlockRows);
    mean_prec->Apply(input_block, result_block);
  }

  return 0;
}

int
Stokhos::MeanBasedPreconditioner::
ApplyInverse(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  int myBlockRows = epetraCijk->numMyRows();

  if (!use_block_apply) {
    EpetraExt::BlockMultiVector sg_input(View, *base_map, Input);
    EpetraExt::BlockMultiVector sg_result(View, *base_map, Result);
    for (int i=0; i<myBlockRows; i++) {
      mean_prec->ApplyInverse(*(sg_input.GetBlock(i)),
                              *(sg_result.GetBlock(i)));
    }
  }

  else {
    int m = Input.NumVectors();
    Epetra_MultiVector input_block(
      View, *base_map, Input.Values(), base_map->NumMyElements(),
      m*myBlockRows);
    Epetra_MultiVector result_block(
      View, *base_map, Result.Values(), base_map->NumMyElements(),
      m*myBlockRows);
    mean_prec->ApplyInverse(input_block, result_block);
  }

  return 0;
}

double
Stokhos::MeanBasedPreconditioner::
NormInf() const
{
  return mean_prec->NormInf();
}


const char*
Stokhos::MeanBasedPreconditioner::
Label () const
{
  return const_cast<char*>(label.c_str());
}

bool
Stokhos::MeanBasedPreconditioner::
UseTranspose() const
{
  return useTranspose;
}

bool
Stokhos::MeanBasedPreconditioner::
HasNormInf() const
{
  return true;
}

const Epetra_Comm &
Stokhos::MeanBasedPreconditioner::
Comm() const
{
  return *sg_comm;
}
const Epetra_Map&
Stokhos::MeanBasedPreconditioner::
OperatorDomainMap() const
{
  return *sg_map;
}

const Epetra_Map&
Stokhos::MeanBasedPreconditioner::
OperatorRangeMap() const
{
  return *sg_map;
}
