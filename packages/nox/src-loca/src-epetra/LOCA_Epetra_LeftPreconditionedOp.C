// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Epetra_config.h"
#include "Epetra_MultiVector.h"
#include "LOCA_Epetra_LeftPreconditionedOp.H"

LOCA::Epetra::LeftPreconditionedOp::LeftPreconditionedOp(
    const Teuchos::RCP<Epetra_Operator>& jacOperator,
    const Teuchos::RCP<Epetra_Operator>& precOperator) :
  label("LOCA::Epetra::LeftPreconditionedOp"),
  J(jacOperator),
  M(precOperator),
  useTranspose(false)
{
}

LOCA::Epetra::LeftPreconditionedOp::~LeftPreconditionedOp()
{
}

int
LOCA::Epetra::LeftPreconditionedOp::SetUseTranspose(bool UseTranspose)
{
  useTranspose = UseTranspose;
  int res_1 = J->SetUseTranspose(UseTranspose);
  int res_2 = M->SetUseTranspose(UseTranspose);

  return res_1 + res_2;
}

int
LOCA::Epetra::LeftPreconditionedOp::Apply(const Epetra_MultiVector& Input,
                      Epetra_MultiVector& Result) const
{
  // Create temporary multivector
  Epetra_MultiVector tmp(Input);
  int res_1, res_2;

  if (!useTranspose) {

    // Compute J*Input
    res_1 = J->Apply(Input, tmp);

    // Compute M^-1*J*Input
    res_2 = M->ApplyInverse(tmp, Result);

  }
  else {

    // Compute M^-T*Input
    res_1 = M->ApplyInverse(Input, tmp);

    // Compute J^T*M^-T*Input
    res_2 = J->Apply(tmp, Result);

  }

  return res_1 + res_2;
}

int
LOCA::Epetra::LeftPreconditionedOp::ApplyInverse(
                    const Epetra_MultiVector& Input,
                    Epetra_MultiVector& Result) const
{
  // Create temporary multivector
  Epetra_MultiVector tmp(Input);
  int res_1, res_2;

  if (!useTranspose) {

    // Compute M*Input
    res_1 = M->Apply(Input, tmp);

    // Compute J^-1*M*Input
    res_2 = J->ApplyInverse(tmp, Result);

  }
  else {

    // Compute J^-T*Input
    res_1 = J->ApplyInverse(Input, tmp);

    // Compute M^T*J^-T*Input
    res_2 = M->Apply(tmp, Result);

  }

  return res_1 + res_2;
}

double
LOCA::Epetra::LeftPreconditionedOp::NormInf() const
{
  return J->NormInf() + 1.0/M->NormInf();
}


const char*
LOCA::Epetra::LeftPreconditionedOp::Label () const
{
  return const_cast<char*>(label.c_str());
}

bool
LOCA::Epetra::LeftPreconditionedOp::UseTranspose() const
{
  return useTranspose;
}

bool
LOCA::Epetra::LeftPreconditionedOp::HasNormInf() const
{
  return J->HasNormInf() && M->HasNormInf();
}

const Epetra_Comm &
LOCA::Epetra::LeftPreconditionedOp::Comm() const
{
  return J->Comm();
}
const Epetra_Map&
LOCA::Epetra::LeftPreconditionedOp::OperatorDomainMap() const
{
  return J->OperatorDomainMap();
}

const Epetra_Map&
LOCA::Epetra::LeftPreconditionedOp::OperatorRangeMap() const
{
  return M->OperatorDomainMap();
}
