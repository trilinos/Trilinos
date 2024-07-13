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
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "LOCA_Epetra_IdentityOp.H"

LOCA::Epetra::IdentityOp::IdentityOp(
             const Teuchos::RCP<const Epetra_Comm>& comm_,
             const Teuchos::RCP<const Epetra_Map>& map_) :
  label("LOCA::Epetra::IdentityOp"),
  comm(comm_),
  map(map_),
  useTranspose(false)
{
}

LOCA::Epetra::IdentityOp::~IdentityOp()
{
}

int
LOCA::Epetra::IdentityOp::SetUseTranspose(bool UseTranspose)
{
  useTranspose = UseTranspose;
  return 0;
}

int
LOCA::Epetra::IdentityOp::Apply(const Epetra_MultiVector& Input,
                Epetra_MultiVector& Result) const
{
  Result.Scale(1.0, Input);

  return 0;
}

int
LOCA::Epetra::IdentityOp::ApplyInverse(const Epetra_MultiVector& Input,
                       Epetra_MultiVector& Result) const
{
  Result.Scale(1.0, Input);

  return 0;
}

double
LOCA::Epetra::IdentityOp::NormInf() const
{
  return 1.0;
}


const char*
LOCA::Epetra::IdentityOp::Label () const
{
  return const_cast<char*>(label.c_str());
}

bool
LOCA::Epetra::IdentityOp::UseTranspose() const
{
  return useTranspose;
}

bool
LOCA::Epetra::IdentityOp::HasNormInf() const
{
  return true;
}

const Epetra_Comm &
LOCA::Epetra::IdentityOp::Comm() const
{
  return *comm;
}
const Epetra_Map&
LOCA::Epetra::IdentityOp::OperatorDomainMap() const
{
  return *map;
}

const Epetra_Map&
LOCA::Epetra::IdentityOp::OperatorRangeMap() const
{
  return *map;
}
