// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ProblemLOCAPrototype.H"

ProblemLOCAPrototype::ProblemLOCAPrototype()
{
}

ProblemLOCAPrototype::~ProblemLOCAPrototype()
{
}


bool ProblemLOCAPrototype::
PrintSolutionFile( const std::string & /* fileName */, const Epetra_Vector & /* x */,
  const Teuchos::ParameterList & /* xParams */)
{
  return true;
}

bool ProblemLOCAPrototype::
SetSolutionFileParameters(const Epetra_Vector & /* x */)
{
  return true;
}

Teuchos::RCP <Teuchos::ParameterList> ProblemLOCAPrototype::
GetSolutionFileParameters()
{
  // Returns an empty list
  Teuchos::RCP <Teuchos::ParameterList> emptyList =
    Teuchos::rcp (new Teuchos::ParameterList());

  return emptyList;
}

bool ProblemLOCAPrototype::
SetContinuationFileParameters(const Epetra_Vector & /* x */)
{
  return true;
}

Teuchos::RCP <Teuchos::ParameterList> ProblemLOCAPrototype::
GetContinuationFileParameters()
{
  // Returns an empty list
  Teuchos::RCP <Teuchos::ParameterList> emptyList =
    Teuchos::rcp (new Teuchos::ParameterList());

  return emptyList;
}

bool ProblemLOCAPrototype::
UpdateContinuationFile( const std::string & /* fileName */,
  const int & /* idStep */,
  const Teuchos::ParameterList & /* continuationFileParams */)
{
  return true;
}

bool ProblemLOCAPrototype::
ComputePeriodicDirectionDerivative(const Epetra_Vector & /* x */,
    Epetra_Vector & /* dx */)
{
  return true;
}
