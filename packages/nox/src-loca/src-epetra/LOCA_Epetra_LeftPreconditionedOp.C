// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

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
