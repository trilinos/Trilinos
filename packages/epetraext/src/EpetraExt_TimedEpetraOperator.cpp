//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "Epetra_config.h"
#include "EpetraExt_TimedEpetraOperator.hpp"

EpetraExt::Epetra_Timed_Operator::Epetra_Timed_Operator(const Teuchos::RCP<Epetra_Operator>& A_) 
  : A(A_)
{
ApplyTimer = Teuchos::rcp(new Teuchos::Time("apply timer",false));
ApplyInverseTimer = Teuchos::rcp(new Teuchos::Time("apply inverse timer",false));
}

EpetraExt::Epetra_Timed_Operator::~Epetra_Timed_Operator()
{
}

int 
EpetraExt::Epetra_Timed_Operator::SetUseTranspose(bool useTranspose) 
{
  int success;
  success = A->SetUseTranspose(useTranspose);
  return success;
}

int 
EpetraExt::Epetra_Timed_Operator::Apply(const Epetra_MultiVector& Input, 
				   Epetra_MultiVector& Result) const
{
  int success;
  ApplyTimer->start();
  success = A->Apply(Input,Result);
  ApplyTimer->stop();
  return success;
}

int 
EpetraExt::Epetra_Timed_Operator::ApplyInverse(const Epetra_MultiVector& Input, 
					  Epetra_MultiVector& Result) const
{
  int success;
  ApplyInverseTimer->start();
  success = A->ApplyInverse(Input,Result);
  ApplyInverseTimer->stop();
  return success;
}

double 
EpetraExt::Epetra_Timed_Operator::NormInf() const
{
  return A->NormInf();
}


const char* 
EpetraExt::Epetra_Timed_Operator::Label () const
{
  return A->Label();
}
  
bool 
EpetraExt::Epetra_Timed_Operator::UseTranspose() const
{
  return A->UseTranspose();
}

bool 
EpetraExt::Epetra_Timed_Operator::HasNormInf() const
{
  return A->HasNormInf();
}

const Epetra_Comm & 
EpetraExt::Epetra_Timed_Operator::Comm() const
{
  return A->Comm();
}
const Epetra_Map& 
EpetraExt::Epetra_Timed_Operator::OperatorDomainMap() const
{
  return A->OperatorDomainMap();
}

const Epetra_Map& 
EpetraExt::Epetra_Timed_Operator::OperatorRangeMap() const
{
  return A->OperatorRangeMap();
}
