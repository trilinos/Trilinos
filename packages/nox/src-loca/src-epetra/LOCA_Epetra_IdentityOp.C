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
