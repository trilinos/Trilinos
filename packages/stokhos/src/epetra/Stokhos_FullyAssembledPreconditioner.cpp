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

#include "Stokhos_FullyAssembledPreconditioner.hpp"
#include "Stokhos_FullyAssembledOperator.hpp"

Stokhos::FullyAssembledPreconditioner::
FullyAssembledPreconditioner(
   const Teuchos::RCP<Stokhos::AbstractPreconditionerFactory>& prec_factory_,
   const Teuchos::RCP<Teuchos::ParameterList>& params_) :
  label("Stokhos Fully Assembled Preconditioner"),
  prec_factory(prec_factory_),
  prec()
{
}

Stokhos::FullyAssembledPreconditioner::
~FullyAssembledPreconditioner()
{
}

void
Stokhos::FullyAssembledPreconditioner::
setupPreconditioner(const Teuchos::RCP<Stokhos::SGOperator>& sg_op,
                    const Epetra_Vector& x)
{

  Teuchos::RCP<Stokhos::FullyAssembledOperator > fa_op =
    Teuchos::rcp_dynamic_cast<Stokhos::FullyAssembledOperator>(sg_op, true);
  prec = prec_factory->compute(Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(fa_op));
  label = std::string("Stokhos Fully Assembled Preconditioner:\n") +
    std::string("               ***** ") +
    std::string(prec->Label());
}

int
Stokhos::FullyAssembledPreconditioner::
SetUseTranspose(bool UseTheTranspose)
{
  return prec->SetUseTranspose(UseTheTranspose);
}

int
Stokhos::FullyAssembledPreconditioner::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  return prec->Apply(Input, Result);
}

int
Stokhos::FullyAssembledPreconditioner::
ApplyInverse(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  return prec->ApplyInverse(Input, Result);
}

double
Stokhos::FullyAssembledPreconditioner::
NormInf() const
{
  return prec->NormInf();
}


const char*
Stokhos::FullyAssembledPreconditioner::
Label () const
{
  return const_cast<char*>(label.c_str());
}

bool
Stokhos::FullyAssembledPreconditioner::
UseTranspose() const
{
  return prec->UseTranspose();
}

bool
Stokhos::FullyAssembledPreconditioner::
HasNormInf() const
{
  return prec->HasNormInf();
}

const Epetra_Comm &
Stokhos::FullyAssembledPreconditioner::
Comm() const
{
  return prec->Comm();
}
const Epetra_Map&
Stokhos::FullyAssembledPreconditioner::
OperatorDomainMap() const
{
  return prec->OperatorDomainMap();
}

const Epetra_Map&
Stokhos::FullyAssembledPreconditioner::
OperatorRangeMap() const
{
  return prec->OperatorRangeMap();
}
