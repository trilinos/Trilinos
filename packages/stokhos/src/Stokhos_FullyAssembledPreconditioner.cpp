// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
    std::string("		***** ") + 
    std::string(prec->Label());
}

int 
Stokhos::FullyAssembledPreconditioner::
SetUseTranspose(bool UseTranspose) 
{
  return prec->SetUseTranspose(UseTranspose);
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
