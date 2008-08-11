// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
