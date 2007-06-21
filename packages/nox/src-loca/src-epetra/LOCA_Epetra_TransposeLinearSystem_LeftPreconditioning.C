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

#include "LOCA_Epetra_TransposeLinearSystem_LeftPreconditioning.H"
#include "NOX_Epetra_LinearSystem.H"
#include "NOX_Epetra_Vector.H"
#include "LOCA_Epetra_LeftPreconditionedOp.H"
#include "LOCA_Epetra_IdentityOp.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"

LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::
LeftPreconditioning(
	     const Teuchos::RCP<LOCA::GlobalData>& global_data,
	     const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
	     const Teuchos::RCP<NOX::Epetra::LinearSystem>& linsys_) :
  globalData(global_data),
  linsys(linsys_),
  jac(),
  prec()
{
  // This strategy makes sense only if the linear system defines a 
  // preconditioner
  if (!linsys->hasPreconditioner())
    globalData->locaErrorCheck->throwError(
	  string("LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning") + 
	  string("::LeftPreconditioning()"),
	  string("Left preconditioning for transpose solve is valid only") + 
	  string(" when the linear system defines a preconditioner"));
}

LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::
~LeftPreconditioning()
{
}

bool
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::
applyJacobianTransposeInverse(Teuchos::ParameterList &params, 
			      const NOX::Epetra::Vector &input, 
			      NOX::Epetra::Vector &result)
{  
  // Create preconditioned operator
  Teuchos::RCP<Epetra_Operator> left_prec_jac = 
    Teuchos::rcp(new LOCA::Epetra::LeftPreconditionedOp(jac, prec));

  // Replace Jacobian operator with transposed, left preconditioned op
  linsys->setJacobianOperatorForSolve(left_prec_jac);

  // Create identity operator as a right preconditioner
  Teuchos::RCP<Epetra_Operator> identity_prec = 
    Teuchos::rcp(new LOCA::Epetra::IdentityOp(
			   Teuchos::rcp(&(jac->Comm()), false), 
			   Teuchos::rcp(&(jac->OperatorDomainMap()), false)));

  // Replace preconditioner with identity
  linsys->setPrecOperatorForSolve(identity_prec);

  // Precondition the RHS
  NOX::Epetra::Vector prec_input(input);
  prec->ApplyInverse(input.getEpetraVector(), prec_input.getEpetraVector());

  // Solve the system
  bool res = linsys->applyJacobianInverse(params, prec_input, result);

  return res;
}

bool
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::
createJacobianTranspose()
{
  // Set Jacobian to use the transpose
  jac = linsys->getJacobianOperator();
  jac->SetUseTranspose(true);

  return true;
}

bool
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::
createTransposePreconditioner(const NOX::Epetra::Vector& x, 
			      Teuchos::ParameterList& p)
{
  // Compute the preconditioner and set it to use the transpose
  bool res1 = linsys->destroyPreconditioner();
  linsys->setJacobianOperatorForSolve(jac);
  bool res2 = linsys->createPreconditioner(x, p, true);
  prec = linsys->getGeneratedPrecOperator();
  prec->SetUseTranspose(true);

  return res1 && res2;
}

Teuchos::RCP<Epetra_Operator> 
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::
getJacobianTransposeOperator()
{
  return jac;
}

Teuchos::RCP<Epetra_Operator> 
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::
getTransposePreconditioner()
{
  return prec;
}

void
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::
setJacobianTransposeOperator(
	       const Teuchos::RCP<Epetra_Operator>& new_jac_trans)
{
  jac = new_jac_trans;
}

void
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::
setTransposePreconditioner(
	      const Teuchos::RCP<Epetra_Operator>& new_prec_trans)
{
  prec = new_prec_trans;
}
