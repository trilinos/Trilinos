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
	  std::string("LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning") + 
	  std::string("::LeftPreconditioning()"),
	  std::string("Left preconditioning for transpose solve is valid only") + 
	  std::string(" when the linear system defines a preconditioner"));
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
