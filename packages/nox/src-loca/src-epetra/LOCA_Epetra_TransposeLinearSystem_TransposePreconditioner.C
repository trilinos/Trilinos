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

#include "LOCA_Epetra_TransposeLinearSystem_TransposePreconditioner.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Epetra_LinearSystem.H"
#include "NOX_Epetra_Scaling.H"
#include "Epetra_Operator.h"

LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::
TransposePreconditioner(
	     const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	     const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams,
	     const Teuchos::RefCountPtr<NOX::Epetra::LinearSystem>& linsys_) :
  globalData(global_data),
  linsys(linsys_),
  jac(),
  prec(),
  scaling_trans()
{
  // Get transpose scaling object
  if (solverParams->isParameter("Transpose Scaling"))
    scaling_trans = (*solverParams).INVALID_TEMPLATE_QUALIFIER 
      get< Teuchos::RefCountPtr<NOX::Epetra::Scaling> >("Transpose Scaling");
}

LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::
~TransposePreconditioner()
{
}

bool
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::
applyJacobianTransposeInverse(Teuchos::ParameterList &params, 
			      const NOX::Epetra::Vector &input, 
			      NOX::Epetra::Vector &result)
{  
  // Replace operators with transposed operators
  linsys->setJacobianOperatorForSolve(jac);
  if (linsys->hasPreconditioner())
    linsys->setPrecOperatorForSolve(prec);

  // Set the transpose scaling object if we have one
  Teuchos::RefCountPtr<NOX::Epetra::Scaling> scaling_orig;
  if (scaling_trans != Teuchos::null) {
    scaling_orig = linsys->getScaling();
    linsys->resetScaling(scaling_trans);
  }

  // Solve the system
  bool res = linsys->applyJacobianInverse(params, input, result);

  // Set the original scaling object back in the linear system
  if (scaling_trans != Teuchos::null)
    linsys->resetScaling(scaling_orig);

  return res;
}

bool
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::
computeJacobianTranspose(const NOX::Epetra::Vector& x)
{
  // Compute the Jacobian and set it to use the transpose
  bool res = linsys->computeJacobian(x);
  jac = linsys->getJacobianOperator();
  jac->SetUseTranspose(true);

  return res;
}

bool
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::
createTransposePreconditioner(const NOX::Epetra::Vector& x, 
			      Teuchos::ParameterList& p)
{
  // We're done if the linear system doesn't define a preconditioner
  if (!linsys->hasPreconditioner())
    return true;

  // Destroy any previous preconditioner
  bool res1 = linsys->destroyPreconditioner();

  // Set the transpose scaling object if we have one
  Teuchos::RefCountPtr<NOX::Epetra::Scaling> scaling_orig;
  if (scaling_trans != Teuchos::null) {
    scaling_orig = linsys->getScaling();
    linsys->resetScaling(scaling_trans);
  }

  // Compute the preconditioner and set it to use the transpose
  linsys->setJacobianOperatorForSolve(jac);
  bool res2 = linsys->createPreconditioner(x, p, true);
  prec = linsys->getGeneratedPrecOperator();
  prec->SetUseTranspose(true);

  // Set the original scaling object back in the linear system
  if (scaling_trans != Teuchos::null)
    linsys->resetScaling(scaling_orig);

  return res1 && res2;
}

Teuchos::RefCountPtr<Epetra_Operator> 
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::
getJacobianTransposeOperator()
{
  return jac;
}

Teuchos::RefCountPtr<Epetra_Operator> 
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::
getTransposePreconditioner()
{
  return prec;
}

void
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::
setJacobianTransposeOperator(
	       const Teuchos::RefCountPtr<Epetra_Operator>& new_jac_trans)
{
  jac = new_jac_trans;
}

void
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::
setTransposePreconditioner(
	      const Teuchos::RefCountPtr<Epetra_Operator>& new_prec_trans)
{
  prec = new_prec_trans;
}
