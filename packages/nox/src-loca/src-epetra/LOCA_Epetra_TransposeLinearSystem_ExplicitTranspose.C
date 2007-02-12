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

#include "LOCA_Epetra_TransposeLinearSystem_ExplicitTranspose.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Epetra_LinearSystem.H"
#include "NOX_Epetra_Scaling.H"
#include "Epetra_RowMatrix.h"

LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::
ExplicitTranspose(
	     const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	     const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams,
	     const Teuchos::RefCountPtr<NOX::Epetra::LinearSystem>& linsys_) :
  globalData(global_data),
  linsys(linsys_),
  jac_trans(),
  prec_trans(),
  scaling_trans(),
  transposer(true) // make data contiguous
{
  // Get transpose scaling object
  if (solverParams->isParameter("Transpose Scaling"))
    scaling_trans = (*solverParams).INVALID_TEMPLATE_QUALIFIER 
      get< Teuchos::RefCountPtr<NOX::Epetra::Scaling> >("Transpose Scaling");
}

LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::
~ExplicitTranspose()
{
}

bool
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::
applyJacobianTransposeInverse(Teuchos::ParameterList &params, 
			      const NOX::Epetra::Vector &input, 
			      NOX::Epetra::Vector &result)
{  
  // Replace operators with transposed operators
  linsys->setJacobianOperatorForSolve(jac_trans);
  if (linsys->hasPreconditioner())
    linsys->setPrecOperatorForSolve(prec_trans);

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
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::
createJacobianTranspose()
{
  // Get Jacobian
  Teuchos::RefCountPtr<Epetra_RowMatrix> jac = 
    Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(linsys->getJacobianOperator());

  if (jac == Teuchos::null)
    globalData->locaErrorCheck->throwError(
	 string("LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::") +
	 string("createJacobianTranspose()"), 
	 string("Jacobian operator must be of type Epetra_RowMatrix for ") +
	 string("Explicit Transpose method"));

  // Form transpose if we haven't already, otherwise just migrate data
  if (jac_trans == Teuchos::null)
    jac_trans = Teuchos::rcp(&(transposer(*jac)), false);
  else
    transposer.fwd();

  return true;
}

bool
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::
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
  
  // Set Jacobian-transpose operator and compute preconditioner
  linsys->setJacobianOperatorForSolve(jac_trans);
  bool res2 = linsys->createPreconditioner(x, p, true);
  prec_trans = linsys->getGeneratedPrecOperator();

  // Set the original scaling object back in the linear system
  if (scaling_trans != Teuchos::null)
    linsys->resetScaling(scaling_orig);

  return res1 && res2;
}

Teuchos::RefCountPtr<Epetra_Operator> 
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::
getJacobianTransposeOperator()
{
  return jac_trans;
}

Teuchos::RefCountPtr<Epetra_Operator> 
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::
getTransposePreconditioner()
{
  return prec_trans;
}

void
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::
setJacobianTransposeOperator(
	       const Teuchos::RefCountPtr<Epetra_Operator>& new_jac_trans)
{
  jac_trans = new_jac_trans;
}

void
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::
setTransposePreconditioner(
	      const Teuchos::RefCountPtr<Epetra_Operator>& new_prec_trans)
{
  prec_trans = new_prec_trans;
}
