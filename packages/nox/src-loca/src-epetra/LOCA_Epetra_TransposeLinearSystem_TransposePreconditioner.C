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

#include "LOCA_Epetra_TransposeLinearSystem_TransposePreconditioner.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Epetra_LinearSystem.H"
#include "NOX_Epetra_Scaling.H"
#include "Epetra_Operator.h"

LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::
TransposePreconditioner(
         const Teuchos::RCP<LOCA::GlobalData>& global_data,
         const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
         const Teuchos::RCP<NOX::Epetra::LinearSystem>& linsys_) :
  globalData(global_data),
  linsys(linsys_),
  jac(),
  prec(),
  scaling_trans()
{
  // Get transpose scaling object
  if (solverParams->isParameter("Transpose Scaling"))
    scaling_trans = (*solverParams).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RCP<NOX::Epetra::Scaling> >("Transpose Scaling");
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
  Teuchos::RCP<NOX::Epetra::Scaling> scaling_orig;
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
createJacobianTranspose()
{
  // Set Jacobian to use the transpose
  jac = linsys->getJacobianOperator();
  jac->SetUseTranspose(true);

  return true;
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
  Teuchos::RCP<NOX::Epetra::Scaling> scaling_orig;
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

Teuchos::RCP<Epetra_Operator>
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::
getJacobianTransposeOperator()
{
  return jac;
}

Teuchos::RCP<Epetra_Operator>
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::
getTransposePreconditioner()
{
  return prec;
}

void
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::
setJacobianTransposeOperator(
           const Teuchos::RCP<Epetra_Operator>& new_jac_trans)
{
  jac = new_jac_trans;
}

void
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::
setTransposePreconditioner(
          const Teuchos::RCP<Epetra_Operator>& new_prec_trans)
{
  prec = new_prec_trans;
}
