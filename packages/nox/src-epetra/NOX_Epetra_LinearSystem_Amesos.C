/*
//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
*/


#include "NOX_Epetra_LinearSystem_Amesos.H"

#ifdef HAVE_NOX_AMESOS

#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "NOX_Epetra_Interface_Jacobian.H"

NOX::Epetra::LinearSystemAmesos::
LinearSystemAmesos(
  Teuchos::ParameterList& printingParams, 
  Teuchos::ParameterList& linearSolverParams, 
  const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
  const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac, 
  const Teuchos::RCP<Epetra_Operator>& J,
  const NOX::Epetra::Vector& cloneVector,
  const Teuchos::RCP<NOX::Epetra::Scaling> s):
  amesosProblem(Teuchos::null),
  amesosSolver(Teuchos::null),
  factory(),
  isValidFactorization(false),
  jacInterfacePtr(iJac),
  jacPtr(J),
  leftHandSide(Teuchos::rcp(new Epetra_Vector(cloneVector.getEpetraVector()))),
  rightHandSide(Teuchos::rcp(new Epetra_Vector(cloneVector.getEpetraVector()))),
  scaling(s),
  timer(cloneVector.getEpetraVector().Comm()),
  utils(printingParams)
{
  amesosProblem = Teuchos::rcp(new Epetra_LinearProblem(
				      dynamic_cast<Epetra_CrsMatrix *>(jacPtr.get()),
				      leftHandSide.get(),
				      rightHandSide.get()));

  Amesos_BaseSolver * tmp = factory.Create(linearSolverParams.get("Amesos Solver","Amesos_Klu"), 
      *amesosProblem);
  TEUCHOS_TEST_FOR_EXCEPTION ( tmp == 0, Teuchos::Exceptions::InvalidParameterValue, 
      "Invalid Amesos Solver: " << linearSolverParams.get<string>("Amesos Solver"));
  amesosSolver = Teuchos::rcp(tmp);

  amesosSolver->SetParameters(linearSolverParams);
}

NOX::Epetra::LinearSystemAmesos::
~LinearSystemAmesos()
{
}

bool NOX::Epetra::LinearSystemAmesos::
applyJacobian(const NOX::Epetra::Vector& input, 
      		     NOX::Epetra::Vector& result) const
{
  jacPtr->SetUseTranspose(false);
  int status = jacPtr->Apply(input.getEpetraVector(), 
				  result.getEpetraVector());
  return (status == 0);
}

bool NOX::Epetra::LinearSystemAmesos::
applyJacobianTranspose(const NOX::Epetra::Vector& input, 
      			      NOX::Epetra::Vector& result) const
{
  jacPtr->SetUseTranspose(true);
  int status = jacPtr->Apply(input.getEpetraVector(), 
				  result.getEpetraVector());
  jacPtr->SetUseTranspose(false);

  return (status == 0);
}

bool NOX::Epetra::LinearSystemAmesos::
applyJacobianInverse(Teuchos::ParameterList &params, 
      			    const NOX::Epetra::Vector &input, 
      			    NOX::Epetra::Vector &result)
{
  double startTime = timer.WallTime();

  *rightHandSide = input.getEpetraVector();
  int err = 0;
  bool status = true;
  if (!isValidFactorization) 
  {
    err = amesosSolver->SymbolicFactorization();
    if (err > 0) 
      status = false;

    err = amesosSolver->NumericFactorization();
    if (err > 0) 
      status = false;

    if (status) isValidFactorization = true;

  }

  err = amesosSolver->Solve();
  if (err > 0) status = false;

  result.getEpetraVector() = *leftHandSide;

  double endTime = timer.WallTime();
  if (utils.isPrintType(Utils::LinearSolverDetails))
    utils.out() << "\n       Time required for one linear solve : " 
         << (endTime - startTime) << " (sec.)" << std::endl;;

  return status;
}

bool NOX::Epetra::LinearSystemAmesos::
applyRightPreconditioning(bool useTranspose,
      			      Teuchos::ParameterList& params, 
      			      const NOX::Epetra::Vector& input, 
      			      NOX::Epetra::Vector& result) const
{
  return false;
}

Teuchos::RCP<NOX::Epetra::Scaling> NOX::Epetra::LinearSystemAmesos::
getScaling()
{
  return scaling;
}

void NOX::Epetra::LinearSystemAmesos::
resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& s)
{
  scaling = s;
  return;
}

bool NOX::Epetra::LinearSystemAmesos::
computeJacobian(const NOX::Epetra::Vector& x)
{
  bool success = jacInterfacePtr->computeJacobian(x.getEpetraVector(), 
						  *jacPtr);
  if (success) isValidFactorization = false;
  return success;
}

bool NOX::Epetra::LinearSystemAmesos::
createPreconditioner(const NOX::Epetra::Vector& x, 
      			    Teuchos::ParameterList& p,
      			    bool recomputeGraph) const
{
  return false;
}

bool NOX::Epetra::LinearSystemAmesos::
destroyPreconditioner() const
{
  return false;
}

bool NOX::Epetra::LinearSystemAmesos::
recomputePreconditioner(const NOX::Epetra::Vector& x, 
      		Teuchos::ParameterList& linearSolverParams) const
{
  return false;
}

NOX::Epetra::LinearSystem::PreconditionerReusePolicyType 
NOX::Epetra::LinearSystemAmesos::
getPreconditionerPolicy(bool advanceReuseCounter)
{
  return PRPT_REUSE;
}

bool NOX::Epetra::LinearSystemAmesos::
isPreconditionerConstructed() const
{
  return false;
}

bool NOX::Epetra::LinearSystemAmesos::
hasPreconditioner() const
{
  return false;
}

Teuchos::RCP<const Epetra_Operator> NOX::Epetra::LinearSystemAmesos::
getJacobianOperator() const
{
  return jacPtr;
}

Teuchos::RCP<Epetra_Operator> NOX::Epetra::LinearSystemAmesos::
getJacobianOperator()
{
  return jacPtr;
}

Teuchos::RCP<const Epetra_Operator> NOX::Epetra::LinearSystemAmesos::
getGeneratedPrecOperator() const
{
  return Teuchos::null;
}

Teuchos::RCP<Epetra_Operator> NOX::Epetra::LinearSystemAmesos::
getGeneratedPrecOperator()
{
  return Teuchos::null;
}

void NOX::Epetra::LinearSystemAmesos::
setJacobianOperatorForSolve(const 
      	 Teuchos::RCP<const Epetra_Operator>& solveJacOp)
{
  jacPtr = Teuchos::rcp_const_cast<Epetra_Operator>(solveJacOp);
  isValidFactorization = false;
  return;
}

void NOX::Epetra::LinearSystemAmesos::
setPrecOperatorForSolve(const 
      	 Teuchos::RCP<const Epetra_Operator>& solvePrecOp)
{
  return;
}

#endif //HAVE_NOX_AMESOS
