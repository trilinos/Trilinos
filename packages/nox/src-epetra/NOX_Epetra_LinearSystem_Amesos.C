// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#include "NOX_Epetra_LinearSystem_Amesos.H"

#ifdef HAVE_NOX_AMESOS

#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "NOX_Epetra_Interface_Jacobian.H"

NOX::Epetra::LinearSystemAmesos::
LinearSystemAmesos(
  Teuchos::ParameterList& printingParams,
  Teuchos::ParameterList& linearSolverParams,
  const Teuchos::RCP<NOX::Epetra::Interface::Required>& /* iReq */,
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
      "Invalid Amesos Solver: " << linearSolverParams.get<std::string>("Amesos Solver"));
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
applyJacobianInverse(Teuchos::ParameterList &/* params */,
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
applyRightPreconditioning(bool /* useTranspose */,
                        Teuchos::ParameterList& /* params */,
                        const NOX::Epetra::Vector& /* input */,
                        NOX::Epetra::Vector& /* result */) const
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
createPreconditioner(const NOX::Epetra::Vector& /* x */,
                      Teuchos::ParameterList& /* p */,
                      bool /* recomputeGraph */) const
{
  return false;
}

bool NOX::Epetra::LinearSystemAmesos::
destroyPreconditioner() const
{
  return false;
}

bool NOX::Epetra::LinearSystemAmesos::
recomputePreconditioner(const NOX::Epetra::Vector& /* x */,
              Teuchos::ParameterList& /* linearSolverParams */) const
{
  return false;
}

NOX::Epetra::LinearSystem::PreconditionerReusePolicyType
NOX::Epetra::LinearSystemAmesos::
getPreconditionerPolicy(bool /* advanceReuseCounter */)
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
           Teuchos::RCP<const Epetra_Operator>& /* solvePrecOp */)
{
  return;
}

#endif //HAVE_NOX_AMESOS
