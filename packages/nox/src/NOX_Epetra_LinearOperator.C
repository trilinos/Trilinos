// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

// Different RowMatrix implementations
//#include "NOX_Epetra_RowMatrix_ExplicitJacobian.H"

#include <string>
#include "NOX_Epetra_Group.H"

#include "NOX_Epetra_LinearOperator.H" 
#include "NOX_Parameter_List.H"
#include "NOX_Epetra_SharedJacobian.H"
#include "NOX_Epetra_MatrixFree.H"
#include "NOX_Epetra_FiniteDifference.H"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"


using namespace NOX;
using namespace NOX::Epetra;

//Constructor  - Jacobian ONLY!
LinearOperator::LinearOperator(const Parameter::List& p, Interface& i, Epetra_RowMatrix& j) :
  precType(p.getParameter("Preconditioning Matrix Type","Use Jacobian")),
  userInterface(i),
  jacobian(j)
{
  // Set the Jacobian matrix type
  jacType = getJacobianType();

  // Create the shared jacobian
  sharedJacobian = new SharedJacobian(j);

  // Set the preconditioner if required
  if (precType == "None") {
    prec = 0;
  }
  else if (precType == "Use Jacobian") {
    prec = &jacobian;
  }
  else if (precType == "User Supplied Action on Vector") {
    prec = 0;
  }
  else {
    cout << "ERROR:NOX::Epetra::LinearOperator - Preconditioning Matrix Type "
	 << "requested is invalid for this constructor!" << endl;
    throw 1;
  }
}

// Constructor
LinearOperator::LinearOperator(const Parameter::List& p, Interface& i, 
			       Epetra_RowMatrix& j, Epetra_RowMatrix& m) :
  precType(p.getParameter("Preconditioning Matrix Type", 
			  "User Supplied Matrix")),
  userInterface(i),
  jacobian(j),
  prec(&m)
{
  // Set the Jacobian matrix type
  jacType = getJacobianType();

  sharedJacobian = new SharedJacobian(j,m);
}


// Destructor
LinearOperator::~LinearOperator()
{
  delete sharedJacobian;
}

void LinearOperator::setAztecOptions(const Parameter::List& p, AztecOO& aztec)
{

  // Preconditioning Matrix Type 
    if (p.isParameterEqual("Preconditioning Matrix Type", "None")) 
      aztec.SetAztecOption(AZ_precond, AZ_none);
     
  // Preconditioning options
  if (p.isParameter("Preconditioning")) {
    
    string prec;
    p.getParameter("Preconditioning",prec);
    
    if (prec == "None") 
      aztec.SetAztecOption(AZ_precond, AZ_none);
    
    else if (prec == "ilut") 
      aztec.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
    
    else if (prec == "Polynomial") {
      aztec.SetAztecOption(AZ_precond,AZ_Neumann);
      aztec.SetAztecOption(AZ_poly_ord,3);
    }
    else {
      cout << "ERROR: NOX::Epetra::LinearOperator::setAztecOptions" << endl
	   << "Preconditioning flag is invalid!" << endl;
      throw;
    }
  }

  // Frequency of linear solve residual output
  if (p.isParameter("Iteration Output Frequency")) {
    int i = p.getParameter("Iteration Output Frequency", AZ_last);
    aztec.SetAztecOption(AZ_output, i);
  }
}
// derived  
bool LinearOperator::solveLinearSystem(const Parameter::List& params, 
				       Group* g)
{
 
  // Get the Jacobian 
  /* (Have to get non-const version which requires reasserting
     ownership. This is due to a flaw in Epetra_LinearProblem). */
  Epetra_RowMatrix& Jacobian = g->sharedJacobian.getJacobian(g);

  // Create Epetra problem for the linear solve
  Epetra_LinearProblem Problem(&Jacobian, 
			       &(g->NewtonVector.getEpetraVector()), 
			       &(g->RHSVector.getEpetraVector()));

  // Set the default Problem parameters to "hard" (this sets Aztec deafuats
  // during the AztecOO instantiation)
  Problem.SetPDL(hard);

  // Create aztec problem
  AztecOO aztec(Problem);

  // Set specific Aztec parameters requested by NOX
  setAztecOptions(params, aztec);

  int maxit = params.getParameter("Max Iterations", 400);
  double tol = params.getParameter("Tolerance", 1.0e-6);

  // Solve Aztex problem
  aztec.Iterate(maxit, tol);

  // For now return true.  We should put a check in here.
  return true;
}

SharedJacobian& LinearOperator::getSharedJacobian()
{
  return *sharedJacobian;
}

Interface& LinearOperator::getUserInterface()
{
  return userInterface;
}

 
bool LinearOperator::computeRHS(const Epetra_Vector& x, Epetra_Vector& RHS)
{
  return userInterface.computeRHS(x, RHS);
}
  
bool LinearOperator::computeJacobian(const Epetra_Vector& x, Epetra_RowMatrix& Jac)
{
  bool status = false;
  
  if (jacType == "User Supplied") {
    status = userInterface.computeJacobian(x, Jac);
  }
  else if (jacType == "Finite Difference") {
    status = (dynamic_cast<FiniteDifference&>(jacobian)).computeJacobian(x, Jac);
  }
  else if (jacType == "Matrix-Free") {
    status = (dynamic_cast<MatrixFree&>(Jac)).computeJacobian(x, Jac);
  }

  return status;
}

bool LinearOperator::computePreconditioner(Epetra_RowMatrix& M)
{
  return userInterface.computePreconditioner(M);
}

bool LinearOperator::preconditionVector(Epetra_Vector& y)
{
  return userInterface.preconditionVector(y);
}

string LinearOperator::getJacobianType()
{
  Epetra_RowMatrix* test = 0;
  
  // Is it a "Matrix-Free" Jacobian?
  test = dynamic_cast<MatrixFree*>(&jacobian);
  if (test != 0) return "Matrix-Free"; 
  
  // Is it a "Finite Difference" Jacobian?
  test = dynamic_cast<FiniteDifference*>(&jacobian);
  if (test != 0) return "Finite Difference"; 
  
  // Otherwise it must be a "User Supplied" Epetra_RowMatrix.
  return "User Supplied";
}

