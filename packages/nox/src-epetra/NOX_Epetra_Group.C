// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "NOX_Epetra_Group.H"	// class definition
#include "NOX_Epetra_Interface.H"
#include "NOX_Epetra_SharedOperator.H"
#include "NOX_Parameter_List.H"
#include "NOX_Epetra_MatrixFree.H"
#include "NOX_Epetra_FiniteDifference.H"
#include "NOX_Epetra_Operator.H"

// External include files - linking to Aztec00 and Epetra in Trilinos
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Vector.h" 
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

using namespace NOX;
using namespace NOX::Epetra;

Group::Group(const Parameter::List& params, Interface& i, 
	     Epetra_Vector& x, Epetra_Operator& J):
  xVectorPtr(new Vector(x, DeepCopy)), // deep copy x     
  xVector(*xVectorPtr),    
  RHSVectorPtr(new Vector(x, ShapeCopy)), // new vector of same size
  RHSVector(*RHSVectorPtr), 
  gradVectorPtr(new Vector(x, ShapeCopy)), // new vector of same size
  gradVector(*gradVectorPtr), 
  NewtonVectorPtr(new Vector(x, ShapeCopy)), // new vector of same size
  NewtonVector(*NewtonVectorPtr), 
  tmpVectorPtr(0),
  sharedJacobianPtr(new SharedOperator(J)), // pass J to SharedJacobian
  sharedJacobian(*sharedJacobianPtr), // pass J to SharedJacobian
  sharedPreconditionerPtr(0),  // separate preconditioner is not used in this ctor
  sharedPreconditioner(*sharedJacobianPtr),  // point to the Jacobian
  userInterface(i)
{
  // Set all isValid flags to false
  resetIsValid();

  // Set the operators
  jacobianOperatorType = getOperatorType(J);
  preconditionerOperatorType = None; // This ctor is for one operator only!

  // Set the requested preconditioning.  Defaults to "None".
  preconditioner = params.getParameter("Preconditioning", "None");

  // Make sure the correct underlying objects were supplied for the requested
  // preconditioning options.
  checkOperatorConsistency();
}

Group::Group(const Parameter::List& params, Interface& i, 
	     Vector& x, Epetra_Operator& J):
  xVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(x.clone(DeepCopy))),
  xVector(*xVectorPtr),    
  RHSVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(x.clone(ShapeCopy))),
  RHSVector(*RHSVectorPtr), 
  gradVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(x.clone(ShapeCopy))),
  gradVector(*gradVectorPtr), 
  NewtonVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(x.clone(ShapeCopy))),
  NewtonVector(*NewtonVectorPtr), 
  tmpVectorPtr(0),
  sharedJacobianPtr(new SharedOperator(J)), // pass J to SharedJacobian
  sharedJacobian(*sharedJacobianPtr), // pass J to SharedJacobian
  sharedPreconditionerPtr(0),  // separate preconditioner is not used in this ctor
  sharedPreconditioner(*sharedJacobianPtr),  // point to the Jacobian
  userInterface(i)
{
  // Set all isValid flags to false
  resetIsValid();

  // Set the operators
  jacobianOperatorType = getOperatorType(J);
  preconditionerOperatorType = None; // This ctor is for one operator only!

  // Set the requested preconditioning.  Defaults to "None".
  preconditioner = params.getParameter("Preconditioning", "None");

  // Make sure the correct underlying objects were supplied for the requested
  // preconditioning options.
  checkOperatorConsistency();
}

Group::Group(const Parameter::List& params, Interface& i, 
	     Epetra_Vector& x, Epetra_Operator& J, Epetra_Operator& M):
  xVectorPtr(new Vector(x, DeepCopy)), // deep copy x     
  xVector(*xVectorPtr),    
  RHSVectorPtr(new Vector(x, ShapeCopy)), // new vector of same size
  RHSVector(*RHSVectorPtr), 
  gradVectorPtr(new Vector(x, ShapeCopy)), // new vector of same size
  gradVector(*gradVectorPtr), 
  NewtonVectorPtr(new Vector(x, ShapeCopy)), // new vector of same size
  NewtonVector(*NewtonVectorPtr), 
  tmpVectorPtr(0),
  sharedJacobianPtr(new SharedOperator(J)), // pass J to SharedOperator
  sharedJacobian(*sharedJacobianPtr), // create reference from pointer
  sharedPreconditionerPtr(new SharedOperator(M)), // pass M to SharedOperator
  sharedPreconditioner(*sharedPreconditionerPtr), // pass M to SharedOperator
  userInterface(i)
{
  // Set all isValid flags to false
  resetIsValid();

  // Set the operators
  jacobianOperatorType = getOperatorType(J);
  preconditionerOperatorType = getOperatorType(M); 

  // Set the requested preconditioning.  Defaults to "None".
  preconditioner = params.getParameter("Preconditioning", "None");

  // Make sure the correct underlying objects were supplied for the requested
  // preconditioning options.
  checkOperatorConsistency();
}

Group::Group(const Parameter::List& params, Interface& i, 
	     Vector& x, Epetra_Operator& J, Epetra_Operator& M):
  xVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(x.clone(DeepCopy))),
  xVector(*xVectorPtr),    
  RHSVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(x.clone(ShapeCopy))),
  RHSVector(*RHSVectorPtr), 
  gradVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(x.clone(ShapeCopy))),
  gradVector(*gradVectorPtr), 
  NewtonVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(x.clone(ShapeCopy))),
  NewtonVector(*NewtonVectorPtr), 
  tmpVectorPtr(0),
  sharedJacobianPtr(new SharedOperator(J)), // pass J to SharedOperator
  sharedJacobian(*sharedJacobianPtr), // create reference from pointer
  sharedPreconditionerPtr(new SharedOperator(M)), // pass M to SharedOperator
  sharedPreconditioner(*sharedPreconditionerPtr), // pass M to SharedOperator
  userInterface(i)
{
  // Set all isValid flags to false
  resetIsValid();

  // Set the operators
  jacobianOperatorType = getOperatorType(J);
  preconditionerOperatorType = getOperatorType(M); 

  // Set the requested preconditioning.  Defaults to "None".
  preconditioner = params.getParameter("Preconditioning", "None");

  // Make sure the correct underlying objects were supplied for the requested
  // preconditioning options.
  checkOperatorConsistency();
}

Group::Group(const Group& source, CopyType type) :
  xVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(source.xVector.clone(type))),
  xVector(*xVectorPtr),    
  RHSVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(source.RHSVector.clone(type))),
  RHSVector(*RHSVectorPtr), 
  gradVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(source.gradVector.clone(type))),
  gradVector(*gradVectorPtr), 
  NewtonVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(source.NewtonVector.clone(type))),
  NewtonVector(*NewtonVectorPtr), 
  tmpVectorPtr(0),
  sharedJacobianPtr(0),
  sharedJacobian(source.sharedJacobian),
  sharedPreconditionerPtr(0),
  sharedPreconditioner(source.sharedPreconditioner),
  userInterface(source.userInterface),
  jacobianOperatorType(source.jacobianOperatorType),
  preconditionerOperatorType(source.preconditionerOperatorType),
  preconditioner(source.preconditioner)
{
 
  switch (type) {
    
  case DeepCopy:
    
    isValidRHS = source.isValidRHS;
    isValidGrad = source.isValidGrad;
    isValidNewton = source.isValidNewton;
    isValidJacobian = source.isValidJacobian;
    isValidNormNewtonSolveResidual = source.isValidNormNewtonSolveResidual;
    normRHS = source.normRHS;
    normNewtonSolveResidual = source.normNewtonSolveResidual;
    
    // New copy takes ownership of the shared Jacobian
    if (isValidJacobian)
      sharedJacobian.getOperator(this);

    break;

  case ShapeCopy:
    resetIsValid();
    break;

  default:
    cerr << "ERROR: Invalid ConstructorType for group copy constructor." << endl;
    throw "NOX Error";
  }

}

Group::~Group() 
{
  delete tmpVectorPtr;
  delete sharedJacobianPtr;
  delete sharedPreconditionerPtr;
  delete NewtonVectorPtr;
  delete gradVectorPtr;
  delete RHSVectorPtr;
  delete xVectorPtr;
}

void Group::resetIsValid() //private
{
  isValidRHS = false;
  isValidJacobian = false;
  isValidGrad = false;
  isValidNewton = false;
  isValidNormNewtonSolveResidual = false;
}

void Group::setAztecOptions(const Parameter::List& p, AztecOO& aztec) const
{
  // Set the Aztec Solver
  string linearSolver = p.getParameter("Aztec Solver", "GMRES");
  if (linearSolver == "CG")
    aztec.SetAztecOption(AZ_solver, AZ_cg);
  else if (linearSolver == "GMRES")
    aztec.SetAztecOption(AZ_solver, AZ_gmres);
  else if (linearSolver == "CGS")
    aztec.SetAztecOption(AZ_solver, AZ_cgs);
  else if (linearSolver == "TFQMR")
    aztec.SetAztecOption(AZ_solver, AZ_tfqmr);
  else if (linearSolver == "BiCGStab")
    aztec.SetAztecOption(AZ_solver, AZ_bicgstab);
  else if (linearSolver == "LU")
    aztec.SetAztecOption(AZ_solver, AZ_lu);
  else {
    cout << "ERROR: NOX::Epetra::Group::setAztecOptions" << endl
	 << "\"Aztec Solver\" parameter \"" << linearSolver 
	 <<  "\" is invalid!" << endl;
    throw "NOX Error";
  }
 
  // Preconditioning Matrix Type 
  if (preconditioner == "None")
    aztec.SetAztecOption(AZ_precond, AZ_none);

  // Preconditioning where AztecOO inverts the Preconditioning Matrix
  if ((preconditioner == "AztecOO: Jacobian Matrix") || 
      (preconditioner == "AztecOO: User RowMatrix")) {
    
    string aztecPreconditioner = p.getParameter("Aztec Preconditioner", "ilu");

    if (aztecPreconditioner == "ilu") {
      aztec.SetAztecOption(AZ_precond, AZ_dom_decomp);
      aztec.SetAztecOption(AZ_overlap, p.getParameter("Overlap", 0));
      aztec.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
      aztec.SetAztecOption(AZ_graph_fill, p.getParameter("Graph Fill", 0));
    }
    else if (aztecPreconditioner == "ilut") { 
      aztec.SetAztecOption(AZ_precond, AZ_dom_decomp);
      aztec.SetAztecOption(AZ_overlap, p.getParameter("Overlap", 0));
      aztec.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
      aztec.SetAztecParam(AZ_drop, p.getParameter("Drop Tolerance", 0.0));
      aztec.SetAztecParam(AZ_ilut_fill, p.getParameter("Fill Factor", 1.0));
    }
    else if (aztecPreconditioner == "Jacobi") {
      aztec.SetAztecOption(AZ_precond, AZ_Jacobi);
      aztec.SetAztecOption(AZ_poly_ord, p.getParameter("Steps", 3));
    }
    else if (aztecPreconditioner == "Symmetric Gauss-Siedel") {
      aztec.SetAztecOption(AZ_precond, AZ_sym_GS);
      aztec.SetAztecOption(AZ_poly_ord, p.getParameter("Steps", 3));
    }
    else if (aztecPreconditioner == "Polynomial") {
      aztec.SetAztecOption(AZ_precond, AZ_Neumann);
      aztec.SetAztecOption(AZ_poly_ord, p.getParameter("Polynomial Order", 3));
    }
    else if (aztecPreconditioner == "Least-squares Polynomial") {
      aztec.SetAztecOption(AZ_precond, AZ_ls);
      aztec.SetAztecOption(AZ_poly_ord, p.getParameter("Polynomial Order", 3));
    }
    else {
      cout << "ERROR: NOX::Epetra::Group::setAztecOptions" << endl
	   << "\"Aztec Preconditioner\" parameter \"" << aztecPreconditioner
	   << "\" is invalid!" << endl;
      throw "NOX Error";
    }

  }

  // Turns on RCM reordering in conjunction with domain decomp preconditioning
  // default is "Disabled" = no reordering
  string RcmReordering = p.getParameter("RCM Reordering", "Disabled");
  if (RcmReordering == "Enabled") 
    aztec.SetAztecOption(AZ_reorder, 1);
  else if (RcmReordering == "Disabled")
    aztec.SetAztecOption(AZ_reorder, 0);
  else {
    cout << "ERROR: NOX::Epetra::Group::setAztecOptions" << endl
	 << "\"RCM Reordering\" parameter \"" << RcmReordering
	 << "\" is invalid!" << endl;
    throw "NOX Error";
  }
    
  // Gram-Schmidt orthogonalization procedure
  string orthog = p.getParameter("Orthogonalization", "Classical");
  if (orthog == "Classical") 
    aztec.SetAztecOption(AZ_orthog, AZ_classic);
  else if (RcmReordering == "Modified")
    aztec.SetAztecOption(AZ_orthog, AZ_modified);
  else {
    cout << "ERROR: NOX::Epetra::Group::setAztecOptions" << endl
	 << "\"Orthogonalization\" parameter \"" << orthog
	 << "\" is invalid!" << endl;
    throw "NOX Error";
  }

  // Size of the krylov space
  aztec.SetAztecOption(AZ_kspace, p.getParameter("Size of Krylov Subspace", 300));

  // Convergence criteria to use in the linear solver
  string convCriteria = p.getParameter("Convergence Criteria", "r0");
  if (convCriteria == "r0") 
    aztec.SetAztecOption(AZ_conv, AZ_r0);
  else if (convCriteria == "rhs")
    aztec.SetAztecOption(AZ_conv, AZ_rhs);
  else if (convCriteria == "Anorm")
    aztec.SetAztecOption(AZ_conv, AZ_Anorm);
  else if (convCriteria == "no scaling")
    aztec.SetAztecOption(AZ_conv, AZ_noscaled);
  else if (convCriteria == "sol")
    aztec.SetAztecOption(AZ_conv, AZ_sol);
  else {
    cout << "ERROR: NOX::Epetra::Group::setAztecOptions" << endl
	 << "\"Convergence Criteria\" parameter \"" << convCriteria
	 << "\" is invalid!" << endl;
    throw "NOX Error";
  }

  // Frequency of linear solve residual output
  aztec.SetAztecOption(AZ_output, 
		       p.getParameter("Output Frequency", AZ_last));
}

Abstract::Group* Group::clone(CopyType type) const 
{
  Group* newgrp = new Group(*this, type);
  return newgrp;
}

Abstract::Group& Group::operator=(const Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

Abstract::Group& Group::operator=(const Group& source)
{
  // Copy the xVector
  xVector = source.xVector;

  // Copy reference to sharedJacobian
  sharedJacobian = source.sharedJacobian;

  // Update the isValidVectors
  isValidRHS = source.isValidRHS;
  isValidGrad = source.isValidGrad;
  isValidNewton = source.isValidNewton;
  isValidJacobian = source.isValidJacobian;
  isValidNormNewtonSolveResidual = source.isValidNormNewtonSolveResidual;

  // Only copy vectors that are valid
  if (isValidRHS) {
    RHSVector = source.RHSVector;
    normRHS = source.normRHS;
  }

  if (isValidGrad)
    gradVector = source.gradVector;

  if (isValidNewton)
    NewtonVector = source.NewtonVector;

  if (isValidNormNewtonSolveResidual)
    normNewtonSolveResidual = source.normNewtonSolveResidual;

  // If valid, this takes ownership of the shared Jacobian
  if (isValidJacobian)
    sharedJacobian.getOperator(this);
    
  // Copy linear solver options
  jacobianOperatorType = source.jacobianOperatorType;
  preconditionerOperatorType = source.preconditionerOperatorType;
  preconditioner = source.preconditioner;

  return *this;
}

bool Group::setX(const Abstract::Vector& y) 
{
  return setX(dynamic_cast<const Vector&> (y));
}

bool Group::setX(const Vector& y) 
{
  resetIsValid();
  xVector = y;
  return true;
}

bool Group::computeX(const Abstract::Group& grp, 
		     const Abstract::Vector& d, 
		     double step) 
{
  // Cast to appropriate type, then call the "native" computeX
  const Group& epetragrp = dynamic_cast<const Group&> (grp);
  const Vector& epetrad = dynamic_cast<const Vector&> (d);
  return computeX(epetragrp, epetrad, step); 
}

bool Group::computeX(const Group& grp, const Vector& d, double step) 
{
  resetIsValid();
  xVector.update(1.0, grp.xVector, step, d);
  return true;
}

bool Group::computeF() 
{
  if (isF())
    return true;

  bool status = false;
  
  status = userInterface.computeF(xVector.getEpetraVector(), RHSVector.getEpetraVector());

  if (status == false) {
    cout << "ERROR: Epetra::Group::computeF() - fill failed!!!"
	 << endl;
    throw "NOX Error: Fill Failed";
  } 

  normRHS = RHSVector.norm();

  isValidRHS = true;

  return status;
}

bool Group::computeJacobian() 
{
  // Skip if the Jacobian is already valid
  if (isJacobian())
    return true;

  // Take ownership of the Jacobian and get a reference to the underlying operator
  Epetra_Operator& Jacobian = sharedJacobian.getOperator(this);

  // Fill the Jacobian 
  bool status = false;
  
  if ((jacobianOperatorType == EpetraOperator) || 
      (jacobianOperatorType == EpetraRowMatrix)) {
    status = userInterface.computeJacobian(xVector.getEpetraVector(), Jacobian);
  }
  else if (jacobianOperatorType == NoxOperator) {
    status = (dynamic_cast<Operator&>(Jacobian)).compute(xVector.getEpetraVector(), &Jacobian);
  }
  else if (jacobianOperatorType == NoxFiniteDifferenceRowMatrix) {
    status = (dynamic_cast<FiniteDifference&>(Jacobian)).computeJacobian(xVector.getEpetraVector(), Jacobian);
  }
  else if (jacobianOperatorType == NoxMatrixFreeOperator) {
    status = (dynamic_cast<MatrixFree&>(Jacobian)).computeJacobian(xVector.getEpetraVector(), Jacobian);
  }
  else {
    cout << "ERROR: NOX::Epetra::Group::computeJacobian() - jacobianOperatorType flag is "
	 << "not valid!" << endl;
    throw "NOX Error";
  }

  if (status == false) {
    cout << "ERROR: Epetra::Group::computeJacobian() - fill failed!!!"
	 << endl;
    throw "NOX Error: Fill Failed";
  } 

  // Update status of Jacobian wrt solution vector
  isValidJacobian = true;

  return status;
}

bool Group::computeGradient() 
{
  if (isGradient())
    return true;
  
  if (!isF()) {
    cerr << "ERROR: NOX::Epetra::Group::computeGradient() - RHS is out of date wrt X!" << endl;
    throw "NOX Error";
  }

  if (!isJacobian()) {
    cerr << "ERROR: NOX::Epetra::Group::computeGradient() - Jacobian is out of date wrt X!" << endl;
    throw "NOX Error";
  }
  
  // Get a reference to the Jacobian (it's validity was checked above)
  const Epetra_Operator& Jacobian = sharedJacobian.getOperator(this);

  // Compute grad = Jacobian^T * RHS.
  bool UseTranspose = true;
  bool NoTranspose = false;
  (const_cast<Epetra_Operator&>(Jacobian)).SetUseTranspose(UseTranspose);
  Jacobian.Apply(RHSVector.getEpetraVector(), gradVector.getEpetraVector());
  (const_cast<Epetra_Operator&>(Jacobian)).SetUseTranspose(NoTranspose);

  // Update state
  isValidGrad = true;

  // Return result
  return true;
}

bool Group::computeNewton(NOX::Parameter::List& p) 
{
  if (isNewton())
    return true;

  if (!isF()) {
    cerr << "ERROR: NOX::Epetra::Group::computeNewton() - invalid RHS" << endl;
    throw "NOX Error";
  }

  if (!isJacobian()) {
    cerr << "ERROR: NOX::Epetra::Group::computeNewton() - invalid Jacobian" << endl;
    throw "NOX Error";
  }
  
  // Create Epetra problem for the linear solve
  applyJacobianInverse(p, RHSVector, NewtonVector);

  // Scale soln by -1
  NewtonVector.scale(-1.0);

  // Update state
  isValidNewton = true;

  // Compute the 2-norm of the Newton solve residual ||Js+f||
  computeNormNewtonSolveResidual();

  // Return solution
  return true;
}

bool Group::applyJacobian(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& epetrainput = dynamic_cast<const Vector&> (input);
  Vector& epetraresult = dynamic_cast<Vector&> (result);
  return applyJacobian(epetrainput, epetraresult);
}

bool Group::applyJacobian(const Vector& input, Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian()) 
    return false;

  // Get a reference to the Jacobian (it's validity was checked above)
  const Epetra_Operator& Jacobian = sharedJacobian.getOperator();

  // Apply the Jacobian
  bool NoTranspose = false;
  (const_cast <Epetra_Operator&>(Jacobian)).SetUseTranspose(NoTranspose);
  Jacobian.Apply(input.getEpetraVector(), result.getEpetraVector());

  return true;
}

bool Group::applyJacobianInverse (Parameter::List &p, const Abstract::Vector &input, Abstract::Vector &result) const
{
  const Vector& epetraInput = dynamic_cast<const Vector&>(input);
  Vector& epetraResult = dynamic_cast<Vector&>(result);
  return applyJacobianInverse(p, epetraInput, epetraResult);
}

bool Group::applyJacobianInverse (Parameter::List &p, const Vector &input, Vector &result) const
{
  // Get the Jacobian 
  /* Have to get non-const version which requires reasserting
   * ownership. 
   */
  Epetra_Operator& Jacobian = sharedJacobian.getOperator(this);

  // Epetra_LinearProblem does not allow for const vectors
  Vector& nonConstInput = const_cast<Vector&>(input);

  // Create Epetra problem for the linear solve
  Epetra_LinearProblem Problem(&Jacobian, 
  			       &(result.getEpetraVector()), 
			       &(nonConstInput.getEpetraVector()));

  // Set the default Problem parameters to "hard" (this sets Aztec defaults
  // during the AztecOO instantiation)
  Problem.SetPDL(hard);

  // Scale the problem if requested
  string scalingOption = p.getParameter("Scaling", "None");
  if (scalingOption == "Row Sum") {
    if (tmpVectorPtr == 0) 
      tmpVectorPtr = new Epetra_Vector(xVector.getEpetraVector());

    // make sure the Jacobian is an Epetra_RowMatrix, otherwise we can't scale
    Epetra_RowMatrix* test = 0;
    test = dynamic_cast<Epetra_RowMatrix*>(&Jacobian);
    if (test == 0) {
      cout << "ERROR: NOX::Epetra::Group::applyJacobianInverse() - "
	   << "For \"Row Sum\" scaling, the Jacobian must be an "
	   << "Epetra_RowMatrix derived object!" << endl;
      throw "NOX Error";
    }

    test->InvRowSums(*tmpVectorPtr);
    Problem.LeftScale(*tmpVectorPtr);
  }
  else if (scalingOption == "None") {
    // Do nothing
  }
  else {
    // Throw an error, the requested scaling option is not vaild
    cout << "ERROR: NOX::Epetra::Group::applyJacobianInverse() - "
	 << " The parameter chosen for \"Scaling\" is not valid!"
	 << endl;
    throw "NOX Error";
  }
  
  // Create aztec problem
  AztecOO aztec(Problem);    

  // Set specific Aztec parameters requested by NOX
  setAztecOptions(p, aztec);

  // Compute and set the Preconditioner in AztecOO if needed
  if (preconditioner != "None") {
    computePreconditioner(aztec);
  }

  int maxit = p.getParameter("Max Iterations", 400);
  double tol = p.getParameter("Tolerance", 1.0e-6);

  // Solve Aztec problem
  int aztecStatus = aztec.Iterate(maxit, tol);

  // Unscale the linear problem
  if (scalingOption == "Row Sum") {
    tmpVectorPtr->Reciprocal(*tmpVectorPtr);
    Problem.LeftScale(*tmpVectorPtr);
  }

  if (aztecStatus != 0) 
    return false;
  
  // Set the output parameters
  NOX::Parameter::List& outputList = p.sublist("Output");
  int linearIters = outputList.getParameter("Number of Linear Iterations", 0);
  outputList.setParameter("Number of Linear Iterations", 
		 (linearIters + aztec.NumIters()));
  outputList.setParameter("True Unscaled Residual", aztec.TrueResidual());

  return true;
}

bool Group::applyJacobianDiagonalInverse(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& epetrainput = dynamic_cast<const Vector&> (input);
  Vector& epetraresult = dynamic_cast<Vector&> (result);
  return applyJacobianDiagonalInverse(epetrainput, epetraresult);
}

bool Group::applyJacobianDiagonalInverse(const Vector& input, Vector& result) const
{
  if (!isJacobian()) 
    return false;

  /* To extract the diagonal inverse, we must have a real matrix 
   * (we can NOT be using matrix-free operators).  Thus it must 
   * be an Epetra_RowMatrix.  Check for this and if not, throw an 
   * error.
   */
  // Get a reference to the Jacobian operator 
  const Epetra_Operator& tmpJacobian = sharedJacobian.getOperator();
  // Try and cast it to an Epetra_RowMatrix
  const Epetra_RowMatrix* testRowMatrix = dynamic_cast<const Epetra_RowMatrix*>(&tmpJacobian);
  if (testRowMatrix == NULL) {
    cout << "ERROR: NOX::Epetra::Group::applyJacobianDiagonalInverse) - "
	 << "the Jacobian operator is NOT a matrix!" << endl;
    throw "NOX Error";
  }

  // Convert the Epetra_RowMatrix this to a reference
  const Epetra_RowMatrix& Jacobian = *testRowMatrix;

  // Get epetra reference to the result vector
  Epetra_Vector& r = result.getEpetraVector();

  // Allocate the extra tmpVectorPtr if necessary
  if (tmpVectorPtr == NULL)
    tmpVectorPtr = new Epetra_Vector(r.Map());

  // Get the reference to the temporary vector
  Epetra_Vector& tmpVector = *tmpVectorPtr;

  // Put a copy of the diagonal of the Jacobian into tmpVector
  int retcode = Jacobian.ExtractDiagonalCopy(tmpVector);
  
  // Take element-wise absolute value of diagonal vector
  retcode = tmpVector.Abs(tmpVector);
  
  // Check minimum absolute value of diagonal vector
  double minAbsValue = 0;
  retcode = tmpVector.MinValue(&minAbsValue);

  if(minAbsValue <= 1.e-6) // This minimum threshold can be adjusted
  {
    cout << "Poor scaling on Jacobian diagonal (min abs value: " <<
             minAbsValue << " ) --> NO nonlinear Preconditioning !!" << endl;
    return false; 
  }
  
  // Check if ExtractDiagonalCopy is supported
  if (retcode != 0)
    return false;

  // Calculate r = input ./ tmpVector (./ is element-by-element divide)
  retcode = r.ReciprocalMultiply(1.0, tmpVector, input.getEpetraVector(), 0.0);

  // Check if this worked
  if (retcode != 0)
    return false;

  return true;
}

bool Group::applyJacobianTranspose(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& epetrainput = dynamic_cast<const Vector&> (input);
  Vector& epetraresult = dynamic_cast<Vector&> (result);
  return applyJacobianTranspose(epetrainput, epetraresult);
}

bool Group::applyJacobianTranspose(const Vector& input, Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian()) 
    return false;

  // Get a reference to the Jacobian (it's validity was check above)
  const Epetra_Operator& Jacobian = sharedJacobian.getOperator();

  // Apply the Jacobian
  bool UseTranspose = true;
  bool NoTranspose = false;
  (const_cast<Epetra_Operator&>(Jacobian)).SetUseTranspose(UseTranspose);
  Jacobian.Apply(input.getEpetraVector(), result.getEpetraVector());
  (const_cast<Epetra_Operator&>(Jacobian)).SetUseTranspose(NoTranspose);

  return true;
}



bool Group::isF() const 
{   
  return isValidRHS;
}

bool Group::isJacobian() const 
{  
  return ((sharedJacobian.isOwner(this)) && (isValidJacobian));
}

bool Group::isGradient() const 
{   
  return isValidGrad;
}

bool Group::isNewton() const 
{   
  return isValidNewton;
}

bool Group::isNormNewtonSolveResidual() const 
{   
  return isValidNormNewtonSolveResidual;
}

const Abstract::Vector& Group::getX() const 
{
  return xVector;
}

const Abstract::Vector& Group::getF() const 
{  
  if (!isF()) {
    cerr << "ERROR: NOX::Epetra::Group::getF() - invalid RHS" << endl;
    throw "NOX Error";
  }
    
  return RHSVector;
}

double Group::getNormF() const
{
  if (!isF()) {
    cerr << "ERROR: NOX::Epetra::Group::getNormF() - invalid RHS" << endl;
    throw "NOX Error";
  }
    
  return normRHS;
}

const Abstract::Vector& Group::getGradient() const 
{ 
  if (!isGradient()) {
    cerr << "ERROR: NOX::Epetra::Group::getGradient() - invalid gradient" << endl;
    throw "NOX Error";
  }
    
  return gradVector;
}

const Abstract::Vector& Group::getNewton() const 
{
  if (!isNewton()) {
    cerr << "ERROR: NOX::Epetra::Group::getNewton() - invalid Newton vector" << endl;
    throw "NOX Error";
  }
    
  return NewtonVector;
}

double NOX::Epetra::Group::getNormNewtonSolveResidual() const
{
  // Make sure value is not already calculated
  if (isValidNormNewtonSolveResidual) 
    return normNewtonSolveResidual;

  // Otherwise throw an error since a Newton direction has not been calculated
  // wrt this solution group
  /*
  cout << "ERROR: NOX::Epetra::Group::getNormNewtonSolveResidual() - Group has "
       << "not performed a Newton solve corresponding to this solution vector!"
       << endl;
  throw "NOX Error";
  */
  return -1.0;
}  

SharedOperator& Group::getSharedJacobian()
{
  return sharedJacobian;
}

SharedOperator& Group::getSharedPreconditioner()
{
  return sharedPreconditioner;
}

Interface& Group::getUserInterface()
{
  return userInterface;
}

bool Group::checkOperatorConsistency() 
{

  // Checks that the requested preconditioning option (preconditioner) has 
  // the correct objects avalable 
  Epetra_RowMatrix* test = 0;

  if (preconditioner == "AztecOO: Jacobian Matrix") {

    // Get a reference to the Jacobian 
    Epetra_Operator& Jacobian = sharedJacobian.getOperator(this);

    // The Jacobian Operator MUST be an Eptra_RowMatrix
    test = dynamic_cast<Epetra_RowMatrix*>(&Jacobian);

    if (test == 0) {
      cout << "ERROR: NOX::Epetra::Group::setLinearSolver() - The flag "
	   << "\"Preconditioning\" with value \"" << preconditioner
	   << "\" requires an Epetra_RowMatrix object for the Jacobian!" 
	   << endl;
      cout << "Matrx-Free Jacobians can not be used for preconditioners!"
	   << endl;
      throw "NOX Error";
    }
  }
  else if (preconditioner == "AztecOO: User RowMatrix") {
    
    // Make sure a separate operator was supplied by the user
    if (preconditionerOperatorType == None) {
      cout << "ERROR: NOX::Epetra::Group::setLinearSolver() - The flag "
	   << "\"Preconditioning\" with value \"" << preconditioner
	   << "\" requires a NOX::Epetra::Group constructor with "
	   << "separate objects for the Jacobian and preconditioner!" << endl
	   << endl;
      throw "NOX Error";
    }

    // The Preconditioning Operator MUST BE an Epetra_RowMatrix
    Epetra_Operator& Prec = sharedPreconditionerPtr->getOperator(this);
    test = dynamic_cast<Epetra_RowMatrix*>(&Prec);
    if (test == 0) {
      cout << "ERROR: NOX::Epetra::Group::setLinearSolver() - The flag "
	   << "\"Preconditioning\" with value \"" << preconditioner
	   << "\" requires an Epetra_RowMatrix object for the "
	   << "preconditioner!" << endl;
      throw "NOX Error";
    }
  }
  else if (preconditioner == "User Supplied Preconditioner") {

    // Make sure the operator was supplied by the user
    if (preconditionerOperatorType == None) {
      cout << "ERROR: NOX::Epetra::Group::setLinearSolver() - The flag "
	   << "\"Preconditioning\" with value \"" << preconditioner
	   << "\" requires a NOX::Epetra::Group constructor with "
	   << "separate objects for the Jacobian and preconditioner!" << endl; 
      throw "NOX Error";
    }
    
  }
  else if (preconditioner == "None") {
    // Do nothing. This is a valid argument and won't cause a throw but we still
    // want to catch mis-spelled "Preconditioner" flags (see next "else").
  } 
  else {
    // An invalid choice was found 
    cout << "ERROR: NOX::Epetra::Group::setLinearSolver() - The flag "
	 << "\"" << preconditioner << "\" is not a valid choice for the " 
	 << "parameter \"Preconditioning\"!" << endl;
      throw "NOX Error";
  }
  
  // All checks passed without a throw!
  return true;
}

bool Group::computePreconditioner(AztecOO& aztec) const
{
  
  if (preconditioner == "AztecOO: Jacobian Matrix") {
    // Do nothing, the Jacobian has already been evaluated at the
    // current solution.
  }
  else if (preconditioner == "AztecOO: User RowMatrix") {
    
    Epetra_RowMatrix& precMatrix = 
      dynamic_cast<Epetra_RowMatrix&>(sharedPreconditioner.getOperator(this));
    
    if (preconditionerOperatorType == NoxFiniteDifferenceRowMatrix) {
      (dynamic_cast<FiniteDifference&>(precMatrix))
	.computePreconditioner(xVector.getEpetraVector(), precMatrix);
    }
    else 
      userInterface.computePrecMatrix(xVector.getEpetraVector(), precMatrix);
    
    aztec.SetPrecMatrix(&precMatrix);
  }
  else if (preconditioner == "User Supplied Preconditioner") {
    
    Epetra_Operator& precOperator = sharedPreconditioner.getOperator(this);

    // For user supplied Precs the user may supply them in one of two ways:
    if (preconditionerOperatorType == NoxOperator) {
      // (1) Use a NOX::Epetra::Operator derived object
      const Epetra_Operator* jacobianPtr = 0; 
      if (isJacobian())
	jacobianPtr = &sharedJacobian.getOperator();

      (dynamic_cast<NOX::Epetra::Operator&>(precOperator))
	.compute(xVector.getEpetraVector(), jacobianPtr);
    }
    else {
      // (2) Supply an Epetra_Operator derived class and implement 
      // the NOX::Epetra::Interface::computePreconditioner() method
      // to evaluate it at the current solution.
      userInterface.computePreconditioner(xVector.getEpetraVector(),
					  precOperator);
    }
    aztec.SetPrecOperator(&precOperator); 
  }

  return true;
}

bool Group::computeNormNewtonSolveResidual ()
{
  // Make sure value is not already calculated
  if (isValidNormNewtonSolveResidual) 
    return true;

  // Make sure NewtonVector and RHSVector are valid
  // We could return false, but for now we will throw errors
  if (!isValidRHS) {
    cerr << "ERROR: NOX::Epetra::Group::computeNormNewtonSolveResidual() - invalid RHS" 
	 << endl;
    throw "NOX Error";
  }
  if (!isValidNewton) {
    cerr << "ERROR: NOX::Epetra::Group::computeNormNewtonSolveResidual() - invalid "
	 << "Newton direction" << endl;
    throw "NOX Error";
  }
  
  // Allocate the tmpVectorPtr if not already done (deleted in ~Group)   
  if (tmpVectorPtr == 0) {
    tmpVectorPtr = new Epetra_Vector(RHSVector.getEpetraVector());
  }
  Vector tmpNoxVector(*tmpVectorPtr, ShapeCopy); 

  this->applyJacobian(NewtonVector, tmpNoxVector);    
  tmpNoxVector.update(1.0, RHSVector, 1.0);
  normNewtonSolveResidual = tmpNoxVector.norm();
  
  isValidNormNewtonSolveResidual = true;
  
  return true;
}

Group::OperatorType Group::getOperatorType(const Epetra_Operator& Op)
{
  // NOTE: The order in which the following tests occur is important!
  const Epetra_Operator* testOperator = 0;
  
  // Is it a NOX::Epetra::MatrixFree Operator?
  testOperator = dynamic_cast<const MatrixFree*>(&Op);
  if (testOperator != 0) 
    return NoxMatrixFreeOperator; 
  
  // Is it a NOX::Epetra::Finite Difference Operator?
  testOperator = dynamic_cast<const FiniteDifference*>(&Op);
  if (testOperator != 0) 
    return NoxFiniteDifferenceRowMatrix; 
  
  // Is it a NOX::Epetra::Operator ?
  testOperator = dynamic_cast<const Operator*>(&Op);
  if (testOperator != 0) 
    return NoxOperator; 

  // Is it an Epetra_RowMatrix ?
  testOperator = dynamic_cast<const Epetra_RowMatrix*>(&Op);
  if (testOperator != 0) 
    return EpetraRowMatrix; 

  // Otherwise it must be an Epetra_Operator!
  return EpetraOperator;
}
