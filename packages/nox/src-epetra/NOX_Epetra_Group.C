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
#include "NOX_Epetra_FiniteDifferenceColoring.H"
#include "NOX_Epetra_Operator.H"
#include "NOX_Utils.H"

// External include files - linking to Aztec00 and Epetra in Trilinos
#include "AztecOO.h"
#include "AztecOO_Operator.h"
#include "AztecOO_StatusTest.h"
#include "AztecOO_StatusTestCombo.h"
#include "AztecOO_StatusTestMaxIters.h"
#include "AztecOO_StatusTestResNorm.h"
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Vector.h" 
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"

using namespace NOX;
using namespace NOX::Epetra;

Group::Group(Parameter::List& printParams, 
	     Parameter::List& params, Interface& i, 
	     Epetra_Vector& x, Epetra_Operator& J):
  utils(printParams),
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
  userInterface(i),
  aztecSolver(0),
  ifpackGraph(0),
  ifpackPreconditioner(0)
{
  // Set all isValid flags to false
  resetIsValid();

  // Set the operators
  jacobianOperatorType = getOperatorType(J);
  preconditionerOperatorType = None; // This ctor is for one operator only!

  // Set the requested preconditioning.  Defaults to "None".
  preconditioner = params.getParameter("Preconditioning", "None");
  isValidPreconditioner = false;

  // Make sure the correct underlying objects were supplied for the requested
  // preconditioning options.
  checkOperatorConsistency();
}

Group::Group(Parameter::List& printParams, 
	     Parameter::List& params, Interface& i, 
	     Vector& x, Epetra_Operator& J):
  utils(printParams),
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
  userInterface(i),
  aztecSolver(0),
  ifpackGraph(0),
  ifpackPreconditioner(0)
{
  // Set all isValid flags to false
  resetIsValid();

  // Set the operators
  jacobianOperatorType = getOperatorType(J);
  preconditionerOperatorType = None; // This ctor is for one operator only!

  // Set the requested preconditioning.  Defaults to "None".
  preconditioner = params.getParameter("Preconditioning", "None");
  isValidPreconditioner = false;

  // Make sure the correct underlying objects were supplied for the requested
  // preconditioning options.
  checkOperatorConsistency();
}

Group::Group(Parameter::List& printParams, 
	     Parameter::List& params, Interface& i, 
	     Epetra_Vector& x, Epetra_Operator& J, Epetra_Operator& M):
  utils(printParams),
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
  userInterface(i),
  aztecSolver(0),
  ifpackGraph(0),
  ifpackPreconditioner(0)
{
  // Set all isValid flags to false
  resetIsValid();

  // Set the operators
  jacobianOperatorType = getOperatorType(J);
  preconditionerOperatorType = getOperatorType(M); 

  // Set the requested preconditioning.  Defaults to "None".
  preconditioner = params.getParameter("Preconditioning", "None");
  isValidPreconditioner = false;

  // Make sure the correct underlying objects were supplied for the requested
  // preconditioning options.
  checkOperatorConsistency();
}

Group::Group(Parameter::List& printParams, 
	     Parameter::List& params, Interface& i, 
	     Vector& x, Epetra_Operator& J, Epetra_Operator& M):
  utils(printParams),
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
  userInterface(i),
  aztecSolver(0),
  ifpackGraph(0),
  ifpackPreconditioner(0)
{
  // Set all isValid flags to false
  resetIsValid();

  // Set the operators
  jacobianOperatorType = getOperatorType(J);
  preconditionerOperatorType = getOperatorType(M); 

  // Set the requested preconditioning.  Defaults to "None".
  preconditioner = params.getParameter("Preconditioning", "None");
  isValidPreconditioner = false;

  // Make sure the correct underlying objects were supplied for the requested
  // preconditioning options.
  checkOperatorConsistency();
}

Group::Group(const Group& source, CopyType type) :
  utils(source.utils),
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
  preconditioner(source.preconditioner),
  aztecSolver(0),
  ifpackGraph(0),
  ifpackPreconditioner(0)
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

    // New copy does NOT take ownership of the shared preconditioner nor
    // any preconditioning objects
    isValidPreconditioner = false;

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
  destroyPreconditioner();
  delete aztecSolver;
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
  // Never reset isValidPreconditioner unless you first 
  // call destroyPreconditioner().  I have commented this out since
  // resetIsValid() implarts no knowledge that the preconditioner is
  // also destroyed.
  //destroyPrecondtioner();
  //isValidPreconditioner = false;
  return;
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

  // Set the ill-conditioning threshold for the upper hessenberg matrix
  if (p.isParameter("Ill-Conditioning Threshold")) {
    aztec.SetAztecParam(AZ_ill_cond_thresh, 
			p.getParameter("Ill-Conditioning Threshold", 1.0e+11));
  }

  // Frequency of linear solve residual output
  aztec.SetAztecOption(AZ_output, 
		       p.getParameter("Output Frequency", AZ_last));

  // Print a summary of the aztec options if "Details" is enabled
  if (utils.isPrintType(Utils::LinearSolverDetails)) {
    aztec.CheckInput();
  }
  return;
}

Abstract::Group* Group::clone(CopyType type) const 
{
  NOX::Abstract::Group* newgrp = new NOX::Epetra::Group(*this, type);
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
    
  // Never take ownership of the shared preconditioner
  isValidPreconditioner = false;
    
  // Copy linear solver options
  jacobianOperatorType = source.jacobianOperatorType;
  preconditionerOperatorType = source.preconditionerOperatorType;
  preconditioner = source.preconditioner;

  return *this;
}

void Group::setX(const Abstract::Vector& y) 
{
  setX(dynamic_cast<const Vector&> (y));
  return;
}

void Group::setX(const Vector& y) 
{
  resetIsValid();
  destroyAztecSolver();
  xVector = y;
  return;
}

void Group::computeX(const Abstract::Group& grp, 
		     const Abstract::Vector& d, 
		     double step) 
{
  // Cast to appropriate type, then call the "native" computeX
  const Group& epetragrp = dynamic_cast<const Group&> (grp);
  const Vector& epetrad = dynamic_cast<const Vector&> (d);
  computeX(epetragrp, epetrad, step); 
  return;
}

void Group::computeX(const Group& grp, const Vector& d, double step) 
{
  resetIsValid();
  destroyAztecSolver();
  xVector.update(1.0, grp.xVector, step, d);
  return;
}

Abstract::Group::ReturnType Group::computeF() 
{
  if (isF())
    return Abstract::Group::Ok;

  bool status = false;
  
  status = userInterface.computeF(xVector.getEpetraVector(), RHSVector.getEpetraVector());

  if (status == false) {
    cout << "ERROR: Epetra::Group::computeF() - fill failed!!!"
	 << endl;
    throw "NOX Error: Fill Failed";
  } 

  normRHS = RHSVector.norm();

  isValidRHS = true;

  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType Group::computeJacobian() 
{
  // Skip if the Jacobian is already valid
  if (isJacobian())
    return Abstract::Group::Ok;

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

  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType Group::computeGradient() 
{
  if (isGradient())
    return Abstract::Group::Ok;
  
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
  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType Group::computeNewton(NOX::Parameter::List& p) 
{
  if (isNewton())
    return Abstract::Group::Ok;

  if (!isF()) {
    cerr << "ERROR: NOX::Epetra::Group::computeNewton() - invalid RHS" << endl;
    throw "NOX Error";
  }

  if (!isJacobian()) {
    cerr << "ERROR: NOX::Epetra::Group::computeNewton() - invalid Jacobian" << endl;
    throw "NOX Error";
  }

  Abstract::Group::ReturnType status;
  
  // Create Epetra problem for the linear solve
  status = applyJacobianInverse(p, RHSVector, NewtonVector);

  // Scale soln by -1
  NewtonVector.scale(-1.0);

  // Update state EVEN IF LINEAR SOLVE FAILED
  // We still may want to use the vector even it it just missed it's 
  isValidNewton = true;

  // Compute the 2-norm of the linear solve residual ||Js+f||
  computeNormNewtonSolveResidual();

  // return status of the linear solver
  if( status != Abstract::Group::Ok ) 
    return status;

  // Return solution
  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType Group::applyJacobian(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& epetrainput = dynamic_cast<const Vector&> (input);
  Vector& epetraresult = dynamic_cast<Vector&> (result);
  return applyJacobian(epetrainput, epetraresult);
}

Abstract::Group::ReturnType Group::applyJacobian(const Vector& input, Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian()) 
    return Abstract::Group::BadDependency;

  // Get a reference to the Jacobian (it's validity was checked above)
  const Epetra_Operator& Jacobian = sharedJacobian.getOperator();

  // Apply the Jacobian
  bool NoTranspose = false;
  (const_cast <Epetra_Operator&>(Jacobian)).SetUseTranspose(NoTranspose);
  // The next call could be made to return a Abstract::Group::ReturnType
  Jacobian.Apply(input.getEpetraVector(), result.getEpetraVector());

  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType Group::applyJacobianInverse (Parameter::List &p, const Abstract::Vector &input, Abstract::Vector &result) const
{
  const Vector& epetraInput = dynamic_cast<const Vector&>(input);
  Vector& epetraResult = dynamic_cast<Vector&>(result);
  return applyJacobianInverse(p, epetraInput, epetraResult);
}

Abstract::Group::ReturnType Group::applyJacobianInverse (Parameter::List &p, const Vector &input, Vector &result) const
{
  // Get the non-const versions of the Jacobian and input vector
  // Epetra_LinearProblem requires non-const versions so we can perform
  // scaling of the linear problem.
  Epetra_Operator& Jacobian = sharedJacobian.getOperator(this);
  Vector& nonConstInput = const_cast<Vector&>(input);
  
  // Zero out the delta X of the linear problem
  result.init(0.0);

  // Create Epetra linear problem object for the linear solve
  Epetra_LinearProblem Problem(&Jacobian, 
  			       &(result.getEpetraVector()), 
			       &(nonConstInput.getEpetraVector()));

  // ************* Begin linear system scaling *******************

  // Get the norm of the RHS BEFORE scaling the linear system.  This is 
  // required for adjustable forcing terms when used in combination with
  // row sum scaling. 
  double normF = input.norm();

  string scalingOption = p.getParameter("Scaling", "None");
  if (scalingOption == "Row Sum") {
    if (tmpVectorPtr == 0) 
      tmpVectorPtr = new Epetra_Vector(xVector.getEpetraVector());

    // Make sure the Jacobian is an Epetra_RowMatrix, otherwise we can't 
    // perform a row sum scale!
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
    
    // Get the row sum vector back (right now it is the inverse of the row sum)
    // We need this to (1) unscale the system and (2) for special convergence 
    // test if row sum scaling and adjustable forcing terms are used in 
    // combination.
    tmpVectorPtr->Reciprocal(*tmpVectorPtr);

    if (utils.isPrintProcessAndType(Utils::Details))
      cout << endl << "       Linear Problem Scaling: Row Sum" << endl;

  }
  else if (scalingOption == "None") {
    if (utils.isPrintProcessAndType(Utils::Details))
      cout << endl << "       Linear Problem Scaling: None" << endl;
  }
  else {
    // Throw an error, the requested scaling option is not vaild
    cout << "ERROR: NOX::Epetra::Group::applyJacobianInverse() - "
	 << " The parameter chosen for \"Scaling\" is not valid!"
	 << endl;
    throw "NOX Error";
  } 
  // ************* End linear system scaling *******************
  
  // Set the default Problem parameters to "hard" (this sets Aztec defaults
  // during the AztecOO instantiation)
  Problem.SetPDL(hard);

  // Create the solver.  To be safe, first remove the old solver/preconditionr 
  // that was created by applyRightPreconditioner().
  destroyAztecSolver();
  aztecSolver = new AztecOO(Problem);
  
  // Set specific Aztec parameters requested by NOX
  setAztecOptions(p, *aztecSolver);
  
  // Compute and set the Preconditioner in AztecOO if needed
  createPreconditioner(ImplicitConstruction, p);

  // Get linear solver convergence parameters
  int maxit = p.getParameter("Max Iterations", 400);
  double tol = p.getParameter("Tolerance", 1.0e-6);
  
  // If we are using an adjustable forcing term combined with Row Sum 
  // scaling, the convergence test in the linear solver must be modified.
  // Don't use the default ones in Aztec!
  AztecOO_StatusTestMaxIters* maxIters = 0;
  AztecOO_StatusTestResNorm* forcingResNorm = 0;
  AztecOO_StatusTestCombo* convTest = 0;
 
  // If an adjustable forcing term is being used, we need to explicitly 
  // enforce the use of a particular type of convergence test.  Additionally,
  // if scaling of the linear system is used, the status test must account
  // for this!
  bool usingAdjustableForcingTerm =  
    p.getParameter("Using Adjustable Forcing Term", false);
  
  if ((usingAdjustableForcingTerm) && (scalingOption == "Row Sum")) {

    maxIters = new AztecOO_StatusTestMaxIters(maxit);

    forcingResNorm = new AztecOO_StatusTestResNorm(Jacobian, 
		     result.getEpetraVector(), input.getEpetraVector(), tol);

    forcingResNorm->DefineResForm(AztecOO_StatusTestResNorm::Explicit,
    			  AztecOO_StatusTestResNorm::TwoNorm,
    			  tmpVectorPtr);

    forcingResNorm->DefineScaleForm(AztecOO_StatusTestResNorm::UserProvided,
				    AztecOO_StatusTestResNorm::TwoNorm, 0,
				    normF);

    convTest = new AztecOO_StatusTestCombo(AztecOO_StatusTestCombo::OR, 
					   *maxIters, *forcingResNorm);

    aztecSolver->SetStatusTest(convTest);
  }
  else if (usingAdjustableForcingTerm) {
    
    maxIters = new AztecOO_StatusTestMaxIters(maxit);

    forcingResNorm = new AztecOO_StatusTestResNorm(Jacobian, 
		     result.getEpetraVector(), input.getEpetraVector(), tol);

    forcingResNorm->DefineResForm(AztecOO_StatusTestResNorm::Implicit,
    			  AztecOO_StatusTestResNorm::TwoNorm);

    forcingResNorm->DefineScaleForm(AztecOO_StatusTestResNorm::NormOfRHS,
				    AztecOO_StatusTestResNorm::TwoNorm);

    convTest = new AztecOO_StatusTestCombo(AztecOO_StatusTestCombo::OR, 
					   *maxIters, *forcingResNorm);

    aztecSolver->SetStatusTest(convTest);

  }

  // Solve Aztec problem
  int aztecStatus = -1;
  if (p.getParameter("Use Adaptive Linear Solve", false)) {
    aztecSolver->SetUseAdaptiveDefaultsTrue();
    aztecStatus = aztecSolver->AdaptiveIterate(maxit, 
		   p.getParameter("Max Adaptive Solve Attempts", 5), tol);
  }
  else
    aztecStatus = aztecSolver->Iterate(maxit, tol);

  // Unscale the linear problem
  if (scalingOption == "Row Sum") {
    Problem.LeftScale(*tmpVectorPtr);
  }

  // Set the output parameters in the "Output" sublist
  NOX::Parameter::List& outputList = p.sublist("Output");
  int prevLinIters = outputList.getParameter("Total Number of Linear Iterations", 0);
  int curLinIters = 0;
  double achievedTol = -1.0;
  if (usingAdjustableForcingTerm) {
    curLinIters = maxIters->GetNumIters();
    achievedTol = forcingResNorm->GetTestValue();
  }
  else {
    curLinIters = aztecSolver->NumIters();
    achievedTol = aztecSolver->ScaledResidual();
  }
  outputList.setParameter("Number of Linear Iterations", curLinIters);
  outputList.setParameter("Total Number of Linear Iterations", (prevLinIters + curLinIters));
  outputList.setParameter("Achieved Tolerance", achievedTol);

  // Delete the special convergence tests
  delete convTest;
  delete maxIters;
  delete forcingResNorm;

  // Delete the solver and reset to NULL.
  destroyAztecSolver();

  if (aztecStatus != 0) 
    return Abstract::Group::NotConverged;
  
  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType Group::applyJacobianTranspose(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& epetrainput = dynamic_cast<const Vector&> (input);
  Vector& epetraresult = dynamic_cast<Vector&> (result);
  return applyJacobianTranspose(epetrainput, epetraresult);
}

Abstract::Group::ReturnType Group::applyJacobianTranspose(const Vector& input, Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian()) 
    return Abstract::Group::BadDependency;

  // Get a reference to the Jacobian (it's validity was check above)
  const Epetra_Operator& Jacobian = sharedJacobian.getOperator();

  // Apply the Jacobian
  bool UseTranspose = true;
  bool NoTranspose = false;
  (const_cast<Epetra_Operator&>(Jacobian)).SetUseTranspose(UseTranspose);
  Jacobian.Apply(input.getEpetraVector(), result.getEpetraVector());
  (const_cast<Epetra_Operator&>(Jacobian)).SetUseTranspose(NoTranspose);

  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType Group::applyRightPreconditioning(Parameter::List& params,
				      const Abstract::Vector& input, 
				      Abstract::Vector& result) const
{
  const Vector& epetraInput = dynamic_cast<const Vector&>(input);
  Vector& epetraResult = dynamic_cast<Vector&>(result);
  
  return applyRightPreconditioning(params, epetraInput, epetraResult);
}

Abstract::Group::ReturnType Group::applyRightPreconditioning(Parameter::List& params,
				      const Vector& input, 
				      Vector& result) const
{
  int errorCode = 1;

  if (preconditioner == "None") {
    result = input;
    return Abstract::Group::Ok;
  }
  else if ((preconditioner == "AztecOO: Jacobian Matrix") || 
	   (preconditioner == "AztecOO: User RowMatrix")) {

    // RPP: We can not directly access Aztec preconditioners.
    // A cheesy way to apply an aztec preconditioner to an arbitrary 
    // vector is to call a Aztec.iterate() but only take one GMRES iteration.
    // This does NOT give the exact preconditioner, but it is a good
    // approximation.  We implement this here but highly recommend the 
    // use of IFPACK preconditioners if available!  

    // Create the temporary vector if necessary 
    if (tmpVectorPtr == 0) {
      tmpVectorPtr = new Epetra_Vector(RHSVector.getEpetraVector());
    }

    // Zero out the temporary vector
    tmpVectorPtr->PutScalar(0.0);

    // Create the solver
    if (aztecSolver == 0)
      aztecSolver = new AztecOO();

    // If a preconditioner is not constructed, create one.
    if (!isPreconditioner()) {

      // Set the preconditioning options
      setAztecOptions(params, *aztecSolver);

      // Get the Jacobian operator/matrix
      Epetra_Operator* matrix = &(sharedJacobian.getOperator(this));
      
      // Try to cast to a Epetra_RowMatrix
      Epetra_RowMatrix* rowMatrix= dynamic_cast<Epetra_RowMatrix*>(matrix);

      if (rowMatrix != 0) 
	aztecSolver->SetUserMatrix(rowMatrix);
      else
	aztecSolver->SetUserOperator(matrix);

      aztecSolver->SetLHS(tmpVectorPtr);
      aztecSolver->SetRHS(tmpVectorPtr);

      createPreconditioner(ExplicitConstruction, params);
    }

    // Turn off printing in Aztec when using applyRightPreconditioner
    aztecSolver->SetAztecOption(AZ_output,AZ_none);

    // Get the number of iterations in the preconditioner
    int numIters = params.getParameter("Preconditioner Iterations", 1);
    
    AztecOO_Operator prec(aztecSolver, numIters);
    
    errorCode = prec.ApplyInverse(input.getEpetraVector(), 
				  result.getEpetraVector());
  }
  else if ((preconditioner == "IFPACK: Jacobian Matrix") ||
	   (preconditioner == "IFPACK: User RowMatrix")){

    if (!isPreconditioner())
      createPreconditioner(ExplicitConstruction, params);

    errorCode = ifpackPreconditioner->ApplyInverse(input.getEpetraVector(), 
						   result.getEpetraVector());

  }
  else if (preconditioner == "User Supplied Preconditioner") {

    if (!isPreconditioner())
      createPreconditioner(ExplicitConstruction, params);

    Epetra_Operator& prec = sharedPreconditionerPtr->getOperator(this);

    errorCode = prec.ApplyInverse(input.getEpetraVector(), 
				  result.getEpetraVector());
  }
  else { 
    cout << "ERROR: NOX::Epetra::Group::applyRightPreconditioning() - "
	 << "Parameter \"preconditioner\" is not vaild for this method"
	 << endl;
    throw "NOX Error";
  }

  if (errorCode == 0) 
    return Abstract::Group::Ok;
  else
    return Abstract::Group::Failed;
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

bool Group::isPreconditioner() const 
{  
  return isValidPreconditioner;
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

Abstract::Group::ReturnType NOX::Epetra::Group::getNormLastLinearSolveResidual(double& residual) const
{
  // Make sure value is not already calculated
  if (isValidNormNewtonSolveResidual) {
    residual = normNewtonSolveResidual;
    return NOX::Abstract::Group::Ok;
  }
  
  // Otherwise give warning since a Newton direction has not been calculated
  // wrt this solution group
  if (utils.isPrintProcessAndType(Utils::Warning)) {
    cout << "ERROR: NOX::Epetra::Group::getNormLastLinearSolveResidual() - "
	 << "Group has not performed a Newton solve corresponding to this "
	 << "solution vector!" << endl;
  }
  return NOX::Abstract::Group::BadDependency;
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

  if ((preconditioner == "AztecOO: Jacobian Matrix") || 
      (preconditioner == "IFPACK: Jacobian Matrix")) {

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
  else if ((preconditioner == "AztecOO: User RowMatrix") || 
	   (preconditioner == "IFPACK: User RowMatrix")){
    
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
    // Do nothing. This is a valid argument and won't cause a throw but 
    // we still want to catch mis-spelled "Preconditioner" flags (see 
    // next "else").
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

bool Group::createPreconditioner(PrecConstructionType ct,
				 Parameter::List& p) const
{
  if (preconditioner == "None") {
    return true;
  }
  else if (preconditioner == "AztecOO: Jacobian Matrix") {
    // The Jacobian has already been evaluated at the current solution.
    // Just enforce explicit constuction if requested
    if (ct == ExplicitConstruction) {
      double conditionNumberEstimate;
      aztecSolver->ConstructPreconditioner(conditionNumberEstimate);
      isValidPreconditioner = true;
    }
  }
  else if (preconditioner == "AztecOO: User RowMatrix") {
    
    Epetra_RowMatrix& precMatrix = 
      dynamic_cast<Epetra_RowMatrix&>(sharedPreconditioner.getOperator(this));
    
    if (preconditionerOperatorType == NoxFiniteDifferenceRowMatrix) {
      (dynamic_cast<FiniteDifference&>(precMatrix))
	.computePreconditioner(xVector.getEpetraVector(), precMatrix);
    }
    else {
      userInterface.computePrecMatrix(xVector.getEpetraVector(), precMatrix);
    }    
    
    aztecSolver->SetPrecMatrix(&precMatrix);

    if (ct == ExplicitConstruction) {
      double conditionNumberEstimate;
      aztecSolver->ConstructPreconditioner(conditionNumberEstimate);
      isValidPreconditioner = true;
    }

  }
  else if (preconditioner == "IFPACK: Jacobian Matrix") {
    
    createIfpackPreconditioner(p);
    isValidPreconditioner = true;
    if (aztecSolver != 0)
      aztecSolver->SetPrecOperator(ifpackPreconditioner);
    
  }
  else if (preconditioner == "IFPACK: User RowMatrix") {

    Epetra_RowMatrix& precMatrix = 
      dynamic_cast<Epetra_RowMatrix&>(sharedPreconditioner.getOperator(this));
    
    if (preconditionerOperatorType == NoxFiniteDifferenceRowMatrix) {
      (dynamic_cast<FiniteDifference&>(precMatrix))
	.computePreconditioner(xVector.getEpetraVector(), precMatrix);
    }
    else {
      userInterface.computePrecMatrix(xVector.getEpetraVector(), precMatrix);
    }    

    createIfpackPreconditioner(p);
    isValidPreconditioner = true;
    if (aztecSolver != 0)
      aztecSolver->SetPrecOperator(ifpackPreconditioner);

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
    if (aztecSolver != 0)
      aztecSolver->SetPrecOperator(&precOperator);
    isValidPreconditioner = true; 

  }
  else {
    // The "preconditioner" string is invalid
    cout << "ERROR: NOX::Epetra::Group::createPreconditioner() - the choice "
	 << "for \"Preconditioning\" is invalid" << endl;
    throw "NOX Error";
  }

  return true;
}

bool Group::createIfpackPreconditioner(Parameter::List& p) const
{
  //for ifpack we need a VBR or CRS matrix to get the correct graph
  const Epetra_Operator* op = 0;

  // Get the correct matrix for preconditioning
  if (preconditioner == "IFPACK: Jacobian Matrix") {
    op = &sharedJacobian.getOperator();
  }
  else if (preconditioner == "IFPACK: User RowMatrix") {
    op = &sharedPreconditioner.getOperator();
  }
  else {
    cout << "ERROR: NOX::Epetra::Group::createIfpackPreconditioner() - "
	 << " preconditioner choice is not valid!" << endl;
    throw "NOX Error";
  }

  //check to see if it is a VBR matrix
  const Epetra_VbrMatrix* vbr= dynamic_cast<const Epetra_VbrMatrix*>(op);
  if (vbr != 0) {
    ifpackGraph = new Ifpack_IlukGraph(vbr->Graph(),
				       p.getParameter("Fill Factor", 1),
				       p.getParameter("Overlap", 0));
    ifpackGraph->ConstructFilledGraph();
    ifpackPreconditioner = new Ifpack_CrsRiluk(*ifpackGraph);
    ifpackPreconditioner->InitValues(*vbr);
    ifpackPreconditioner->Factor();
    return true;
  }

  // check to see if it is a Crs matrix
  const Epetra_CrsMatrix* crs= dynamic_cast<const Epetra_CrsMatrix*>(op);
  if (crs != 0) {

    if (ifpackGraph != 0) cout << "Ifpack Graph NOT NULL" << endl;
    if (ifpackPreconditioner != 0) cout << "Ifpack Prec NOT NULL" << endl;


    ifpackGraph = new Ifpack_IlukGraph(crs->Graph(),
				       p.getParameter("Fill Factor", 1),
				       p.getParameter("Overlap", 0));
    ifpackGraph->ConstructFilledGraph();
    ifpackPreconditioner = new Ifpack_CrsRiluk(*ifpackGraph);
    ifpackPreconditioner->InitValues(*crs);
    ifpackPreconditioner->Factor();
    return true;
  }
  
  // If we made it this far, this routine should not have been called
  // in the first place.  An incorrect prec matrix object was supplied and
  // should have been caught in the checkOperatorConsistency() method.
  cout << "ERROR: NOX::Epetra::Group::createIfpackPreconditioner() - "
       << " preconditioner object is not a CRS or VBR matrix!" << endl;
  throw "NOX Error";

  return false;
}

bool Group::destroyPreconditioner() const
{
  // For Aztec preconditioners.  
  if (isPreconditioner()) {
    if ((preconditioner == "AztecOO: Jacobian Matrix") ||
	(preconditioner == "AztecOO: User RowMatrix")) {
      if (aztecSolver != 0)
      aztecSolver->DestroyPreconditioner();
    }
  }

  // For IFPACK Preconditioners
  delete ifpackPreconditioner;
  ifpackPreconditioner = 0;
  delete ifpackGraph;
  ifpackGraph = 0;
   
  isValidPreconditioner = false;

  return true;
}

bool Group::destroyAztecSolver() const {
  destroyPreconditioner();
  delete aztecSolver;
  aztecSolver = 0;
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
