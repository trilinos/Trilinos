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
#include "NOX_Epetra_Preconditioner.H"

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
  xVector(x), // deep copy x     
  RHSVector(x, ShapeCopy), // new vector of same size
  gradVector(x, ShapeCopy), // new vector of same size
  NewtonVector(x, ShapeCopy), // new vector of same size
  tmpVectorPtr(0),
  sharedJacobianPtr(new SharedOperator(J)), // pass J to SharedJacobian
  sharedJacobian(*sharedJacobianPtr), // pass J to SharedJacobian
  sharedPreconditionerPtr(0),  // A separate preconditioner is not supplied
  userInterface(i)
{
  resetIsValid();

  setLinearSolver(params);
}

Group::Group(const Parameter::List& params, Interface& i, 
	     Epetra_Vector& x, Epetra_Operator& J, Epetra_Operator& M):
  xVector(x), // deep copy x     
  RHSVector(x, ShapeCopy), // new vector of same size
  gradVector(x, ShapeCopy), // new vector of same size
  NewtonVector(x, ShapeCopy), // new vector of same size
  tmpVectorPtr(NULL),
  sharedJacobianPtr(new SharedOperator(J)), // pass J to SharedOperator
  sharedJacobian(*sharedJacobianPtr), // create reference from pointer
  sharedPreconditionerPtr(new SharedOperator(M)), // pass M to SharedOperator
  userInterface(i)
{
  resetIsValid();

  setLinearSolver(params);
}

Group::Group(const Group& source, CopyType type) :
  xVector(source.xVector.getEpetraVector(), type), 
  RHSVector(source.RHSVector.getEpetraVector(), type), 
  gradVector(source.gradVector.getEpetraVector(), type), 
  NewtonVector(source.NewtonVector.getEpetraVector(), type),
  tmpVectorPtr(NULL),
  sharedJacobianPtr(NULL),
  sharedJacobian(source.sharedJacobian),
  sharedPreconditionerPtr(source.sharedPreconditionerPtr),
  userInterface(source.userInterface),
  jacType(source.jacType),
  precOption(source.precOption),
  precType(source.precType),
  aztecPrec(source.aztecPrec)
{
 
  switch (type) {
    
  case DeepCopy:
    
    isValidRHS = source.isValidRHS;
    isValidGrad = source.isValidGrad;
    isValidNewton = source.isValidNewton;
    isValidJacobian = source.isValidJacobian;
    isValidPreconditioner = source.isValidPreconditioner;
    normRHS = source.normRHS;
    
    // New copy takes ownership of the shared Jacobian
    if (isValidJacobian)
      sharedJacobian.getOperator(this);

    // New copy takes ownership of the shared Preconditioning Matrix
    if (isValidPreconditioner)
      sharedPreconditionerPtr->getOperator(this);

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
}

void Group::resetIsValid() //private
{
  isValidRHS = false;
  isValidJacobian = false;
  isValidGrad = false;
  isValidNewton = false;
  isValidPreconditioner = false;
}

void Group::setAztecOptions(const Parameter::List& p, AztecOO& aztec)
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
  if (precOption == "None")
    aztec.SetAztecOption(AZ_precond, AZ_none);

  // Preconditioning where AztecOO inverts the Preconditioning Matrix
  if ((precOption == "AztecOO: Jacobian Matrix") || 
      (precOption == "AztecOO: User RowMatrix")) {
    
    string aztecPrec;
    aztecPrec = p.getParameter("Aztec Preconditioner", "ilu");
    
    if (aztecPrec == "ilu") {
      aztec.SetAztecOption(AZ_precond, AZ_dom_decomp);
      aztec.SetAztecOption(AZ_overlap, p.getParameter("Overlap", 0));
      aztec.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
      aztec.SetAztecOption(AZ_graph_fill, p.getParameter("Graph Fill", 0));
    }
    else if (aztecPrec == "ilut") { 
      aztec.SetAztecOption(AZ_precond, AZ_dom_decomp);
      aztec.SetAztecOption(AZ_overlap, p.getParameter("Overlap", 0));
      aztec.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
      aztec.SetAztecParam(AZ_drop, p.getParameter("Drop Tolerance", 0.0));
      aztec.SetAztecParam(AZ_ilut_fill, p.getParameter("Fill Factor", 1.0));
    }
    else if (aztecPrec == "Jacobi") {
      aztec.SetAztecOption(AZ_precond, AZ_Jacobi);
      aztec.SetAztecOption(AZ_poly_ord, p.getParameter("Steps", 3));
    }
    else if (aztecPrec == "Symmetric Gauss-Siedel") {
      aztec.SetAztecOption(AZ_precond, AZ_sym_GS);
      aztec.SetAztecOption(AZ_poly_ord, p.getParameter("Steps", 3));
    }
    else if (aztecPrec == "Polynomial") {
      aztec.SetAztecOption(AZ_precond, AZ_Neumann);
      aztec.SetAztecOption(AZ_poly_ord, p.getParameter("Polynomial Order", 3));
    }
    else if (aztecPrec == "Least-squares Polynomial") {
      aztec.SetAztecOption(AZ_precond, AZ_ls);
      aztec.SetAztecOption(AZ_poly_ord, p.getParameter("Polynomial Order", 3));
    }
    else {
      cout << "ERROR: NOX::Epetra::Group::setAztecOptions" << endl
	   << "\"Aztec Preconditioner\" parameter \"" << aztecPrec
	   << "\" is invalid!" << endl;
      throw "NOX Error";
    }

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
  isValidJacobian = source.isValidPreconditioner;

  // Only copy vectors that are valid
  if (isValidRHS) {
    RHSVector = source.RHSVector;
    normRHS = source.normRHS;
  }

  if (isValidGrad)
    gradVector = source.gradVector;

  if (isValidNewton)
    NewtonVector = source.NewtonVector;

  // If valid, this takes ownership of the shared Jacobian
  if (isValidJacobian)
    sharedJacobian.getOperator(this);
    
  // If valid, this takes ownership of the shared preconditioning matrix
  if (isValidPreconditioner)
    sharedJacobian.getOperator(this);
    
  // Copy linear solver options
  jacType = source.jacType;
  precOption = source.precOption;
  precType = source.precType;
  aztecPrec = source.aztecPrec;

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
  // Need to add a check here to find out if computing the Jacobian was successful.  
  bool status = false;
  
  if (jacType == EpetraOperatorJac) {
    status = userInterface.computeJacobian(xVector.getEpetraVector(), Jacobian);
  }
  else if (jacType == NoxFiniteDifferenceJac) {
    status = (dynamic_cast<FiniteDifference&>(Jacobian)).computeJacobian(xVector.getEpetraVector(), Jacobian);
  }
  else if (jacType == NoxMatrixFreeJac) {
    status = (dynamic_cast<MatrixFree&>(Jacobian)).computeJacobian(xVector.getEpetraVector(), Jacobian);
  }
  else {
    cout << "ERROR: NOX::Epetra::Group::computeJacobian() - jacType flag is "
	 << "not valid!" << endl;
    throw "NOX Error";
  }

  if (status == false) {
    cout << "ERROR: Epetra::Group::computeJacobian() - fill failed!!!"
	 << endl;
    throw "NOX Error: Fill Failed";
  } 

  // Update status
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
  
  // Get the Jacobian 
  /* (Have to get non-const version which requires reasserting
   * ownership. This is due to the design of Epetra_LinearProblem class). 
   */
  Epetra_Operator& Jacobian = sharedJacobian.getOperator(this);

  // Create Epetra problem for the linear solve
  Epetra_LinearProblem Problem(&Jacobian, 
  			       &(NewtonVector.getEpetraVector()), 
			       &(RHSVector.getEpetraVector()));

  // Set the default Problem parameters to "hard" (this sets Aztec defaults
  // during the AztecOO instantiation)
  Problem.SetPDL(hard);

  // Create aztec problem
  AztecOO aztec(Problem);    

  // Set specific Aztec parameters requested by NOX
  setAztecOptions(p, aztec);

  // Compute and set the Preconditioner in AztecOO if needed
  if (precOption != "None") {
    computePreconditioner(aztec);
  }

  int maxit = p.getParameter("Max Iterations", 400);
  double tol = p.getParameter("Tolerance", 1.0e-6);

  // Solve Aztec problem
  aztec.Iterate(maxit, tol);

  // Scale soln by -1
  NewtonVector.scale(-1.0);

  // Update state
  isValidNewton = true;

  // Return solution
  return true;
}

bool Group::computePreconditioner()
{
  cout << "NOX::Epetra::Group::computePreconditioner() NOT Implemented yet!!" << endl;
  exit(0);

  /* RPP: This will come online soon 

  // Skip if the preconditioning matrix is already valid
  if (isPreconditioner())
    return true;

  // Take ownership of the Preconditioner and get a reference to the underlying operator
  Epetra_Operator& prec = sharedPreconditioner.getOperator(this);

  // Fill the Preconditioning matrix

  //RPP: Move over the computePreconditioner() method into here!!

  */
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

bool Group::applyPreconditionerInverse(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& epetrainput = dynamic_cast<const Vector&> (input);
  Vector& epetraresult = dynamic_cast<Vector&> (result);
  return applyPreconditionerInverse(epetrainput, epetraresult);
}

bool Group::applyPreconditionerInverse(const Vector& input, Vector& result) const
{
  return applyJacobianDiagonalInverse(input, result);

  /* RPP: Stuff below will come online soon 
  if (isPreconditioner()) {
    NOX::Epetra::Preconditioner* prec = dynamic_cast<NOX::Epetra::Preconditoner*>(&(sharedPreconditioner.getOperator()));
    return prec.ApplyInverse(input.getEpetraVector(), result.getEpetraVector());
  }
  else {
    cout << "ERROR: NOX::Epetra::Group::applyPreconditionerInverse() - is only "
	 << "implemented for a User Supplied Preconditioner through the "
	 << "NOX::Epetra::Preconditioner interface!" << endl;
    throw "NOX Error";
  }
  return false;
  */
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


SharedOperator& Group::getSharedJacobian()
{
  return sharedJacobian;
}

SharedOperator& Group::getSharedPreconditioner()
{
  return *sharedPreconditionerPtr;
}

Interface& Group::getUserInterface()
{
  return userInterface;
}

bool Group::setLinearSolver(const Parameter::List& params) 
{

  // Set the type of the Jacobian operator 
  jacType = getJacobianType();

  // Set the requested preconditioning option.  Defaults to "None".
  precOption = params.getParameter("Preconditioning", "None");

  // Check that the requested preconditioning option (precOption) has 
  // the correct objects avalable 
  Epetra_RowMatrix* test = 0;

  if (precOption == "AztecOO: Jacobian Matrix") {

    // Get a reference to the Jacobian 
    Epetra_Operator& Jacobian = sharedJacobian.getOperator(this);

    // The Jacobian Operator MUST be an Eptra_RowMatrix
    test = dynamic_cast<Epetra_RowMatrix*>(&Jacobian);

    if (test == 0) {
      cout << "ERROR: NOX::Epetra::Group::setLinearSolver() - The flag "
	   << "\"Preconditioning\" with value \"" << precOption
	   << "\" requires an Epetra_RowMatrix object for the Jacobian!" 
	   << endl;
      cout << "Matrx-Free Jacobians can not be used for preconditioners!"
	   << endl;
      throw "NOX Error";
    }
  }
  else if (precOption == "AztecOO: User RowMatrix") {
    
    // Set the type of Preconditioning operator
    precType = getPrecType();

    // Get a reference to the Preconditioner 
    /* NOTE:SharedJacobian.getPrec() print a WARNING if no 
     * preconditioner was supplied.  We'll throw an error and explain why. 
     */   
    Epetra_Operator& Prec = sharedPreconditionerPtr->getOperator(this);

    // The Preconditioning Operator MUST BE an Epetra_RowMatrix
    test = dynamic_cast<Epetra_RowMatrix*>(&Prec);
    
    if (test == 0) {
      cout << "ERROR: NOX::Epetra::Group::setLinearSolver() - The flag "
	   << "\"Preconditioning\" with value \"" << precOption
	   << "\" requires a NOX::Epetra::Group constructor with "
	   << "separate objects for the Jacobian and preconditioner!" << endl;
      cout << "Make sure the preconditioner is an Epetra_RowMatrix object!"
	   << endl;
      throw "NOX Error";
    }
  }
  else if (precOption == "User Supplied Preconditioner") {

    // Set the type of Preconditioning operator
    precType = getPrecType();

    // Get a reference to the Preconditioner
    /* NOTE:SharedJacobian.getPrec() print a WARNING if no 
     * preconditioner was supplied.  We'll throw an error and explain why. 
     */   
    const Epetra_Operator& Prec = sharedPreconditionerPtr->getOperator();
    if (&Prec == 0) {
      cout << "ERROR: NOX::Epetra::Group::setLinearSolver() - The flag "
	   << "\"Preconditioning\" with value \"" << precOption
	   << "\" requires a NOX::Epetra::Group constructor with "
	   << "separate objects for the Jacobian and preconditioner!" << endl; 
      throw "NOX Error";
    }
    
  }
  else if (precOption == "None") {
    // Do nothing. This is a valid argument and won't cause a throw - 
    // see next "else" statement.
  } 
  else {
    // An invalid choice was found 
    cout << "ERROR: NOX::Epetra::Group::setLinearSolver() - The flag "
	 << "\"" << precOption << "\" is not a valid choice for the " 
	 << "parameter \"Preconditioning\"!" << endl;
      throw "NOX Error";
  }
  
  // All checks passed without a throw!
  return true;
}

bool Group::computePreconditioner(AztecOO& aztec)
{
  
  if (precOption == "AztecOO: Jacobian Matrix") {
    // Do nothing, the Jacobian has already been evaluated at the
    // current solution.
  }
  else if (precOption == "AztecOO: User RowMatrix") {
    
    Epetra_RowMatrix& precMatrix = 
      dynamic_cast<Epetra_RowMatrix&>(sharedPreconditionerPtr->getOperator(this));
    
    if (precType == NoxFiniteDifferencePrec) {
      (dynamic_cast<FiniteDifference&>(precMatrix))
	.computePreconditioner(xVector.getEpetraVector(), precMatrix);
    }
    else 
      userInterface.computePreconditioner(xVector.getEpetraVector(), precMatrix);
    
    aztec.SetPrecMatrix(&precMatrix);
  }
  else if (precOption == "User Supplied Preconditioner") {
    
    Epetra_Operator& precOperator = sharedPreconditionerPtr->getOperator(this);

    // Fro user supplied Precs the user may supply them in one of two ways:
    if (precType == NoxPreconditioner) {
      // (1) Use a NOX::Epetra::Preconditioner derived object
      const Epetra_Operator* jacobianPtr = 0; 
      if (isJacobian())
	jacobianPtr = &sharedJacobian.getOperator();

      (dynamic_cast<NOX::Epetra::Preconditioner&>(precOperator))
	.computePreconditioner(xVector.getEpetraVector(), jacobianPtr);
    }
    else {
      // (2) Supply an Epetra_Operator derived class and implement 
      // the NOX::Epetra::Interface::computePreconditioner() method
      userInterface.computePreconditioner(xVector.getEpetraVector(),
					  precOperator);
    }
    aztec.SetPrecOperator(&precOperator); 
  }

  return true;
}

Group::JacType Group::getJacobianType()
{
  
  // Get a reference to the Jacobian 
  Epetra_Operator& Jacobian = sharedJacobian.getOperator(this);

  Epetra_Operator* test = 0;
  
  // Is it a Matrix-Free Jacobian?
  test = dynamic_cast<MatrixFree*>(&Jacobian);
  if (test != 0) 
    return NoxMatrixFreeJac; 
  
  // Is it a Finite Difference Jacobian?
  test = dynamic_cast<FiniteDifference*>(&Jacobian);
  if (test != 0) 
    return NoxFiniteDifferenceJac; 
  
  // Otherwise it must be a user supplied Epetra_Operator or Epetra_RowMatrix.
  return EpetraOperatorJac;
  
}

Group::PrecType Group::getPrecType()
{
  // Get a reference to the preconditioner 
  Epetra_Operator& precMatrix = sharedPreconditionerPtr->getOperator(this);

  Epetra_Operator* test = 0;
  
  // Is it a "Finite Difference" object?
  test = dynamic_cast<FiniteDifference*>(&precMatrix);
  if (test != 0) 
    return NoxFiniteDifferencePrec; 
  
  // Is it a NOX::Epetra::Preconditioner object?
  test = dynamic_cast<NOX::Epetra::Preconditioner*>(&precMatrix);
  if (test != 0) 
    return NoxPreconditioner; 
  
  // Otherwise it must be a user supplied Epetra_Operator or  Epetra_RowMatrix.
  return EpetraOperatorPrec;
}
