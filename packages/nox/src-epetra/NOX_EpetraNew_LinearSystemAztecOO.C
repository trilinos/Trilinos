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

#include "NOX_EpetraNew_LinearSystemAztecOO.H"	// class definition

// NOX includes
#include "NOX_EpetraNew_Interface_Required.H"
#include "NOX_EpetraNew_Interface_Jacobian.H"
#include "NOX_EpetraNew_Interface_Preconditioner.H"
#include "NOX_EpetraNew_MatrixFree.H"
#include "NOX_EpetraNew_FiniteDifference.H"
#include "NOX_Parameter_List.H"
#include "NOX_EpetraNew_Scaling.H"
#include "NOX_Utils.H"

// External include files for Epetra, Aztec00, and Ifpack
#include "Epetra_Map.h"
#include "Epetra_Vector.h" 
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "AztecOO_Operator.h"
#include "AztecOO_StatusTest.h"
#include "AztecOO_StatusTestCombo.h"
#include "AztecOO_StatusTestMaxIters.h"
#include "AztecOO_StatusTestResNorm.h"
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"

#ifdef HAVE_NOX_ML_EPETRA
#include "Teuchos_ParameterList.hpp"
#endif

#include <typeinfo>

//***********************************************************************
NOX::EpetraNew::LinearSystemAztecOO::
LinearSystemAztecOO(NOX::Parameter::List& printParams, 
		    NOX::Parameter::List& linearSolverParams, 
		    NOX::EpetraNew::Interface::Required& iReq, 
		    Epetra_Vector& cloneVector,
		    NOX::EpetraNew::Scaling* s):
  utils(printParams),
  jacInterfacePtr(0),
  jacType(EpetraOperator),
  jacPtr(0),
  ownsJacOperator(false),
  precInterfacePtr(0),
  precType(EpetraOperator),
  precPtr(0),
  ownsPrecOperator(false),
  precMatrixSource(UseJacobian),
  aztecSolverPtr(new AztecOO()),
  ifpackGraphPtr(0),
  ifpackPreconditionerPtr(0),
#ifdef HAVE_NOX_ML_EPETRA
  MLPreconditionerPtr(0),
#endif
  scaling(s),
  tmpVectorPtr(new Epetra_Vector(cloneVector)),
  conditionNumberEstimate(0.0),
  isPrecConstructed(false),
  precQueryCounter(0),
  maxAgeOfPrec(1),
  timer(cloneVector.Comm()),
  timeCreatePreconditioner(0.0),
  timeApplyJacbianInverse(0.0)
{
  // Neither Jacobian or Preconditioner are supplied
  createJacobianOperator(linearSolverParams, iReq, cloneVector);
  createPrecOperator(linearSolverParams, iReq, cloneVector);

  reset(linearSolverParams);
}

//***********************************************************************
NOX::EpetraNew::LinearSystemAztecOO::
LinearSystemAztecOO(NOX::Parameter::List& printParams, 
		    NOX::Parameter::List& linearSolverParams,  
		    NOX::EpetraNew::Interface::Required& iReq, 
		    NOX::EpetraNew::Interface::Jacobian& iJac, 
		    Epetra_Operator& jacobian,
		    Epetra_Vector& cloneVector,
		    NOX::EpetraNew::Scaling* s):
  utils(printParams),
  jacInterfacePtr(&iJac),
  jacType(EpetraOperator),
  jacPtr(&jacobian),
  ownsJacOperator(false),
  precInterfacePtr(0),
  precType(EpetraOperator),
  precPtr(0),
  ownsPrecOperator(false),
  precMatrixSource(UseJacobian),
  aztecSolverPtr(new AztecOO()),
  ifpackGraphPtr(0),
  ifpackPreconditionerPtr(0),
#ifdef HAVE_NOX_ML_EPETRA
  MLPreconditionerPtr(0),
#endif
  scaling(s),
  tmpVectorPtr(new Epetra_Vector(cloneVector)),
  conditionNumberEstimate(0.0),
  isPrecConstructed(false),
  precQueryCounter(0),
  maxAgeOfPrec(1),
  timer(cloneVector.Comm()),
  timeCreatePreconditioner(0.0),
  timeApplyJacbianInverse(0.0)
{
  // Jacobian operator is supplied 
  jacType = getOperatorType(*jacPtr);
  // Preconditioner is not supplied
  createPrecOperator(linearSolverParams, iReq, cloneVector);

  reset(linearSolverParams);
}

//***********************************************************************
NOX::EpetraNew::LinearSystemAztecOO::
LinearSystemAztecOO(NOX::Parameter::List& printParams, 
		    NOX::Parameter::List& linearSolverParams, 
		    NOX::EpetraNew::Interface::Required& iReq, 
		    NOX::EpetraNew::Interface::Preconditioner& iPrec, 
		    Epetra_Operator& preconditioner,
		    Epetra_Vector& cloneVector,
		    NOX::EpetraNew::Scaling* s):
  utils(printParams),
  jacInterfacePtr(0),
  jacType(EpetraOperator),
  jacPtr(0),
  ownsJacOperator(false),
  precInterfacePtr(&iPrec),
  precType(EpetraOperator),
  precPtr(&preconditioner),
  ownsPrecOperator(false),
  precMatrixSource(SeparateMatrix),
  aztecSolverPtr(new AztecOO()),
  ifpackGraphPtr(0),
  ifpackPreconditionerPtr(0),
#ifdef HAVE_NOX_ML_EPETRA
  MLPreconditionerPtr(0),
#endif
  scaling(s),
  tmpVectorPtr(new Epetra_Vector(cloneVector)),
  conditionNumberEstimate(0.0),
  isPrecConstructed(false),
  precQueryCounter(0),
  maxAgeOfPrec(1),
  timer(cloneVector.Comm()),
  timeCreatePreconditioner(0.0),
  timeApplyJacbianInverse(0.0)
{
  // Jacobian operator is not supplied
  createJacobianOperator(linearSolverParams, iReq, cloneVector);
  // Preconditioner operator is supplied
  precType = getOperatorType(*precPtr);

  reset(linearSolverParams);
}

//***********************************************************************
NOX::EpetraNew::LinearSystemAztecOO::
LinearSystemAztecOO(NOX::Parameter::List& printParams, 
		    NOX::Parameter::List& linearSolverParams,
		    NOX::EpetraNew::Interface::Jacobian& iJac, 
		    Epetra_Operator& jacobian,
		    NOX::EpetraNew::Interface::Preconditioner& iPrec, 
		    Epetra_Operator& preconditioner,
		    Epetra_Vector& cloneVector,
		    NOX::EpetraNew::Scaling* s):
  utils(printParams),
  jacInterfacePtr(&iJac),
  jacType(EpetraOperator),
  jacPtr(&jacobian),
  ownsJacOperator(false),
  precInterfacePtr(&iPrec),
  precType(EpetraOperator),
  precPtr(&preconditioner),
  ownsPrecOperator(false),
  precMatrixSource(SeparateMatrix),
  aztecSolverPtr(new AztecOO()),
  ifpackGraphPtr(0),
  ifpackPreconditionerPtr(0),
#ifdef HAVE_NOX_ML_EPETRA
  MLPreconditionerPtr(0),
#endif
  scaling(s),
  tmpVectorPtr(new Epetra_Vector(cloneVector)),
  conditionNumberEstimate(0.0),
  isPrecConstructed(false),
  precQueryCounter(0),
  maxAgeOfPrec(1),
  timer(cloneVector.Comm()),
  timeCreatePreconditioner(0.0),
  timeApplyJacbianInverse(0.0)
{
  // Both operators are supplied
  jacType = getOperatorType(*jacPtr);
  precType = getOperatorType(*precPtr);

  reset(linearSolverParams);
}

//***********************************************************************
NOX::EpetraNew::LinearSystemAztecOO::~LinearSystemAztecOO() 
{
  destroyPreconditioner();
  if (ownsJacOperator) {
    delete jacPtr;
    jacPtr = 0;
  }
  if (ownsPrecOperator) {
    delete precPtr;
    precPtr = 0;
  }
  delete aztecSolverPtr;
  delete tmpVectorPtr;
}

//***********************************************************************
void NOX::EpetraNew::LinearSystemAztecOO::
reset(NOX::Parameter::List& linearSolverParams)
{
  // First remove any preconditioner that may still be active
  destroyPreconditioner();

  // Set the requested preconditioning.
  string prec = linearSolverParams.getParameter("Preconditioner", "None");
  if (prec == "AztecOO")
    precAlgorithm = AztecOO_;
  else if (prec == "Ifpack")
    precAlgorithm = Ifpack_;
#ifdef HAVE_NOX_ML_EPETRA
  else if (prec == "ML")
    precAlgorithm = ML_;
#endif
  else if (prec == "User Defined")
    precAlgorithm = UserDefined_;
  else if (prec == "None") 
    precAlgorithm = None_;
  else {
    string errorMessage = "Option for \"Preconditioner\" is invalid!";
    throwError("LinearSystemAztecOO(J)", errorMessage);
  }
    
  // Make sure the correct objects were supplied for the requested
  // preconditioning choice.
  checkPreconditionerValidity();

  zeroInitialGuess = 
    linearSolverParams.getParameter("Zero Initial Guess", false);

  // Place linear solver details in the "Output" sublist of the
  // "Linear Solver" parameter list
  outputSolveDetails = 
    linearSolverParams.getParameter("Output Solver Details", true);

  // The first time a SetProblem is used on the AztecOO solver
  // it sets all aztec options based on the Epetra_LinearProblem
  // options. Subsequent calls do not.  We call this here so we
  // can set our own parameters and not have them overwritten 
  // by the first SetProblem call.
  // RPP: Not any more.  We don't set the solver objects with
  // the problem class.  It cause seg faults when new preconditioners
  // were computed during prec reuse.
  //Epetra_LinearProblem& problem = *(new Epetra_LinearProblem);
  //aztecSolverPtr->SetProblem(problem);

  // Set the Jacobian in the solver. It must be set before
  // a preconditioner can be set.
//   if ((jacType == EpetraRowMatrix) ||
//       (jacType == EpetraVbrMatrix) ||
//       (jacType == EpetraCrsMatrix)) {
//     aztecSolverPtr->SetUserMatrix(dynamic_cast<Epetra_RowMatrix*>(jacPtr));
//   }
//   else
//     aztecSolverPtr->SetUserOperator(jacPtr);

  // Set the major aztec options.  Must be called after the first 
  // SetProblem() call.
  setAztecOptions(linearSolverParams, *aztecSolverPtr);

  maxAgeOfPrec = linearSolverParams.getParameter("Max Age Of Prec", 1);
}

//***********************************************************************
void NOX::EpetraNew::LinearSystemAztecOO::
setAztecOptions(const Parameter::List& p, AztecOO& aztec) const
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
    cout << "ERROR: NOX::EpetraNew::Group::setAztecOptions" << endl
	 << "\"Aztec Solver\" parameter \"" << linearSolver 
	 <<  "\" is invalid!" << endl;
    throw "NOX Error";
  }
 
  // Preconditioning where AztecOO inverts the Preconditioning Matrix
  if (precAlgorithm == AztecOO_) {
    
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
      string errorMessage = "\"Aztec Preconditioner\" parameter is invalid!";
      throwError("setAztecOptions", errorMessage);
    }

  }
  else 
    aztec.SetAztecOption(AZ_precond, AZ_none);
    
  // Turn on RCM reordering in conjunction with domain decomp preconditioning
  // default is "Disabled" = no reordering
  string rcmReordering = p.getParameter("RCM Reordering", "Disabled");
  if (rcmReordering == "Enabled")
    aztec.SetAztecOption(AZ_reorder, 1);
  else if (rcmReordering == "Disabled")
    aztec.SetAztecOption(AZ_reorder, 0);
  else {
    string errorMessage = "\"RCM Reordering\" parameter is invalid!";
    throwError("setAztecOptions", errorMessage);
  }
    
  // Gram-Schmidt orthogonalization procedure
  string orthog = p.getParameter("Orthogonalization", "Classical");
  if (orthog == "Classical") 
    aztec.SetAztecOption(AZ_orthog, AZ_classic);
  else if (orthog == "Modified")
    aztec.SetAztecOption(AZ_orthog, AZ_modified);
  else {
    string errorMessage = "\"Orthogonalization\" parameter is invalid!";
    throwError("setAztecOptions()", errorMessage);
  }

  // Size of the krylov subspace
  aztec.SetAztecOption(AZ_kspace, 
		       p.getParameter("Size of Krylov Subspace", 300));

  // Convergence criteria to use in the linear solver
  string convCriteria = p.getParameter("Convergence Test", "r0");
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
    string errorMessage = "\"Convergence Test\" parameter is invalid!";
    throwError("setAztecOptions()", errorMessage);
  }

  // Set the ill-conditioning threshold for the upper hessenberg matrix
  if (p.isParameter("Ill-Conditioning Threshold")) {
    aztec.SetAztecParam(AZ_ill_cond_thresh, 
			p.getParameter("Ill-Conditioning Threshold", 1.0e+11));
  }

  // Frequency of linear solve residual output
  if (utils.isPrintType(Utils::LinearSolverDetails))
    aztec.SetAztecOption(AZ_output, p.getParameter("Output Frequency", 
						   AZ_last));
  else
    aztec.SetAztecOption(AZ_output, p.getParameter("Output Frequency", 0));

  // Print a summary of the aztec options if "Details" is enabled
  if (utils.isPrintType(Utils::Debug)) {
    //aztec.CheckInput();
  }

  // Some Debugging utilities 
  if (utils.isPrintProcessAndType(Utils::Debug)) {
    cout << "NOX::Epetra::LinearSystemAztecOO Operator Information" << endl;
    cout << "jacType = " << jacType << endl;
    cout << "jacPtr = " << jacPtr << endl;
    cout << "jacInterfacePtr = " << jacInterfacePtr << endl;
    cout << "ownsJacOperator = " << ownsJacOperator << endl;
    cout << "precType = " << precType << endl;
    cout << "precPtr = " << precPtr << endl;
    cout << "precInterfacePtr = " << precInterfacePtr << endl;
    cout << "ownsPrecOperator = " << ownsPrecOperator << endl;
  }
    
  return;
}

//***********************************************************************
bool NOX::EpetraNew::LinearSystemAztecOO::
createJacobianOperator(NOX::Parameter::List& lsParams,
		       NOX::EpetraNew::Interface::Required& iReq, 
		       const Epetra_Vector& cloneVector)
{
  string choice = lsParams.getParameter("Jacobian Operator", "Matrix-Free");

  if (choice == "Matrix-Free") {
    jacPtr = new MatrixFree(iReq, cloneVector);
    jacInterfacePtr = 
      dynamic_cast<NOX::EpetraNew::Interface::Jacobian*>(jacPtr);
    jacType = EpetraOperator;
    ownsJacOperator = true;
  }
  else if (choice == "Finite Difference") {
    jacPtr = new FiniteDifference(iReq, cloneVector);
    jacInterfacePtr = 
      dynamic_cast<NOX::EpetraNew::Interface::Jacobian*>(jacPtr);
    jacType = EpetraRowMatrix;
    ownsJacOperator = true;
  }
  else    
    throwError("createJacobianOperator", 
       "The specified value for parameter \" Jacobian Operator\" is not valid");

  return true;
}

//***********************************************************************
bool NOX::EpetraNew::LinearSystemAztecOO::
createPrecOperator(NOX::Parameter::List& lsParams,
		   NOX::EpetraNew::Interface::Required& iReq, 
		   const Epetra_Vector& cloneVector)
{
  string choice = lsParams.getParameter("Preconditioner Operator", 
					"Use Jacobian");

  if (choice == "Use Jacobian") {
    precPtr = 0;
    precInterfacePtr = 0;
    precType = jacType;
    ownsPrecOperator = false;
    precMatrixSource = UseJacobian;
  }
  else if (choice == "Finite Difference") {
    precPtr = new FiniteDifference(iReq, cloneVector);
    precInterfacePtr = 
      dynamic_cast<NOX::EpetraNew::Interface::Preconditioner*>(precPtr);
    precType = EpetraRowMatrix;
    ownsPrecOperator = true;
    precMatrixSource = SeparateMatrix;
  }
  else    
    throwError("createPreconditionerOperator", 
       "The value for the parameter \" Preconditioner Operator\" is not valid");

  return true;
}

//***********************************************************************
bool NOX::EpetraNew::LinearSystemAztecOO::
applyJacobian(const NOX::Epetra::Vector& input, 
	      NOX::Epetra::Vector& result) const
{
  jacPtr->SetUseTranspose(false);
  int status = jacPtr->Apply(input.getEpetraVector(), 
				  result.getEpetraVector());
  return (status == 0);
}

//***********************************************************************
bool NOX::EpetraNew::LinearSystemAztecOO::
applyJacobianTranspose(const NOX::Epetra::Vector& input, 
		       NOX::Epetra::Vector& result) const
{
  // Apply the Jacobian
  jacPtr->SetUseTranspose(true);
  int status = jacPtr->Apply(input.getEpetraVector(), 
				  result.getEpetraVector());
  jacPtr->SetUseTranspose(false);

  return (status == 0);
}

//***********************************************************************
bool NOX::EpetraNew::LinearSystemAztecOO::
applyJacobianInverse(Parameter::List &p,
		     const NOX::Epetra::Vector& input, 
		     NOX::Epetra::Vector& result)
{
  double startTime = timer.WallTime();

  // Need non-const version of the input vector
  // Epetra_LinearProblem requires non-const versions so we can perform
  // scaling of the linear problem.
  NOX::Epetra::Vector& nonConstInput = const_cast<NOX::Epetra::Vector&>(input);
  
  // Zero out the delta X of the linear problem if requested by user.
  if (zeroInitialGuess)
    result.init(0.0);

  // Create Epetra linear problem object for the linear solve
  Epetra_LinearProblem Problem(solveJacOpPtr, 
  			       &(result.getEpetraVector()), 
			       &(nonConstInput.getEpetraVector()));

  // RPP: Don't use the SetProblem - it breaks things if you redo the prec
  // Jacobian operator is set in the constructor
  //aztecSolverPtr->SetProblem(Problem);
  
  aztecSolverPtr->SetLHS(&(result.getEpetraVector()));
  aztecSolverPtr->SetRHS(&(nonConstInput.getEpetraVector()));

  // ************* Begin linear system scaling *******************
  if (scaling != 0) {
    //scaling->computeScaling(Problem);
    scaling->scaleLinearSystem(Problem);

    if (utils.isPrintProcessAndType(Utils::Details)) {
      cout << *scaling << endl;
    }
  }
  // ************* End linear system scaling *******************

  // Compute and set the Preconditioner in AztecOO if needed
  if (!isPrecConstructed && (precAlgorithm != None_)) {
    throwError("applyJacobianInverse", 
       "Preconditioner is not constructed!  Call createPreconditioner() first.");
  }

  // Get linear solver convergence parameters
  int maxit = p.getParameter("Max Iterations", 400);
  double tol = p.getParameter("Tolerance", 1.0e-6);
  bool reusePrec = p.getParameter("Reuse Preconditioner", false);
  
  if ( precAlgorithm == AztecOO_ )
    if ( !checkPreconditionerReuse() )
      aztecSolverPtr->SetAztecOption(AZ_pre_calc, AZ_reuse);

  int aztecStatus = -1;

  aztecStatus = aztecSolverPtr->Iterate(maxit, tol);
  
  // Unscale the linear system
  if (scaling != 0)
    scaling->unscaleLinearSystem(Problem);

  // Set the output parameters in the "Output" sublist
  if (outputSolveDetails) {
    NOX::Parameter::List& outputList = p.sublist("Output");
    int prevLinIters = 
      outputList.getParameter("Total Number of Linear Iterations", 0);
    int curLinIters = 0;
    double achievedTol = -1.0;
    curLinIters = aztecSolverPtr->NumIters();
    achievedTol = aztecSolverPtr->ScaledResidual();

    outputList.setParameter("Number of Linear Iterations", curLinIters);
    outputList.setParameter("Total Number of Linear Iterations", 
			    (prevLinIters + curLinIters));
    outputList.setParameter("Achieved Tolerance", achievedTol);
  }

  double endTime = timer.WallTime();
  timeApplyJacbianInverse += (endTime - startTime);

  if (aztecStatus != 0) 
    return false;
  
  return true;
}

//***********************************************************************
bool NOX::EpetraNew::LinearSystemAztecOO::
applyRightPreconditioning(bool useTranspose, 
			  Parameter::List& params,
			  const NOX::Epetra::Vector& input, 
			  NOX::Epetra::Vector& result) const
{
  int errorCode = 1;

  // Create the preconditioner if not already done.
  if (!isPrecConstructed) {
    throwError("applyRightPreconditioning", 
	 "Preconditioner is not constructed! Call createPreconditioner() first.");
  }

  if (precAlgorithm == None_) {
    if (&result != &input)
      result = input;
    return true;
  }
  else if (precAlgorithm == AztecOO_) {

    // RPP: We can not directly access Aztec preconditioners.
    // A cheesy way to apply an aztec preconditioner to an arbitrary 
    // vector is to call a Aztec.iterate() but only take one GMRES iteration.
    // This does NOT give the exact preconditioner, but it is a good
    // approximation.  We implement this here but highly recommend the 
    // use of IFPACK preconditioners if available!  

    // Zero out the temporary vector
    tmpVectorPtr->PutScalar(0.0);

    // Turn off printing in Aztec when using applyRightPreconditioner
    aztecSolverPtr->SetAztecOption(AZ_output,AZ_none);

    // Get the number of iterations in the preconditioner
    int numIters = params.getParameter("AztecOO Preconditioner Iterations", 1);
    
    AztecOO_Operator prec(aztecSolverPtr, numIters);
    
    errorCode = prec.ApplyInverse(input.getEpetraVector(), 
				  result.getEpetraVector());
  }
  else if (precAlgorithm == Ifpack_){

    if (useTranspose)
      ifpackPreconditionerPtr->SetUseTranspose(useTranspose);

    errorCode = ifpackPreconditionerPtr->ApplyInverse(input.getEpetraVector(), 
						   result.getEpetraVector());
    // Unset the transpose call
    if (useTranspose)
      ifpackPreconditionerPtr->SetUseTranspose(false);    

  }
  else if (precAlgorithm == UserDefined_) {

    if (useTranspose)
      precPtr->SetUseTranspose(true);

    errorCode = precPtr->ApplyInverse(input.getEpetraVector(), 
				      result.getEpetraVector());
    if (useTranspose)
      precPtr->SetUseTranspose(false);

  }
  else
    throwError("applyRightPreconditioning", 
	       "Parameter \"preconditioner\" is not vaild for this method");

  if (errorCode != 0) 
    return false;
  
  return true;
}

//***********************************************************************
bool NOX::EpetraNew::LinearSystemAztecOO::checkPreconditionerValidity() 
{
  if (precAlgorithm == None_)
    return true;

  else if ((precAlgorithm == AztecOO_) || (precAlgorithm == Ifpack_) ||
           (precAlgorithm == ML_)) {

    // Then we require an Epetra_RowMatrix
    if ((precType == EpetraRowMatrix) ||
	(precType == EpetraVbrMatrix) ||
	(precType == EpetraCrsMatrix)) {
      return true;
    }
    else {
      throwError("checkPreconditionerValidity", "Preconditioner must be an Epetra_RowMatrix derived object for AztecOO and Ifpack preconditioners");
    }

  }
  else if (precAlgorithm == UserDefined_){
    
    // Make sure a separate operator was supplied by the user
    if (precPtr == 0) {
      throwError("checkPreconditionerValidity", "Preconditioiner is NULL!");
    }
    return true;
  }
 
  return true;
}

//***********************************************************************
bool NOX::EpetraNew::LinearSystemAztecOO::
createPreconditioner(Epetra_Vector& x, Parameter::List& p, 
		     bool recomputeGraph) const
{
  double startTime = timer.WallTime();  

  if (precAlgorithm == None_) {
    return true;
  }

  // Apply Scaling
  Epetra_LinearProblem Problem(jacPtr, tmpVectorPtr, tmpVectorPtr);
  if (scaling != 0)
    scaling->scaleLinearSystem(Problem);

  if (utils.isPrintProcessAndType(Utils::LinearSolverDetails))
    cout << "\n       Computing a new precondtioner" << endl;;

  if (precAlgorithm == AztecOO_) {
    
    if (precMatrixSource == UseJacobian) {
      // The Jacobian has already been evaluated at the current solution.
      // Just set and enforce explicit constuction
      aztecSolverPtr->SetPrecMatrix(dynamic_cast<Epetra_RowMatrix*>(jacPtr));
      aztecSolverPtr->ConstructPreconditioner(conditionNumberEstimate);
    }
    else if (precMatrixSource == SeparateMatrix) {
      Epetra_RowMatrix& precMatrix = dynamic_cast<Epetra_RowMatrix&>(*precPtr);
      precInterfacePtr->computePreconditioner(x, &p);    
      aztecSolverPtr->SetPrecMatrix(&precMatrix);
      aztecSolverPtr->ConstructPreconditioner(conditionNumberEstimate);
    }

  }
  else if (precAlgorithm == Ifpack_) {
    
    if (precMatrixSource == UseJacobian) {
      createIfpackPreconditioner(p);
      aztecSolverPtr->SetPrecOperator(ifpackPreconditionerPtr);
    }
    else if (precMatrixSource == SeparateMatrix) {
      
      precInterfacePtr->computePreconditioner(x, &p);
      createIfpackPreconditioner(p);
      aztecSolverPtr->SetPrecOperator(ifpackPreconditionerPtr);
    }

  }
#ifdef HAVE_NOX_ML_EPETRA
  else if (precAlgorithm == ML_) {
    
    if (precMatrixSource == UseJacobian) {
      createMLPreconditioner(p);
      aztecSolverPtr->SetPrecOperator(MLPreconditionerPtr);
    }
    else if (precMatrixSource == SeparateMatrix) {
      
      precInterfacePtr->computePreconditioner(x, &p);
      createMLPreconditioner(p);
      aztecSolverPtr->SetPrecOperator(MLPreconditionerPtr);
    }

  }
#endif
  else if (precAlgorithm == UserDefined_) {

    precInterfacePtr->computePreconditioner(x, &p);
    aztecSolverPtr->SetPrecOperator(precPtr);

  }

  isPrecConstructed = true; 

  // Unscale the linear system
  if (scaling != 0)
    scaling->unscaleLinearSystem(Problem);

  double endTime = timer.WallTime();
  timeCreatePreconditioner += (endTime - startTime);

  return true;
}

//***********************************************************************
bool NOX::EpetraNew::LinearSystemAztecOO::
createIfpackPreconditioner(Parameter::List& p) const
{
  //for ifpack we need a VBR or CRS matrix to get the correct graph

  if (ifpackGraphPtr != 0) 
    throwError("createIfpackPreconditioner", "Ifpack Graph NOT NULL");
  if (ifpackPreconditionerPtr != 0) 
    throwError("createIfpackPreconditioner", "Ifpack Prec NOT NULL");

  //check to see if it is a VBR matrix
  if (precType == EpetraVbrMatrix) {

    Epetra_VbrMatrix* vbr = 0;

    if (precMatrixSource == UseJacobian)
      vbr = dynamic_cast<Epetra_VbrMatrix*>(jacPtr);
    else if (precMatrixSource == SeparateMatrix)
      vbr = dynamic_cast<Epetra_VbrMatrix*>(precPtr);

    if (vbr == 0)
      throwError("createIfpackPreconditioner", 
		 "Dynamic cast to VBR Matrix failed!");

    ifpackGraphPtr = new Ifpack_IlukGraph(vbr->Graph(),
					  p.getParameter("Fill Factor", 1),
					  p.getParameter("Overlap", 0));
    ifpackGraphPtr->ConstructFilledGraph();
    ifpackPreconditionerPtr = new Ifpack_CrsRiluk(*ifpackGraphPtr);
    ifpackPreconditionerPtr->InitValues(*vbr);
    ifpackPreconditionerPtr->Factor();
    return true;
  }

  // check to see if it is a Crs matrix
  if (precType == EpetraCrsMatrix) {

    Epetra_CrsMatrix* crs = 0;

    if (precMatrixSource == UseJacobian)
      crs = dynamic_cast<Epetra_CrsMatrix*>(jacPtr);
    else if (precMatrixSource == SeparateMatrix)
      crs = dynamic_cast<Epetra_CrsMatrix*>(precPtr);

    if (crs == 0)
      throwError("createIfpackPreconditioner", 
		 "Dynamic cast to CRS Matrix failed!");

    ifpackGraphPtr = new Ifpack_IlukGraph(crs->Graph(),
					  p.getParameter("Fill Factor", 1),
					  p.getParameter("Overlap", 0));
    ifpackGraphPtr->ConstructFilledGraph();
    ifpackPreconditionerPtr = new Ifpack_CrsRiluk(*ifpackGraphPtr);
    ifpackPreconditionerPtr->InitValues(*crs);
    ifpackPreconditionerPtr->Factor();
    return true;
  }
  
  // check to see if it is an operator that contains a Crs matrix
  if (precType == EpetraRowMatrix) {

    // The following checks only for FiniteDiffernce based operators and
    // further that the underlying matrix is in Crs format......
    Epetra_CrsMatrix* crs = 0;
    NOX::EpetraNew::FiniteDifference* FDoperator = 0;

    if (precMatrixSource == UseJacobian)
      FDoperator = dynamic_cast<NOX::EpetraNew::FiniteDifference*>(jacPtr);
    else if (precMatrixSource == SeparateMatrix)
      FDoperator = dynamic_cast<NOX::EpetraNew::FiniteDifference*>(precPtr);

    if (FDoperator != 0)
      crs = &(FDoperator->getUnderlyingMatrix());

    if (crs == 0)
      throwError("createIfpackPreconditioner", 
		 "FiniteDifference: Underlying matrix NOT CRS Matrix!");

    ifpackGraphPtr = new Ifpack_IlukGraph(crs->Graph(),
					  p.getParameter("Fill Factor", 1),
					  p.getParameter("Overlap", 0));
    ifpackGraphPtr->ConstructFilledGraph();
    ifpackPreconditionerPtr = new Ifpack_CrsRiluk(*ifpackGraphPtr);
    ifpackPreconditionerPtr->InitValues(*crs);
    ifpackPreconditionerPtr->Factor();
    return true;
  }
  
  // If we made it this far, this routine should not have been called
  // in the first place.  An incorrect prec matrix object was supplied and
  // should have been caught in the checkOperatorValidity() method.
  throwError("createIfpackPreconditioner",
	     "Ifpack requires a CRS or VBR matrix!");

  return false;
}

#ifdef HAVE_NOX_ML_EPETRA
//***********************************************************************
bool NOX::EpetraNew::LinearSystemAztecOO::
createMLPreconditioner(Parameter::List& p) const
{
  //for ML we need a Epetra_RowMatrix

  if (MLPreconditionerPtr != 0) 
    throwError("createMLPreconditioner", "ML Prec NOT NULL");

  // Ensure we have a valid Teuchos parameter list to pass to ML
  const Teuchos::ParameterList* constTeuchosParams = 0;
  if (p.isParameter("ML Teuchos Parameter List")) {
    if (p.isParameterArbitrary("ML Teuchos Parameter List")) {
      constTeuchosParams = dynamic_cast<TeuchosParamsAsNOXArbitrary*>
        (p.getArbitraryParameter("ML Teuchos Parameter List").clone())->getTeuchosParams();
      if ( !constTeuchosParams )
        if (utils.isPrintProcess())
          cout << "ERROR: NOX::EpetraNew::LinearSystemAztecOO::"
               << "createMLPreconditioner() - "
               << "Could not cast parameter in \"ML Teuchos Parameter List\""
               << " to type Teuchos::ParameterList!\n" << endl;
    }
    else {
      if (utils.isPrintProcess())
        cout << "ERROR: NOX::EpetraNew::LinearSystemAztecOO::"
             << "createMLPreconditioner() - "
	     << "The parameter in \"ML Teuchos Parameter List\" "
	     << "must be derived from an arbitrary parameter!" << endl;
    }
  }
  if ( !constTeuchosParams ) {
    if (utils.isPrintProcess())
      cout << "ERROR: NOX::EpetraNew::LinearSystemAztecOO::"
           << "createMLPreconditioner() - "
           << "Could not obtain the required Teuchos::ParameterList "
           << "from \"ML Teuchos Parameter List\" " << endl;
    throw "NOX Error";
  }

  // check to see if it is a Crs matrix
  if (precType == EpetraCrsMatrix) {

    Epetra_CrsMatrix* crs = 0;

    if (precMatrixSource == UseJacobian)
      crs = dynamic_cast<Epetra_CrsMatrix*>(jacPtr);
    else if (precMatrixSource == SeparateMatrix)
      crs = dynamic_cast<Epetra_CrsMatrix*>(precPtr);

    if (crs == 0)
      throwError("createMLPreconditioner", 
		 "Dynamic cast to CRS Matrix failed!");

    MLPreconditionerPtr = new ML_Epetra::MultiLevelPreconditioner(
                                  *crs, *constTeuchosParams, true);
    return true;
  }

  // check to see if it is an operator that contains a Crs matrix
  if (precType == EpetraRowMatrix) {

    // The following checks only for FiniteDiffernce based operators and
    // further that the underlying matrix is in Crs format......
    Epetra_CrsMatrix* crs = 0;
    NOX::EpetraNew::FiniteDifference* FDoperator = 0;

    if (precMatrixSource == UseJacobian)
      FDoperator = dynamic_cast<NOX::EpetraNew::FiniteDifference*>(jacPtr);
    else if (precMatrixSource == SeparateMatrix)
      FDoperator = dynamic_cast<NOX::EpetraNew::FiniteDifference*>(precPtr);

    if (FDoperator != 0)
      crs = &(FDoperator->getUnderlyingMatrix());

    if (crs == 0)
      throwError("createMLPreconditioner", 
		 "FiniteDifference: Underlying matrix NOT CRS Matrix!");

    MLPreconditionerPtr = new ML_Epetra::MultiLevelPreconditioner(
                                  *crs, *constTeuchosParams, true);
    return true;
  }
  
  // If we made it this far, this routine should not have been called
  // in the first place.  An incorrect prec matrix object was supplied and
  // should have been caught in the checkOperatorValidity() method.
  throwError("createMLPreconditioner",
	     "ML requires an Epetra_RowMatrix !");

  return false;
}
#endif

//***********************************************************************
bool NOX::EpetraNew::LinearSystemAztecOO::destroyPreconditioner() const
{ 
  if (isPrecConstructed) {
    // For Aztec preconditioners. 
    if (precAlgorithm == AztecOO_) {
      aztecSolverPtr->DestroyPreconditioner();
    }
    // For IFPACK Preconditioners
    else if (precAlgorithm == Ifpack_) {
      delete ifpackPreconditionerPtr;
      ifpackPreconditionerPtr = 0;
      delete ifpackGraphPtr;
      ifpackGraphPtr = 0;
    }
#ifdef HAVE_NOX_ML_EPETRA
    else if (precAlgorithm == ML_) {
      delete MLPreconditionerPtr;
      MLPreconditionerPtr = 0;
    }
#endif
    if (utils.isPrintProcessAndType(Utils::LinearSolverDetails)) {
      cout << "\n       Destroying preconditioner" << endl;
    }
  }
  isPrecConstructed = false;
  return true;
}

//***********************************************************************
NOX::EpetraNew::LinearSystemAztecOO::OperatorType 
NOX::EpetraNew::LinearSystemAztecOO::getOperatorType(const Epetra_Operator& Op)
{
  //***************
  //*** NOTE: The order in which the following tests occur is important!
  //***************

  const Epetra_Operator* testOperator = 0;
  
  // Is it an Epetra_CrsMatrix ?
  testOperator = dynamic_cast<const Epetra_CrsMatrix*>(&Op);
  if (testOperator != 0) 
    return EpetraCrsMatrix; 

  // Is it an Epetra_VbrMatrix ?
  testOperator = dynamic_cast<const Epetra_VbrMatrix*>(&Op);
  if (testOperator != 0) 
    return EpetraVbrMatrix; 

  // Is it an Epetra_RowMatrix ?
  testOperator = dynamic_cast<const Epetra_RowMatrix*>(&Op);
  if (testOperator != 0) 
    return EpetraRowMatrix; 

  // Otherwise it must be an Epetra_Operator!
  return EpetraOperator;
}

//***********************************************************************
bool NOX::EpetraNew::LinearSystemAztecOO::computeJacobian(Epetra_Vector& x)
{
  bool success = jacInterfacePtr->computeJacobian(x);
  return success;
}

//***********************************************************************
void NOX::EpetraNew::LinearSystemAztecOO::
resetScaling(NOX::EpetraNew::Scaling& scalingObject)
{
  scaling = &scalingObject;
  return;
}

//***********************************************************************
void NOX::EpetraNew::LinearSystemAztecOO::
throwError(const string& functionName, const string& errorMsg) const
{
  if (utils.isPrintProcessAndType(Utils::Error)) {
    cout << "NOX::EpetraNew::LinearSystemAztecOO::" << functionName 
	 << " - " << errorMsg << endl;
  }
  throw "NOX Error";
}

//***********************************************************************
const NOX::EpetraNew::Interface::Jacobian& 
NOX::EpetraNew::LinearSystemAztecOO::getJacobianInterface() const
{
  return *jacInterfacePtr;
}

//***********************************************************************
const NOX::EpetraNew::Interface::Preconditioner& 
NOX::EpetraNew::LinearSystemAztecOO::getPrecInterface() const
{
  return *precInterfacePtr;
}

//***********************************************************************
const Epetra_Operator& 
NOX::EpetraNew::LinearSystemAztecOO::getJacobianOperator() const
{
  return *jacPtr;
}

//***********************************************************************
const Epetra_Operator& 
NOX::EpetraNew::LinearSystemAztecOO::getPrecOperator() const
{
  return *precPtr;
}

//***********************************************************************
const Epetra_Operator&
NOX::EpetraNew::LinearSystemAztecOO::getGeneratedPrecOperator() const
{
  return *(aztecSolverPtr->GetPrecOperator());
}

//***********************************************************************
bool NOX::EpetraNew::LinearSystemAztecOO::checkPreconditionerReuse()
{
  if (!isPrecConstructed)
    return false;

  precQueryCounter++;

  // This allows reuse for the entire nonlinear solve
  if( maxAgeOfPrec == -2 )
    return true;
  
  // This allows one recompute of the preconditioner followed by reuse 
  // for the remainder of the nonlinear solve
  else if( maxAgeOfPrec == -1 ) {
    maxAgeOfPrec = -2;
    return false;
  }
  
  // This is the typical use 
  else
    if( precQueryCounter == maxAgeOfPrec ) {
      precQueryCounter = 0;
      return false;
    }
    else
      return true;
}

//***********************************************************************
double 
NOX::EpetraNew::LinearSystemAztecOO::getTimeCreatePreconditioner() const
{
  return timeCreatePreconditioner;
}

//***********************************************************************
double 
NOX::EpetraNew::LinearSystemAztecOO::getTimeApplyJacobianInverse() const
{
  return timeApplyJacbianInverse;
}

void
NOX::EpetraNew::LinearSystemAztecOO::setJacobianOperatorForSolve(
					 const Epetra_Operator& solveJacOp)
{
  solveJacOpPtr = const_cast<Epetra_Operator*>(&solveJacOp);
  OperatorType solveOpType = getOperatorType(solveJacOp);
  if ((solveOpType == EpetraRowMatrix) ||
      (solveOpType == EpetraVbrMatrix) ||
      (solveOpType == EpetraCrsMatrix)) {
    aztecSolverPtr->SetUserMatrix(dynamic_cast<Epetra_RowMatrix*>(solveJacOpPtr));
  }
  else
    aztecSolverPtr->SetUserOperator(solveJacOpPtr);
}
