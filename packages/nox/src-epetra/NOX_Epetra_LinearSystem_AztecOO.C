
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

#include "NOX_Epetra_LinearSystem_AztecOO.H"	// class definition

// NOX includes
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "NOX_Epetra_MatrixFree.H"
#include "NOX_Epetra_FiniteDifference.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Epetra_Scaling.H"
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
#include "Ifpack.h"
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"
#include "Teuchos_ParameterList.hpp"

// EpetraExt includes for dumping a matrix
// #ifdef HAVE_NOX_DEBUG
#ifdef HAVE_NOX_EPETRAEXT
#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_RowMatrixOut.h"
#endif
// #endif

#ifdef HAVE_NOX_ML_EPETRA
#include "Teuchos_ParameterList.hpp"
#endif

#include <typeinfo>

//***********************************************************************
NOX::Epetra::LinearSystemAztecOO::
LinearSystemAztecOO(
 Teuchos::ParameterList& printParams, 
 Teuchos::ParameterList& linearSolverParams, 
 const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
 const NOX::Epetra::Vector& cloneVector,
 const Teuchos::RCP<NOX::Epetra::Scaling> s):
  utils(printParams),
  jacType(EpetraOperator),
  precType(EpetraOperator),
  precMatrixSource(UseJacobian),
  scaling(s),
  conditionNumberEstimate(0.0),
  isPrecConstructed(false),
  precQueryCounter(0),
  maxAgeOfPrec(1),
  timer(cloneVector.getEpetraVector().Comm()),
  timeCreatePreconditioner(0.0),
  timeApplyJacbianInverse(0.0)
{
  // Allocate solver
  aztecSolverPtr = Teuchos::rcp(new AztecOO());
  tmpVectorPtr = Teuchos::rcp(new NOX::Epetra::Vector(cloneVector));

  // Neither Jacobian or Preconditioner are supplied
  createJacobianOperator(printParams, linearSolverParams, iReq, cloneVector);
  createPrecOperator(printParams, linearSolverParams, iReq, cloneVector);

  reset(linearSolverParams);
}

//***********************************************************************
NOX::Epetra::LinearSystemAztecOO::
LinearSystemAztecOO(
 Teuchos::ParameterList& printParams, 
 Teuchos::ParameterList& linearSolverParams,  
 const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
 const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac, 
 const Teuchos::RCP<Epetra_Operator>& jacobian,
 const NOX::Epetra::Vector& cloneVector,
 const Teuchos::RCP<NOX::Epetra::Scaling> s):
  utils(printParams),
  jacInterfacePtr(iJac),
  jacType(EpetraOperator),
  jacPtr(jacobian),
  precType(EpetraOperator),
  precMatrixSource(UseJacobian),
  scaling(s),
  conditionNumberEstimate(0.0),
  isPrecConstructed(false),
  precQueryCounter(0),
  maxAgeOfPrec(1),
  timer(cloneVector.getEpetraVector().Comm()),
  timeCreatePreconditioner(0.0),
  timeApplyJacbianInverse(0.0)
{
  // Allocate solver
  aztecSolverPtr = Teuchos::rcp(new AztecOO());
  tmpVectorPtr = Teuchos::rcp(new NOX::Epetra::Vector(cloneVector));

  // Jacobian operator is supplied 
  jacType = getOperatorType(*jacPtr);
  // Preconditioner is not supplied
  createPrecOperator(printParams, linearSolverParams, iReq, cloneVector);

  reset(linearSolverParams);
}

//***********************************************************************
NOX::Epetra::LinearSystemAztecOO::
LinearSystemAztecOO(
 Teuchos::ParameterList& printParams, 
 Teuchos::ParameterList& linearSolverParams, 
 const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
 const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec, 
 const Teuchos::RCP<Epetra_Operator>& preconditioner,
 const NOX::Epetra::Vector& cloneVector,
 const Teuchos::RCP<NOX::Epetra::Scaling> s):
  utils(printParams),
  jacType(EpetraOperator),
  precInterfacePtr(iPrec),
  precType(EpetraOperator),
  precPtr(preconditioner),
  precMatrixSource(SeparateMatrix),
  scaling(s),
  conditionNumberEstimate(0.0),
  isPrecConstructed(false),
  precQueryCounter(0),
  maxAgeOfPrec(1),
  timer(cloneVector.getEpetraVector().Comm()),
  timeCreatePreconditioner(0.0),
  timeApplyJacbianInverse(0.0)
{
  // Allocate solver
  aztecSolverPtr = Teuchos::rcp(new AztecOO());
  tmpVectorPtr = Teuchos::rcp(new NOX::Epetra::Vector(cloneVector));

  // Jacobian operator is not supplied
  createJacobianOperator(printParams, linearSolverParams, iReq, cloneVector);
  // Preconditioner operator is supplied
  precType = getOperatorType(*precPtr);

  reset(linearSolverParams);
}

//***********************************************************************
NOX::Epetra::LinearSystemAztecOO::
LinearSystemAztecOO(
 Teuchos::ParameterList& printParams, 
 Teuchos::ParameterList& linearSolverParams,
 const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac, 
 const Teuchos::RCP<Epetra_Operator>& jacobian,
 const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec, 
 const Teuchos::RCP<Epetra_Operator>& preconditioner,
 const NOX::Epetra::Vector& cloneVector,
 const Teuchos::RCP<NOX::Epetra::Scaling> s):
  utils(printParams),
  jacInterfacePtr(iJac),
  jacType(EpetraOperator),
  jacPtr(jacobian),
  precInterfacePtr(iPrec),
  precType(EpetraOperator),
  precPtr(preconditioner),
  precMatrixSource(SeparateMatrix),
  scaling(s),
  conditionNumberEstimate(0.0),
  isPrecConstructed(false),
  precQueryCounter(0),
  maxAgeOfPrec(1),
  timer(cloneVector.getEpetraVector().Comm()),
  timeCreatePreconditioner(0.0),
  timeApplyJacbianInverse(0.0)
{
  // Allocate solver
  aztecSolverPtr = Teuchos::rcp(new AztecOO());
  tmpVectorPtr = Teuchos::rcp(new NOX::Epetra::Vector(cloneVector));

  // Both operators are supplied
  jacType = getOperatorType(*jacPtr);
  precType = getOperatorType(*precPtr);

  reset(linearSolverParams);
}

//***********************************************************************
NOX::Epetra::LinearSystemAztecOO::~LinearSystemAztecOO() 
{
  destroyPreconditioner();
}

//***********************************************************************
void NOX::Epetra::LinearSystemAztecOO::
reset(Teuchos::ParameterList& linearSolverParams)
{

  // First remove any preconditioner that may still be active
  destroyPreconditioner();

  // Set the requested preconditioning.
  std::string prec = linearSolverParams.get("Preconditioner", "None");
  if (prec == "AztecOO")
    precAlgorithm = AztecOO_;
  else if (prec == "Ifpack")
    precAlgorithm = Ifpack_;
  else if (prec == "New Ifpack")
    precAlgorithm = NewIfpack_;
#ifdef HAVE_NOX_ML_EPETRA
  else if (prec == "ML")
    precAlgorithm = ML_;
#endif
  else if (prec == "User Defined")
    precAlgorithm = UserDefined_;
  else if (prec == "None") 
    precAlgorithm = None_;
  else {
    std::string errorMessage = "Option for \"Preconditioner\" is invalid!";
    throwError("reset()", errorMessage);
  }
    
  // Make sure the correct objects were supplied for the requested
  // preconditioning choice.
  checkPreconditionerValidity();

  zeroInitialGuess = 
    linearSolverParams.get("Zero Initial Guess", false);

  manualScaling = 
    linearSolverParams.get("Compute Scaling Manually", true);

  // Place linear solver details in the "Output" sublist of the
  // "Linear Solver" parameter list
  outputSolveDetails = 
    linearSolverParams.get("Output Solver Details", true);

  throwErrorOnPrecFailure = 
    linearSolverParams.get("Throw Error on Prec Failure", true);

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
//     aztecSolverPtr->SetUserMatrix(dynamic_cast<Epetra_RowMatrix*>(jacPtr.get()));
//   }
//   else
//     aztecSolverPtr->SetUserOperator(jacPtr.get());
  
  // Set the major aztec options.  Must be called after the first 
  // SetProblem() call.
  setAztecOptions(linearSolverParams, *aztecSolverPtr);

  // Setup the preconditioner reuse policy
  std::string preReusePolicyName = 
    linearSolverParams.get("Preconditioner Reuse Policy", "Rebuild");
  if (preReusePolicyName == "Rebuild")
    precReusePolicy = PRPT_REBUILD;
  else if (preReusePolicyName == "Recompute")
    precReusePolicy = PRPT_RECOMPUTE;
  else if (preReusePolicyName == "Reuse")
    precReusePolicy = PRPT_REUSE;
  else {
    std::string errorMessage = "Option for \"Preconditioner Reuse Policy\" is invalid! \nPossible options are \"Reuse\", \"Rebuild\", and \"Recompute\".";
    throwError("reset()", errorMessage);
  }
  maxAgeOfPrec = linearSolverParams.get("Max Age Of Prec", 1);
  precQueryCounter = 0;

// #ifdef HAVE_NOX_DEBUG
#ifdef HAVE_NOX_EPETRAEXT
  linearSolveCount = 0;
#endif
// #endif

}

//***********************************************************************
void NOX::Epetra::LinearSystemAztecOO::
setAztecOptions(Teuchos::ParameterList& p, AztecOO& aztec) const
{
  // Set the output and error streams in the aztec solver
  aztec.SetOutputStream(utils.pout(NOX::Utils::LinearSolverDetails));
  aztec.SetErrorStream(utils.perr());

  // Set the Aztec Solver
  std::string linearSolver = p.get("Aztec Solver", "GMRES");
  if (linearSolver == "CG")
    aztec.SetAztecOption(AZ_solver, AZ_cg);
  else if (linearSolver == "GMRES")
    aztec.SetAztecOption(AZ_solver, AZ_gmres);
  else if (linearSolver == "GMRESR")
    aztec.SetAztecOption(AZ_solver, AZ_GMRESR);
  else if (linearSolver == "CGS")
    aztec.SetAztecOption(AZ_solver, AZ_cgs);
  else if (linearSolver == "TFQMR")
    aztec.SetAztecOption(AZ_solver, AZ_tfqmr);
  else if (linearSolver == "BiCGStab")
    aztec.SetAztecOption(AZ_solver, AZ_bicgstab);
  else if (linearSolver == "LU")
    aztec.SetAztecOption(AZ_solver, AZ_lu);
  else {
    utils.out() << "ERROR: NOX::Epetra::Group::setAztecOptions" << std::endl
	 << "\"Aztec Solver\" parameter \"" << linearSolver 
	 <<  "\" is invalid!" << std::endl;
    throw "NOX Error";
  }
 
  // Preconditioning where AztecOO inverts the Preconditioning Matrix
  if (precAlgorithm == AztecOO_) {
    
    std::string aztecPreconditioner = p.get("Aztec Preconditioner", "ilu");

    if (aztecPreconditioner == "ilu") {
      aztec.SetAztecOption(AZ_precond, AZ_dom_decomp);
      aztec.SetAztecOption(AZ_overlap, p.get("Overlap", 0));
      aztec.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
      aztec.SetAztecOption(AZ_graph_fill, p.get("Graph Fill", 0));
    }
    else if (aztecPreconditioner == "ilut") { 
      aztec.SetAztecOption(AZ_precond, AZ_dom_decomp);
      aztec.SetAztecOption(AZ_overlap, p.get("Overlap", 0));
      aztec.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
      aztec.SetAztecParam(AZ_drop, p.get("Drop Tolerance", 0.0));
      aztec.SetAztecParam(AZ_ilut_fill, p.get("Fill Factor", 1.0));
    }
    else if (aztecPreconditioner == "Jacobi") {
      aztec.SetAztecOption(AZ_precond, AZ_Jacobi);
      aztec.SetAztecOption(AZ_poly_ord, p.get("Steps", 3));
    }
    else if (aztecPreconditioner == "Symmetric Gauss-Seidel") {
      aztec.SetAztecOption(AZ_precond, AZ_sym_GS);
      aztec.SetAztecOption(AZ_poly_ord, p.get("Steps", 3));
    }
    else if (aztecPreconditioner == "Polynomial") {
      aztec.SetAztecOption(AZ_precond, AZ_Neumann);
      aztec.SetAztecOption(AZ_poly_ord, p.get("Polynomial Order", 3));
    }
    else if (aztecPreconditioner == "Least-squares Polynomial") {
      aztec.SetAztecOption(AZ_precond, AZ_ls);
      aztec.SetAztecOption(AZ_poly_ord, p.get("Polynomial Order", 3));
    }
    else {
      std::string errorMessage = "\"Aztec Preconditioner\" parameter is invalid!";
      throwError("setAztecOptions", errorMessage);
    }

  }
  else 
    aztec.SetAztecOption(AZ_precond, AZ_none);
    
  // Turn on RCM reordering in conjunction with domain decomp preconditioning
  // default is "Disabled" = no reordering
  std::string rcmReordering = p.get("RCM Reordering", "Disabled");
  if (rcmReordering == "Enabled")
    aztec.SetAztecOption(AZ_reorder, 1);
  else if (rcmReordering == "Disabled")
    aztec.SetAztecOption(AZ_reorder, 0);
  else {
    std::string errorMessage = "\"RCM Reordering\" parameter is invalid!";
    throwError("setAztecOptions", errorMessage);
  }
    
  // Gram-Schmidt orthogonalization procedure
  std::string orthog = p.get("Orthogonalization", "Classical");
  if (orthog == "Classical") 
    aztec.SetAztecOption(AZ_orthog, AZ_classic);
  else if (orthog == "Modified")
    aztec.SetAztecOption(AZ_orthog, AZ_modified);
  else {
    std::string errorMessage = "\"Orthogonalization\" parameter is invalid!";
    throwError("setAztecOptions()", errorMessage);
  }

  // Size of the krylov subspace
  aztec.SetAztecOption(AZ_kspace, 
		       p.get("Size of Krylov Subspace", 300));

  // Convergence criteria to use in the linear solver
  std::string convCriteria = p.get("Convergence Test", "r0");
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
    std::string errorMessage = "\"Convergence Test\" parameter is invalid!";
    throwError("setAztecOptions()", errorMessage);
  }

  // Set the ill-conditioning threshold for the upper hessenberg matrix
  if (p.isParameter("Ill-Conditioning Threshold")) {
    aztec.SetAztecParam(AZ_ill_cond_thresh, 
			p.get("Ill-Conditioning Threshold", 1.0e+11));
  }

  // Frequency of linear solve residual output
  if (utils.isPrintType(Utils::LinearSolverDetails))
    aztec.SetAztecOption(AZ_output, p.get("Output Frequency", 
						   AZ_last));
  else
    aztec.SetAztecOption(AZ_output, p.get("Output Frequency", 0));

  // Print a summary of the aztec options if "Details" is enabled
  if (utils.isPrintType(Utils::Debug)) {
    //aztec.CheckInput();
  }

  // Some Debugging utilities 
//   if (utils.isPrintType(Utils::Debug)) {
//     utils.out() << "NOX::Epetra::LinearSystemAztecOO Operator Information" << std::endl;
//     utils.out() << "jacType = " << jacType << std::endl;
//     utils.out() << "jacPtr = " << jacPtr << std::endl;
//     utils.out() << "jacInterfacePtr = " << jacInterfacePtr << std::endl;
//     utils.out() << "precType = " << precType << std::endl;
//     utils.out() << "precPtr = " << precPtr << std::endl;
//     utils.out() << "precInterfacePtr = " << precInterfacePtr << std::endl;
//   }
    
  return;
}

//***********************************************************************
bool NOX::Epetra::LinearSystemAztecOO::createJacobianOperator(
       Teuchos::ParameterList& printParams,
       Teuchos::ParameterList& lsParams,
       const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
       const NOX::Epetra::Vector& cloneVector)
{
  std::string choice = lsParams.get("Jacobian Operator", "Matrix-Free");

  if (choice == "Matrix-Free") {
    jacPtr = 
      Teuchos::rcp(new MatrixFree(printParams, iReq, cloneVector));
    jacInterfacePtr = 
      Teuchos::rcp_dynamic_cast<NOX::Epetra::Interface::Jacobian>(jacPtr);
    jacType = EpetraOperator;
  }
  else if (choice == "Finite Difference") {
    jacPtr = 
      Teuchos::rcp(new FiniteDifference(printParams, iReq, cloneVector));
    jacInterfacePtr = 
      Teuchos::rcp_dynamic_cast<NOX::Epetra::Interface::Jacobian>(jacPtr);
    jacType = EpetraRowMatrix;
  }
  else    
    throwError("createJacobianOperator", 
       "The specified value for parameter \" Jacobian Operator\" is not valid");

  return true;
}

//***********************************************************************
bool NOX::Epetra::LinearSystemAztecOO::createPrecOperator(
       Teuchos::ParameterList& printParams,
       Teuchos::ParameterList& lsParams,
       const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
       const NOX::Epetra::Vector& cloneVector)
{
  std::string choice = lsParams.get("Preconditioner Operator", 
					"Use Jacobian");

  if (choice == "Use Jacobian") {
    precType = jacType;
    precMatrixSource = UseJacobian;
  }
  else if (choice == "Finite Difference") {
    precPtr = 
      Teuchos::rcp(new FiniteDifference(printParams, iReq, cloneVector));
    precInterfacePtr = Teuchos::rcp_dynamic_cast
      <NOX::Epetra::Interface::Preconditioner>(precPtr);
    precType = EpetraRowMatrix;
    precMatrixSource = SeparateMatrix;
  }
  else    
    throwError("createPreconditionerOperator", 
       "The value for the parameter \" Preconditioner Operator\" is not valid");

  return true;
}

//***********************************************************************
bool NOX::Epetra::LinearSystemAztecOO::
applyJacobian(const NOX::Epetra::Vector& input, 
	      NOX::Epetra::Vector& result) const
{
  jacPtr->SetUseTranspose(false);
  int status = jacPtr->Apply(input.getEpetraVector(), 
				  result.getEpetraVector());
  return (status == 0);
}

//***********************************************************************
bool NOX::Epetra::LinearSystemAztecOO::
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
bool NOX::Epetra::LinearSystemAztecOO::
applyJacobianInverse(Teuchos::ParameterList &p,
		     const NOX::Epetra::Vector& input, 
		     NOX::Epetra::Vector& result)
{
  // AGS: Rare option, similar to Max Iters=1 but twice as fast.
    if ( p.get("Use Preconditioner as Solver", false) ) 
      return applyRightPreconditioning(false, p, input, result);


  double startTime = timer.WallTime();

  // Need non-const version of the input vector
  // Epetra_LinearProblem requires non-const versions so we can perform
  // scaling of the linear problem.
  NOX::Epetra::Vector& nonConstInput = const_cast<NOX::Epetra::Vector&>(input);
  
  // Zero out the delta X of the linear problem if requested by user.
  if (zeroInitialGuess)
    result.init(0.0);

  // Create Epetra linear problem object for the linear solve
  Epetra_LinearProblem Problem(jacPtr.get(), 
  			       &(result.getEpetraVector()), 
			       &(nonConstInput.getEpetraVector()));

  // Set objects in aztec solver.
  // RPP: Don't use "aztecSolverPtr->SetProblem(Problem);", it breaks
  // things if you rebuild the prec.  Also, you MUST set Jac operator
  // before Prec Operator in AztecOO object.
  this->setAztecOOJacobian();
  this->setAztecOOPreconditioner();
  aztecSolverPtr->SetLHS(&(result.getEpetraVector()));
  aztecSolverPtr->SetRHS(&(nonConstInput.getEpetraVector()));

  // ************* Begin linear system scaling *******************
  if ( !Teuchos::is_null(scaling) ) {

    if ( !manualScaling )
      scaling->computeScaling(Problem);
    
    scaling->scaleLinearSystem(Problem);

    if (utils.isPrintType(Utils::Details)) {
      utils.out() << *scaling << std::endl;
    }
  }
  // ************* End linear system scaling *******************

  // Use EpetraExt to dump linear system if debuggging
// #ifdef HAVE_NOX_DEBUG
#ifdef HAVE_NOX_EPETRAEXT

  ++linearSolveCount;
  std::ostringstream iterationNumber;
  iterationNumber << linearSolveCount;
    
  std::string prefixName = p.get("Write Linear System File Prefix", 
				 "NOX_LinSys");
  std::string postfixName = iterationNumber.str();
  postfixName += ".mm";

  if (p.get("Write Linear System", false)) {

    std::string mapFileName = prefixName + "_Map_" + postfixName;
    std::string jacFileName = prefixName + "_Jacobian_" + postfixName;    
    std::string rhsFileName = prefixName + "_RHS_" + postfixName;
    
    Epetra_RowMatrix* printMatrix = NULL;
    printMatrix = dynamic_cast<Epetra_RowMatrix*>(jacPtr.get()); 

    if (printMatrix == NULL) {
      std::cout << "Error: NOX::Epetra::LinearSystemAztecOO::applyJacobianInverse() - "
	   << "Could not cast the Jacobian operator to an Epetra_RowMatrix!"
	   << "Please set the \"Write Linear System\" parameter to false."
	   << std::endl;
      throw "NOX Error";
    }

    EpetraExt::BlockMapToMatrixMarketFile(mapFileName.c_str(), 
					  printMatrix->RowMatrixRowMap()); 
    EpetraExt::RowMatrixToMatrixMarketFile(jacFileName.c_str(), *printMatrix, 
					   "test matrix", "Jacobian XXX");
    EpetraExt::MultiVectorToMatrixMarketFile(rhsFileName.c_str(), 
					     nonConstInput.getEpetraVector());

  }
#endif
// #endif

  // Make sure preconditioner was constructed if requested
  if (!isPrecConstructed && (precAlgorithm != None_)) {
    throwError("applyJacobianInverse", 
       "Preconditioner is not constructed!  Call createPreconditioner() first.");
  }

  // Get linear solver convergence parameters
  int maxit = p.get("Max Iterations", 400);
  double tol = p.get("Tolerance", 1.0e-6);
  
  int aztecStatus = -1;

  // RPP: This is a hack to get aztec to reuse the preconditioner.
  // There is a bug in AztecOO in how it stores old
  // preconditioners. The storage bin is based on the aztec options
  // list.  When we call ConstructPreconditioner, the option for
  // az_pre_calc is set to AZ_calc, so the preconditioner is stored
  // with that option.  But then the routine toggles the AZ_pre_calc
  // flag to AZ_reuse.  When we call iterate, the first solve works,
  // but subsequent solves fail to find the preconditioner.
  // Will try get Alan to fix for release 7.0.
  if (precAlgorithm == AztecOO_ && 
      precReusePolicy == PRPT_REUSE) 
    aztecSolverPtr->SetAztecOption(AZ_pre_calc, AZ_calc);

  aztecStatus = aztecSolverPtr->Iterate(maxit, tol);

  // Unscale the linear system
  if ( !Teuchos::is_null(scaling) )
    scaling->unscaleLinearSystem(Problem);

  // Set the output parameters in the "Output" sublist
  if (outputSolveDetails) {
    Teuchos::ParameterList& outputList = p.sublist("Output");
    int prevLinIters = 
      outputList.get("Total Number of Linear Iterations", 0);
    int curLinIters = 0;
    double achievedTol = -1.0;
    curLinIters = aztecSolverPtr->NumIters();
    achievedTol = aztecSolverPtr->ScaledResidual();

    outputList.set("Number of Linear Iterations", curLinIters);
    outputList.set("Total Number of Linear Iterations", 
			    (prevLinIters + curLinIters));
    outputList.set("Achieved Tolerance", achievedTol);
  }

  // Dump solution of linear system
// #ifdef HAVE_NOX_DEBUG
#ifdef HAVE_NOX_EPETRAEXT
  if (p.get("Write Linear System", false)) {
    std::string lhsFileName = prefixName + "_LHS_" + postfixName;
    EpetraExt::MultiVectorToMatrixMarketFile(lhsFileName.c_str(), 
					   result.getEpetraVector());
  }
#endif
// #endif

  double endTime = timer.WallTime();
  timeApplyJacbianInverse += (endTime - startTime);

  if (aztecStatus != 0) 
    return false;
  
  return true;
}

//***********************************************************************
bool NOX::Epetra::LinearSystemAztecOO::
applyRightPreconditioning(bool useTranspose, 
			  Teuchos::ParameterList& params,
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
    tmpVectorPtr->init(0.0);

    // Turn off printing in Aztec when using applyRightPreconditioner
    aztecSolverPtr->SetAztecOption(AZ_output,AZ_none);

    // Get the number of iterations in the preconditioner
    int numIters = params.get("AztecOO Preconditioner Iterations", 1);
    
    AztecOO_Operator prec(aztecSolverPtr.get(), numIters);
    
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
  else if (precAlgorithm == NewIfpack_){

    if (useTranspose)
      newIfpackPreconditionerPtr->SetUseTranspose(useTranspose);

    errorCode = newIfpackPreconditionerPtr->ApplyInverse(input.getEpetraVector(), 
						   result.getEpetraVector());
    // Unset the transpose call
    if (useTranspose)
      newIfpackPreconditionerPtr->SetUseTranspose(false);    

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

  if (errorCode != 0) {
    std::string msg = "Error - NOX::Epetra::LinearSystemAztecOO::applyRightPreconditioning() - A non-zero error code has been returned from the preconditioner.";
    if (throwErrorOnPrecFailure) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
    }
    else {
      if (utils.isPrintType(NOX::Utils::Warning))
	utils.out() << msg << std::endl;
    }
    return false;
  }
  
  return true;
}

//***********************************************************************
bool NOX::Epetra::LinearSystemAztecOO::checkPreconditionerValidity() 
{
  if (precAlgorithm == None_)
    return true;


  else if ((precAlgorithm == AztecOO_)   || (precAlgorithm == Ifpack_) ||
           (precAlgorithm == NewIfpack_) || (precAlgorithm == ML_)) {

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
    if (Teuchos::is_null(precPtr)) {
      throwError("checkPreconditionerValidity", "Preconditioiner is NULL!");
    }
    return true;
  }
 
  return true;
}

//***********************************************************************
bool NOX::Epetra::LinearSystemAztecOO::
createPreconditioner(const NOX::Epetra::Vector& x, Teuchos::ParameterList& p, 
		     bool recomputeGraph) const
{
  double startTime = timer.WallTime();  

  if (utils.isPrintType(Utils::LinearSolverDetails))
    utils.out() << "\n       Creating a new preconditioner" << std::endl;;

  if (precAlgorithm == None_) {
    return true;
  }

  // If Scaling exists, scale the preconditioner
  Epetra_LinearProblem Problem(jacPtr.get(), 
			       &(tmpVectorPtr->getEpetraVector()), 
			       &(tmpVectorPtr->getEpetraVector()));
  if ( !Teuchos::is_null(scaling) ) {
    if (!manualScaling)
      scaling->computeScaling(Problem);
    
    scaling->scaleLinearSystem(Problem);
  }


  if (precAlgorithm == AztecOO_) {
    
    // RPP: Aztec internal preconditioners can't be created until a
    // call to SetUserMatrix or SetUserOperator is made.  So we set
    // the operator here even though it is reset later.

    if (precMatrixSource == UseJacobian) {
      // The Jacobian has already been evaluated at the current solution.
      // Just set and enforce explicit constuction.
      this->setAztecOOJacobian();
      aztecSolverPtr->SetPrecMatrix(dynamic_cast<Epetra_RowMatrix*>(jacPtr.get()));
      aztecSolverPtr->ConstructPreconditioner(conditionNumberEstimate);
    }
    else if (precMatrixSource == SeparateMatrix) {
      Epetra_RowMatrix& precMatrix = dynamic_cast<Epetra_RowMatrix&>(*precPtr);
      precInterfacePtr->computePreconditioner(x.getEpetraVector(), 
					      *precPtr, &p);
      this->setAztecOOJacobian();
      aztecSolverPtr->SetPrecMatrix(&precMatrix);
      aztecSolverPtr->ConstructPreconditioner(conditionNumberEstimate);
    }

  }
  else if (precAlgorithm == Ifpack_) {
    
    if (precMatrixSource == UseJacobian) {
      createIfpackPreconditioner(p);
      solvePrecOpPtr = ifpackPreconditionerPtr;
      
    }
    else if (precMatrixSource == SeparateMatrix) {
      
      precInterfacePtr->computePreconditioner(x.getEpetraVector(),
					      *precPtr, &p);
      createIfpackPreconditioner(p);
      solvePrecOpPtr = ifpackPreconditionerPtr;
    }

  }
  else if (precAlgorithm == NewIfpack_) {
    
    if (precMatrixSource == UseJacobian) {
      createNewIfpackPreconditioner(p);
      solvePrecOpPtr = newIfpackPreconditionerPtr;
    }
    else if (precMatrixSource == SeparateMatrix) {
      
      precInterfacePtr->computePreconditioner(x.getEpetraVector(),
					      *precPtr, &p);
      createNewIfpackPreconditioner(p);
      solvePrecOpPtr = newIfpackPreconditionerPtr;
    }

  }
#ifdef HAVE_NOX_ML_EPETRA
  else if (precAlgorithm == ML_) {
    
    if (precMatrixSource == UseJacobian) {
      createMLPreconditioner(p);
      solvePrecOpPtr = MLPreconditionerPtr;
    }
    else if (precMatrixSource == SeparateMatrix) {
      
      precInterfacePtr->computePreconditioner(x.getEpetraVector(),
					      *precPtr, &p);
      createMLPreconditioner(p);
      solvePrecOpPtr = MLPreconditionerPtr;
    }

  }
#endif
  else if (precAlgorithm == UserDefined_) {

    precInterfacePtr->computePreconditioner(x.getEpetraVector(),
					    *precPtr, &p);
    solvePrecOpPtr = precPtr;

  }

  isPrecConstructed = true; 

  // Unscale the linear system
  if ( !Teuchos::is_null(scaling) )
    scaling->unscaleLinearSystem(Problem);

  double endTime = timer.WallTime();
  timeCreatePreconditioner += (endTime - startTime);

  if (utils.isPrintType(Utils::LinearSolverDetails))
    utils.out() << "\n       Time required to create preconditioner : " 
         << (endTime - startTime) << " (sec.)" << std::endl;;

  return true;
}

//***********************************************************************
bool NOX::Epetra::LinearSystemAztecOO::
createIfpackPreconditioner(Teuchos::ParameterList& p) const
{
  //for ifpack we need a VBR or CRS matrix to get the correct graph

  if ( !Teuchos::is_null(ifpackGraphPtr) ) 
    throwError("createIfpackPreconditioner", "Ifpack Graph NOT NULL");
  if ( !Teuchos::is_null(ifpackPreconditionerPtr)) 
    throwError("createIfpackPreconditioner", "Ifpack Prec NOT NULL");

  if (utils.isPrintType(Utils::Debug))
    utils.out() << "NOX::Epetra::LinearSolverAztecOO : createIfpackPrecon - \n"
         << "  using Fill Factor --> " << p.get("Fill Factor", 1)
         << std::endl;

  //check to see if it is a VBR matrix
  if (precType == EpetraVbrMatrix) {

    Epetra_VbrMatrix* vbr = 0;

    if (precMatrixSource == UseJacobian)
      vbr = dynamic_cast<Epetra_VbrMatrix*>(jacPtr.get());
    else if (precMatrixSource == SeparateMatrix)
      vbr = dynamic_cast<Epetra_VbrMatrix*>(precPtr.get());

    if (vbr == 0)
      throwError("createIfpackPreconditioner", 
		 "Dynamic cast to VBR Matrix failed!");

    ifpackGraphPtr = 
      Teuchos::rcp(new Ifpack_IlukGraph(vbr->Graph(),
					p.get("Fill Factor", 1),
					p.get("Overlap", 0)));
    ifpackGraphPtr->ConstructFilledGraph();
    ifpackPreconditionerPtr = 
      Teuchos::rcp(new Ifpack_CrsRiluk(*ifpackGraphPtr));
    ifpackPreconditionerPtr->SetAbsoluteThreshold(p.get("Absolute Threshold", 0.0));
    ifpackPreconditionerPtr->SetRelativeThreshold(p.get("Relative Threshold", 1.0));
    int err = ifpackPreconditionerPtr->InitValues(*vbr);
    if (err != 0)
      precError(err, "createIfpackPreconditioner()", "Ifpack", "InitValues");

    err = ifpackPreconditionerPtr->Factor();
    if (err != 0)
      precError(err, "createIfpackPreconditioner()", "Ifpack", "Factor");

  }

  // check to see if it is a Crs matrix
  else if (precType == EpetraCrsMatrix) {

    Epetra_CrsMatrix* crs = 0;

    if (precMatrixSource == UseJacobian)
      crs = dynamic_cast<Epetra_CrsMatrix*>(jacPtr.get());
    else if (precMatrixSource == SeparateMatrix)
      crs = dynamic_cast<Epetra_CrsMatrix*>(precPtr.get());

    if (crs == 0)
      throwError("createIfpackPreconditioner", 
		 "Dynamic cast to CRS Matrix failed!");

    ifpackGraphPtr = 
      Teuchos::rcp(new Ifpack_IlukGraph(crs->Graph(),
					p.get("Fill Factor", 1),
					p.get("Overlap", 0)));
    ifpackGraphPtr->ConstructFilledGraph();
    ifpackPreconditionerPtr = 
      Teuchos::rcp(new Ifpack_CrsRiluk(*ifpackGraphPtr));
    ifpackPreconditionerPtr->SetAbsoluteThreshold(p.get("Absolute Threshold", 0.0));
    ifpackPreconditionerPtr->SetRelativeThreshold(p.get("Relative Threshold", 1.0));
    int err = ifpackPreconditionerPtr->InitValues(*crs);
    if (err != 0)
      precError(err, "createIfpackPreconditioner()", "Ifpack", "InitValues");

    err = ifpackPreconditionerPtr->Factor();
    if (err != 0)
      precError(err, "createIfpackPreconditioner()", "Ifpack", "Factor");

  }
  
  // check to see if it is an operator that contains a Crs matrix
  else if (precType == EpetraRowMatrix) {

    // The following checks only for FiniteDiffernce based operators and
    // further that the underlying matrix is in Crs format......
    Epetra_CrsMatrix* crs = 0;
    NOX::Epetra::FiniteDifference* FDoperator = 0;

    if (precMatrixSource == UseJacobian)
      FDoperator = dynamic_cast<NOX::Epetra::FiniteDifference*>(jacPtr.get());
    else if (precMatrixSource == SeparateMatrix)
      FDoperator = dynamic_cast<NOX::Epetra::FiniteDifference*>(precPtr.get());

    if (FDoperator != 0)
      crs = &(FDoperator->getUnderlyingMatrix());

    if (crs == 0)
      throwError("createIfpackPreconditioner", 
		 "FiniteDifference: Underlying matrix NOT CRS Matrix!");

    ifpackGraphPtr = 
      Teuchos::rcp(new Ifpack_IlukGraph(crs->Graph(),
					p.get("Fill Factor", 1),
					p.get("Overlap", 0)));
    ifpackGraphPtr->ConstructFilledGraph();
    ifpackPreconditionerPtr = 
      Teuchos::rcp(new Ifpack_CrsRiluk(*ifpackGraphPtr));
    ifpackPreconditionerPtr->SetAbsoluteThreshold(p.get("Absolute Threshold", 0.0));
    ifpackPreconditionerPtr->SetRelativeThreshold(p.get("Relative Threshold", 1.0));
    int err = ifpackPreconditionerPtr->InitValues(*crs);
    if (err != 0)
      precError(err, "createIfpackPreconditioner()", "Ifpack", "InitValues");

    err = ifpackPreconditionerPtr->Factor();
    if (err != 0)
      precError(err, "createIfpackPreconditioner()", "Ifpack", "InitValues");

  }

  else {
  
    // If we made it this far, this routine should not have been called
    // in the first place.  An incorrect prec matrix object was supplied and
    // should have been caught in the checkOperatorValidity() method.
    throwError("createIfpackPreconditioner",
	       "Ifpack requires a CRS or VBR matrix!");
    return false;
  }

  // Print condition number estimate
  if (utils.isPrintType(Utils::Details)) {
    double kappa;
    ifpackPreconditionerPtr->Condest(false, kappa);
    utils.out() << 
      "\n       Condition number estimate of preconditioner is " << 
      utils.sciformat(kappa) << std::endl;
  }

  return true;
}

//***********************************************************************
bool NOX::Epetra::LinearSystemAztecOO::
createNewIfpackPreconditioner(Teuchos::ParameterList& p) const
{
  //for ifpack we need a VBR or CRS matrix to get the correct graph

  if ( !Teuchos::is_null(newIfpackPreconditionerPtr) ) 
    throwError("createNewIfpackPreconditioner", "Ifpack Prec NOT NULL");

  // Ensure we have a valid Teuchos parameter list to pass to Ifpack
  Teuchos::ParameterList& teuchosParams = 
    p.sublist("Ifpack");

  if (utils.isPrintType(Utils::Debug)) {
    utils.out() << "NOX::Epetra::LinearSolverAztecOO : createNewIfpackPrecon - \n"
         << "  using Teuchos parameter : " << std::endl;
    teuchosParams.print(utils.out());
  } 

  Ifpack Factory;

  //check to see if it is a row matrix
  if (precType == EpetraVbrMatrix || precType == EpetraCrsMatrix ||
      precType == EpetraRowMatrix) {

    Epetra_RowMatrix* row = 0;

    if (precMatrixSource == UseJacobian)
      row = dynamic_cast<Epetra_RowMatrix*>(jacPtr.get());
    else if (precMatrixSource == SeparateMatrix)
      row = dynamic_cast<Epetra_RowMatrix*>(precPtr.get());

    if (row == 0)
      throwError("createIfpackPreconditioner", 
		 "Dynamic cast to Row Matrix failed!");

    newIfpackPreconditionerPtr = Teuchos::rcp(Factory.Create(
      p.get("Ifpack Preconditioner", "ILU"), 
      row, 
      p.get("Overlap", 0),
      p.get("Override Serial Default", false) ));
    newIfpackPreconditionerPtr->SetParameters(teuchosParams);
    int err = newIfpackPreconditionerPtr->Initialize();
    if (err != 0)
      precError(err, "createNewIfpackPreconditioner()", "Ifpack", "Initialize");
    

    err = newIfpackPreconditionerPtr->Compute();
    if (err != 0)
      precError(err, "createNewIfpackPreconditioner()", "Ifpack", "Compute");

    return true;
  }
  
  // If we made it this far, this routine should not have been called
  // in the first place.  An incorrect prec matrix object was supplied and
  // should have been caught in the checkOperatorValidity() method.
  throwError("createNewIfpackPreconditioner",
	     "Ifpack requires a CRS or VBR matrix!");

  return false;
}

#ifdef HAVE_NOX_ML_EPETRA
//***********************************************************************
bool NOX::Epetra::LinearSystemAztecOO::
createMLPreconditioner(Teuchos::ParameterList& p) const
{
  //for ML we need a Epetra_RowMatrix

  if( !Teuchos::is_null(MLPreconditionerPtr) ) 
    throwError("createMLPreconditioner", "ML Prec NOT NULL");

  // Ensure we have a valid Teuchos parameter list to pass to ML
  Teuchos::ParameterList& teuchosParams = p.sublist("ML");

  if (utils.isPrintType(Utils::Debug)) 
  {
    utils.out() << "NOX::Epetra::LinearSolverAztecOO : createMLPreconditioner - \n"
         << "  using Teuchos parameter : " << std::endl;
    teuchosParams.print(utils.out());
  } 

  //check to see if it is a row matrix
  if (precType == EpetraVbrMatrix || precType == EpetraCrsMatrix ||
      precType == EpetraRowMatrix) {

    Epetra_RowMatrix* rowMatrix = 0;

    if (precMatrixSource == UseJacobian)
      rowMatrix = dynamic_cast<Epetra_RowMatrix*>(jacPtr.get());
    else if (precMatrixSource == SeparateMatrix)
      rowMatrix = dynamic_cast<Epetra_RowMatrix*>(precPtr.get());

    if (rowMatrix == 0)
      throwError("createMLPreconditioner", 
		 "Dynamic cast to Epetra_RowMatrix failed!");

    MLPreconditionerPtr = Teuchos::rcp( new ML_Epetra::MultiLevelPreconditioner(
                                  *rowMatrix, teuchosParams, true) );
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
bool NOX::Epetra::LinearSystemAztecOO::
recomputePreconditioner(const NOX::Epetra::Vector& x, 
			Teuchos::ParameterList& linearSolverParams) const
{  
  if (utils.isPrintType(Utils::LinearSolverDetails)) {
    utils.out() << "\n       Recomputing preconditioner" << std::endl;
  }

  if (precAlgorithm == None_) {
    return true;
  }
  else if (precAlgorithm == AztecOO_) {

    if (precMatrixSource == SeparateMatrix)
      precInterfacePtr->computePreconditioner(x.getEpetraVector(),
					      *precPtr, &linearSolverParams);

  }
  else if (precAlgorithm == Ifpack_) {
    // Original ifpack objects don't have a recompute option
    throwError("recomputePreconditioner",
	       "The preconditioner can only be \"Recomputed\" when using the \"NewIfpack\" or \"ML\" Preconditioner types!");
  }
  else if (precAlgorithm == NewIfpack_) {

    if (precMatrixSource == SeparateMatrix)
      precInterfacePtr->computePreconditioner(x.getEpetraVector(),
					      *precPtr, &linearSolverParams);

    int err = newIfpackPreconditionerPtr->Compute();
    if (err != 0)
      precError(err, "recomputePreconditioner", "Ifpack", "Compute");

  }
#ifdef HAVE_NOX_ML_EPETRA
  else if (precAlgorithm == ML_) {

    if (precMatrixSource == SeparateMatrix)
      precInterfacePtr->computePreconditioner(x.getEpetraVector(),
					      *precPtr, &linearSolverParams);
    
    int err = MLPreconditionerPtr->ReComputePreconditioner();

    if (err != 0)
      precError(err, "recomputePreconditioner", "ML", 
		"ReComputePreconditioner");
    

  }
#endif
  else if (precAlgorithm == UserDefined_) {
    
    if (precMatrixSource == SeparateMatrix)
      precInterfacePtr->computePreconditioner(x.getEpetraVector(),
					      *precPtr, &linearSolverParams);

  }

  return true;
}

//***********************************************************************
bool NOX::Epetra::LinearSystemAztecOO::destroyPreconditioner() const
{ 
  if (isPrecConstructed) {
    // For Aztec preconditioners. 
    if (precAlgorithm == AztecOO_) {
      aztecSolverPtr->DestroyPreconditioner();
    }
    // For IFPACK Preconditioners
    else if (precAlgorithm == Ifpack_) {
      ifpackPreconditionerPtr = Teuchos::null;
      ifpackGraphPtr = Teuchos::null;
    }
    else if (precAlgorithm == NewIfpack_) {
      newIfpackPreconditionerPtr = Teuchos::null;
    }
#ifdef HAVE_NOX_ML_EPETRA
    else if (precAlgorithm == ML_)
      MLPreconditionerPtr = Teuchos::null;
#endif
    if (utils.isPrintType(Utils::LinearSolverDetails)) {
      utils.out() << "\n       Destroying preconditioner" << std::endl;
    }
  }
  isPrecConstructed = false;
  solvePrecOpPtr = Teuchos::null;
  return true;
}

//***********************************************************************
NOX::Epetra::LinearSystemAztecOO::OperatorType 
NOX::Epetra::LinearSystemAztecOO::getOperatorType(const Epetra_Operator& Op)
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
bool NOX::Epetra::LinearSystemAztecOO::
computeJacobian(const NOX::Epetra::Vector& x)
{
  bool success = jacInterfacePtr->computeJacobian(x.getEpetraVector(), 
						  *jacPtr);
  return success;
}

//***********************************************************************
Teuchos::RCP<NOX::Epetra::Scaling>
NOX::Epetra::LinearSystemAztecOO::getScaling()
{
  return scaling;
}

//***********************************************************************
void NOX::Epetra::LinearSystemAztecOO::
resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& scalingObject)
{
  scaling = scalingObject;
  return;
}

//***********************************************************************
void NOX::Epetra::LinearSystemAztecOO::
throwError(const std::string& functionName, const std::string& errorMsg) const
{
  if (utils.isPrintType(Utils::Error)) {
    utils.out() << "NOX::Epetra::LinearSystemAztecOO::" << functionName 
	 << " - " << errorMsg << std::endl;
  }
  throw "NOX Error";
}

//***********************************************************************
Teuchos::RCP<const NOX::Epetra::Interface::Jacobian> 
NOX::Epetra::LinearSystemAztecOO::getJacobianInterface() const
{
  return jacInterfacePtr;
}

//***********************************************************************
Teuchos::RCP<const NOX::Epetra::Interface::Preconditioner> 
NOX::Epetra::LinearSystemAztecOO::getPrecInterface() const
{
  return precInterfacePtr;
}

//***********************************************************************
bool
NOX::Epetra::LinearSystemAztecOO::hasPreconditioner() const
{
  return (precAlgorithm != None_);
}

//***********************************************************************
bool
NOX::Epetra::LinearSystemAztecOO::isPreconditionerConstructed() const
{
  return isPrecConstructed;
}

//***********************************************************************
Teuchos::RCP<const Epetra_Operator>
NOX::Epetra::LinearSystemAztecOO::getJacobianOperator() const
{
  return jacPtr;
}

//***********************************************************************
Teuchos::RCP<Epetra_Operator>
NOX::Epetra::LinearSystemAztecOO::getJacobianOperator()
{
  return jacPtr;
}

//***********************************************************************
Teuchos::RCP<const Epetra_Operator>
NOX::Epetra::LinearSystemAztecOO::getPrecOperator() const
{
  return precPtr;
}

//***********************************************************************
Teuchos::RCP<const Epetra_Operator> 
NOX::Epetra::LinearSystemAztecOO::getGeneratedPrecOperator() const
{
  return solvePrecOpPtr;
}

//***********************************************************************
Teuchos::RCP<Epetra_Operator>
NOX::Epetra::LinearSystemAztecOO::getGeneratedPrecOperator()
{
  return solvePrecOpPtr;
}

//***********************************************************************
NOX::Epetra::LinearSystem::PreconditionerReusePolicyType 
NOX::Epetra::LinearSystemAztecOO::
getPreconditionerPolicy(bool advanceReuseCounter)
{

  if (precReusePolicy == PRPT_REBUILD) 
    return PRPT_REBUILD;

  if (precReusePolicy == PRPT_RECOMPUTE) {
    if (isPrecConstructed) 
      return PRPT_RECOMPUTE;
    else
      return PRPT_REBUILD;
  }

  // Below is for Russell's reuse of preconditioner - this toggles
  // between rebuild and reuse depending on how many times this
  // function has been called.

  if (precReusePolicy == PRPT_REUSE) {
    
    // If preconditioner is not built at all, you must build it
    if (!isPrecConstructed) {
      if (advanceReuseCounter)
	precQueryCounter++;
      return PRPT_REBUILD;
    }
    
    if (utils.isPrintType(Utils::Details)) 
      if (advanceReuseCounter)
	utils.out() << "\n       Preconditioner Reuse: Age of Prec --> " 
		    << precQueryCounter << " / " << maxAgeOfPrec << std::endl;
    
    // This allows reuse for the entire nonlinear solve
    if( maxAgeOfPrec == -2 ) {
      if (advanceReuseCounter)
	precQueryCounter++;
      return PRPT_REUSE;
    }
    
    // This allows one recompute of the preconditioner followed by reuse 
    // for the remainder of the nonlinear solve
    else if( maxAgeOfPrec == -1 ) {
      if (advanceReuseCounter)
	precQueryCounter++;
      maxAgeOfPrec = -2;
      return PRPT_REBUILD;
    }
    
    // This is the typical use 
    else {
      if( precQueryCounter == 0 || precQueryCounter >= maxAgeOfPrec ) {
        if (advanceReuseCounter)
          precQueryCounter = 1;
	return PRPT_REBUILD;
      }
      else {
	if (advanceReuseCounter)
	  precQueryCounter++;
	return PRPT_REUSE;
      }
    }
  } // if (precReusePolicy == PRPT_REUSE)
  
  return PRPT_REBUILD;
}

//***********************************************************************
double 
NOX::Epetra::LinearSystemAztecOO::getTimeCreatePreconditioner() const
{
  return timeCreatePreconditioner;
}

//***********************************************************************
double 
NOX::Epetra::LinearSystemAztecOO::getTimeApplyJacobianInverse() const
{
  return timeApplyJacbianInverse;
}

//***********************************************************************
void
NOX::Epetra::LinearSystemAztecOO::setJacobianOperatorForSolve(
	       const Teuchos::RCP<const Epetra_Operator>& solveJacOp)
{
  jacPtr = Teuchos::rcp_const_cast<Epetra_Operator>(solveJacOp);
  jacType = getOperatorType(*solveJacOp);
  this->setAztecOOJacobian();
}

//***********************************************************************
void
NOX::Epetra::LinearSystemAztecOO::setPrecOperatorForSolve(
	       const Teuchos::RCP<const Epetra_Operator>& solvePrecOp)
{
  solvePrecOpPtr = Teuchos::rcp_const_cast<Epetra_Operator>(solvePrecOp);
  this->setAztecOOPreconditioner();
}

//***********************************************************************
void
NOX::Epetra::LinearSystemAztecOO::setAztecOOJacobian() const
{
  if ((jacType == EpetraRowMatrix) ||
      (jacType == EpetraVbrMatrix) ||
      (jacType == EpetraCrsMatrix)) {
    aztecSolverPtr->SetUserMatrix(dynamic_cast<Epetra_RowMatrix*>(jacPtr.get()), false);
  }
  else
    aztecSolverPtr->SetUserOperator(jacPtr.get());
}

//***********************************************************************
void
NOX::Epetra::LinearSystemAztecOO::setAztecOOPreconditioner() const
{
  if ( !Teuchos::is_null(solvePrecOpPtr) && precAlgorithm != AztecOO_)
    aztecSolverPtr->SetPrecOperator(solvePrecOpPtr.get());
}

//***********************************************************************
void NOX::Epetra::LinearSystemAztecOO::
precError(int error_code, 
	  const std::string& nox_function,
	  const std::string& prec_type,
	  const std::string& prec_function) const
{
  if (error_code != 0) {
    
    std::ostringstream msg;

    if (throwErrorOnPrecFailure) 
      msg << "Error - ";
    else 
      msg << "Warning - ";

    msg << "NOX::Epetra::LinearSystemAztecOO::" << nox_function << " - The " << prec_type << " preconditioner has returned a nonzero error code of " << error_code << " for the function \"" << prec_function << "\".  Please consult the preconditioner documentation for this error type.";
    
    if (throwErrorOnPrecFailure) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
    }
    else
      utils.out() << msg.str() << std::endl; 
  }
}

//***********************************************************************

