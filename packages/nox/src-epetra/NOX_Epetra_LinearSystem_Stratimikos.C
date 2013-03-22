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

#include "NOX_Epetra_LinearSystem_Stratimikos.H"	// class definition

#ifdef HAVE_NOX_STRATIMIKOS

// NOX includes
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "NOX_Epetra_MatrixFree.H"
#include "NOX_Epetra_FiniteDifference.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Epetra_Scaling.H"
#include "NOX_Utils.H"

// External include files for Epetra
#include "Epetra_Map.h"
#include "Epetra_Vector.h" 
#include "Epetra_Operator.h"
#include "Epetra_InvOperator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"

// External include files for Stratimikos
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_OperatorVectorTypes.hpp"

// EpetraExt includes for dumping a matrix
#ifdef HAVE_NOX_EPETRAEXT
#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_RowMatrixOut.h"
#endif

#ifdef HAVE_NOX_ML_EPETRA
#include "Teuchos_ParameterList.hpp"
#endif

#ifdef HAVE_NOX_TEKO
#include "Teko_StratimikosFactory.hpp"
#endif

#include "Teuchos_TimeMonitor.hpp"

#include <typeinfo>

//***********************************************************************
NOX::Epetra::LinearSystemStratimikos::
LinearSystemStratimikos(
 Teuchos::ParameterList& printParams, 
 Teuchos::ParameterList& stratSolverParams,  
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
  timeApplyJacbianInverse(0.0),
  getLinearSolveToleranceFromNox(false)
{
  jacType = getOperatorType(*jacPtr);
  precType = jacType;

  reset(stratSolverParams.sublist("NOX Stratimikos Options"));

  // Allocate solver
  initializeStratimikos(stratSolverParams.sublist("Stratimikos"));
  tmpVectorPtr = Teuchos::rcp(new NOX::Epetra::Vector(cloneVector));
}

//***********************************************************************
NOX::Epetra::LinearSystemStratimikos::
LinearSystemStratimikos(
 Teuchos::ParameterList& printParams, 
 Teuchos::ParameterList& stratSolverParams,
 const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac, 
 const Teuchos::RCP<Epetra_Operator>& jacobian,
 const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec, 
 const Teuchos::RCP<Epetra_Operator>& preconditioner,
 const NOX::Epetra::Vector& cloneVector,
 const bool& precIsAlreadyInverted,
 const Teuchos::RCP<NOX::Epetra::Scaling> s):
  utils(printParams),
  jacInterfacePtr(iJac),
  jacType(EpetraOperator),
  jacPtr(jacobian),
  precInterfacePtr(iPrec),
  precType(EpetraOperator),
  scaling(s),
  conditionNumberEstimate(0.0),
  isPrecConstructed(false),
  precQueryCounter(0),
  maxAgeOfPrec(1),
  timer(cloneVector.getEpetraVector().Comm()),
  timeCreatePreconditioner(0.0),
  timeApplyJacbianInverse(0.0),
  getLinearSolveToleranceFromNox(false)
{
  // Interface for user-defined preconditioning -- 
  // requires flipping of the apply and applyInverse methods
  precPtr = preconditioner;
  if (precIsAlreadyInverted) {
    precMatrixSource = UserDefined_;
  }
  else { // User supplies approximate matrix 
    precMatrixSource = SeparateMatrix;
  }

  // Both operators are supplied
  jacType = getOperatorType(*jacPtr);
  precType = getOperatorType(*precPtr);

  reset(stratSolverParams.sublist("NOX Stratimikos Options"));

  initializeStratimikos(stratSolverParams.sublist("Stratimikos"));
  tmpVectorPtr = Teuchos::rcp(new NOX::Epetra::Vector(cloneVector));
}

//***********************************************************************
NOX::Epetra::LinearSystemStratimikos::~LinearSystemStratimikos() 
{
  destroyPreconditioner();
}

//***********************************************************************
void NOX::Epetra::LinearSystemStratimikos::
initializeStratimikos(Teuchos::ParameterList& stratParams)
{
  Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

#ifdef HAVE_NOX_TEKO
  Teko::addTekoToStratimikosBuilder(linearSolverBuilder);
#endif

  linearSolverBuilder.setParameterList(Teuchos::rcp(&stratParams, false));

  // Create a linear solver factory given information read from the
  // parameter list.
  lowsFactory = linearSolverBuilder.createLinearSolveStrategy("");

  // Setup output stream and the verbosity level
  lowsFactory->setOStream(outputStream);
  lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

  // Initialize the LinearOpWithSolve
  lows = lowsFactory->createOp();

/*
  // Store sublist and std::string name where the Tolerance is set
  // so that inexact Newton can modify it in resetTolerance()
  if ("Belos" == stratParams.get("Linear Solver Type","Belos")) {
  }
*/
  

}
//***********************************************************************
void NOX::Epetra::LinearSystemStratimikos::
reset(Teuchos::ParameterList& noxStratParams)
{

  // First remove any preconditioner that may still be active
  destroyPreconditioner();
    
  zeroInitialGuess = 
    noxStratParams.get("Zero Initial Guess", false);

  manualScaling = 
    noxStratParams.get("Compute Scaling Manually", true);

  // Place linear solver details in the "Output" sublist of the
  // "Linear Solver" parameter list
  outputSolveDetails = 
    noxStratParams.get("Output Solver Details", true);

  throwErrorOnPrecFailure = 
    noxStratParams.get("Throw Error on Prec Failure", true);

  // Setup the preconditioner reuse policy
  std::string preReusePolicyName = 
    noxStratParams.get("Preconditioner Reuse Policy", "Rebuild");
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
  maxAgeOfPrec = noxStratParams.get("Max Age Of Prec", 1);
  precQueryCounter = 0;

  // This needs to be in sync with Inexact Newton option from
  // a different sublist
  getLinearSolveToleranceFromNox = 
    noxStratParams.get("Use Linear Solve Tolerance From NOX", false);

  linearSolveCount = 0;
  linearSolveIters_last = 0;
  linearSolveIters_total = 0;
  linearSolveAchievedTol = 0.0;

  if (noxStratParams.isParameter("Output Stream"))
    outputStream = 
      noxStratParams.get< Teuchos::RCP<Teuchos::FancyOStream> >("Output Stream");
  else
    outputStream = Teuchos::VerboseObjectBase::getDefaultOStream();

}

//***********************************************************************
bool NOX::Epetra::LinearSystemStratimikos::
applyJacobian(const NOX::Epetra::Vector& input, 
	      NOX::Epetra::Vector& result) const
{
  jacPtr->SetUseTranspose(false);
  int status = jacPtr->Apply(input.getEpetraVector(), 
  		  	     result.getEpetraVector());
  return (status == 0);
}

//***********************************************************************
bool NOX::Epetra::LinearSystemStratimikos::
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
bool NOX::Epetra::LinearSystemStratimikos::
applyJacobianInverse(Teuchos::ParameterList &p,
		     const NOX::Epetra::Vector& input, 
		     NOX::Epetra::Vector& result)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  NOX_FUNC_TIME_MONITOR("NOX: Total Linear Solve Time");

  double startTime = timer.WallTime();

  // Need non-const version of the input vector
  // Epetra_LinearProblem requires non-const versions so we can perform
  // scaling of the linear problem.
  NOX::Epetra::Vector& nonConstInput = const_cast<NOX::Epetra::Vector&>(input);
  
  // Zero out the delta X of the linear problem if requested by user.
  if (zeroInitialGuess) result.init(0.0);

  // Wrap Thyra objects around Epetra and NOX objects
  Teuchos::RCP<const Thyra::LinearOpBase<double> > linearOp =
    Thyra::epetraLinearOp(jacPtr);

  // Set the linear Op and  precomputed prec on this lows
  if (precObj == Teuchos::null) 
    Thyra::initializeOp(*lowsFactory, linearOp, lows.ptr());
  else 
    Thyra::initializePreconditionedOp<double>(
      *lowsFactory, linearOp, precObj, lows.ptr());

  Teuchos::RCP<Epetra_Vector> resultRCP =
    Teuchos::rcp(&result.getEpetraVector(), false);
  Teuchos::RCP<Epetra_Vector> inputRCP =
    Teuchos::rcp(&nonConstInput.getEpetraVector(), false);

  Teuchos::RCP<Thyra::VectorBase<double> >
    x = Thyra::create_Vector(resultRCP , linearOp->domain() );
  Teuchos::RCP<const Thyra::VectorBase<double> >
    b = Thyra::create_Vector(inputRCP, linearOp->range() );

  // Alter the convergence tolerance, if Inexact Newton
  Teuchos::RCP<Thyra::SolveCriteria<double> > solveCriteria;
  if (getLinearSolveToleranceFromNox) {
    Thyra::SolveMeasureType solveMeasure(
        Thyra::SOLVE_MEASURE_NORM_RESIDUAL,
        Thyra::SOLVE_MEASURE_NORM_INIT_RESIDUAL );
    solveCriteria = Teuchos::rcp(new Thyra::SolveCriteria<double>(
        solveMeasure, p.get<double>("Tolerance") ) );
  }

  // Solve the linear system for x
  Thyra::SolveStatus<double> status =
    lows->solve(Thyra::NOTRANS, *b, x.ptr(), solveCriteria.ptr());

  // MOVE TO FUNCTION: Update statistics: solves, iters, iters_total, achieved tol
  ++linearSolveCount;
  if (status.extraParameters != Teuchos::null) {
    if (status.extraParameters->isParameter("Belos/Iteration Count")) {
      linearSolveIters_last = status.extraParameters->get<int>("Belos/Iteration Count");
      linearSolveIters_total += linearSolveIters_last;
    }
    if (status.extraParameters->isParameter("Belos/Achieved Tolerance")) 
      linearSolveAchievedTol = status.extraParameters->get<double>("Belos/Achieved Tolerance");
    if (status.extraParameters->isParameter("AztecOO/Iteration Count")) {

      linearSolveIters_last = status.extraParameters->get<int>("AztecOO/Iteration Count");
      linearSolveIters_total += linearSolveIters_last;
    }
    if (status.extraParameters->isParameter("AztecOO/Achieved Tolerance")) 
      linearSolveAchievedTol = status.extraParameters->get<double>("AztecOO/Achieved Tolerance");
  }

  // Dump solution of linear system
#ifdef HAVE_NOX_EPETRAEXT
  if (p.get("Write Linear System", false)) {
    std::ostringstream iterationNumber;
    iterationNumber << linearSolveCount;
    
    std::string prefixName = p.get("Write Linear System File Prefix", 
				   "NOX_LinSys");
    std::string postfixName = iterationNumber.str();
    postfixName += ".mm";

    std::string lhsFileName = prefixName + "_LHS_" + postfixName;
    std::string rhsFileName = prefixName + "_RHS_" + postfixName;
    std::string jacFileName = prefixName + "_Jacobian_" + postfixName;
    EpetraExt::MultiVectorToMatrixMarketFile(lhsFileName.c_str(), 
					     result.getEpetraVector());
    EpetraExt::MultiVectorToMatrixMarketFile(rhsFileName.c_str(), 
					     input.getEpetraVector());

    Epetra_RowMatrix* printMatrix = NULL;
    printMatrix = dynamic_cast<Epetra_RowMatrix*>(jacPtr.get()); 
    if (printMatrix == NULL) {
      std::cout << "Error: NOX::Epetra::LinearSystemAztecOO::applyJacobianInverse() - "
	   << "Could not cast the Jacobian operator to an Epetra_RowMatrix!"
	   << "Please set the \"Write Linear System\" parameter to false."
	   << std::endl;
      throw "NOX Error";
    }
    EpetraExt::RowMatrixToMatrixMarketFile(jacFileName.c_str(), *printMatrix, 
					   "test matrix", "Jacobian XXX");
  }
#endif

  //Release RCPs
  x = Teuchos::null; b = Teuchos::null; 
  resultRCP = Teuchos::null; inputRCP = Teuchos::null; 


  double endTime = timer.WallTime();
  timeApplyJacbianInverse += (endTime - startTime);

  return true;
}

//***********************************************************************
bool NOX::Epetra::LinearSystemStratimikos::
applyRightPreconditioning(bool useTranspose, 
			  Teuchos::ParameterList& params,
			  const NOX::Epetra::Vector& input, 
			  NOX::Epetra::Vector& result) const
{

cout << " NOX::Epetra::LinearSystemStratimikos::applyRightPreconditioning\n"
     << " NOT IMPLEMENTED " << std::endl;
return false;
}

//***********************************************************************
bool NOX::Epetra::LinearSystemStratimikos::
createPreconditioner(const NOX::Epetra::Vector& x, Teuchos::ParameterList& p, 
		     bool recomputeGraph) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  NOX_FUNC_TIME_MONITOR("NOX: Total Preconditioner Generation Time");

  double startTime = timer.WallTime();  

  if (utils.isPrintType(Utils::LinearSolverDetails))
    utils.out() << "\n       Creating a new preconditioner" << std::endl;;


//  Teuchos::RCP<Thyra::PreconditionerBase<double> > precObj;

  if (precMatrixSource == UseJacobian) {
    RCP<Thyra::PreconditionerFactoryBase<double> > precFactory =
      lowsFactory->getPreconditionerFactory();

    if (precFactory != Teuchos::null) {
      precObj = precFactory->createPrec();

      // Wrap Thyra objects around Epetra Op
      RCP<const Thyra::LinearOpBase<double> > linearOp =
        Thyra::epetraLinearOp(jacPtr);
  
      RCP<const Thyra::LinearOpSourceBase<double> > losb =
           rcp(new Thyra::DefaultLinearOpSource<double>(linearOp));

      // Computation of prec (e.g. ilu) happens here:
      precFactory->initializePrec(losb, precObj.get());

      // Get underlying Epetra operator
      Teuchos::RCP<Thyra::LinearOpBase<double> > pop;
      pop = precObj->getNonconstRightPrecOp();
      if (pop == Teuchos::null)
      	pop = precObj->getNonconstUnspecifiedPrecOp();
      solvePrecOpPtr = 
      	Teuchos::rcp_dynamic_cast<Thyra::EpetraLinearOp>(pop,true)->epetra_op();
    }
    else // no preconditioner
      precObj = Teuchos::null;

  }
  else if (precMatrixSource == SeparateMatrix) {

    // Compute matrix for use in preconditioning
    precInterfacePtr->computePreconditioner(x.getEpetraVector(), 
				      *precPtr, &p);

    RCP<Thyra::PreconditionerFactoryBase<double> > precFactory =
      lowsFactory->getPreconditionerFactory();

    if (precFactory != Teuchos::null) {
      precObj = precFactory->createPrec();

      // Send Prec Matrix to 
      RCP<const Thyra::LinearOpBase<double> > precOp =
        Thyra::epetraLinearOp(precPtr);

      RCP<const Thyra::LinearOpSourceBase<double> > losb =
           rcp(new Thyra::DefaultLinearOpSource<double>(precOp));

      // Computation of prec (e.g. ilu) happens here:
      precFactory->initializePrec(losb, precObj.get());

      // Get underlying Epetra operator
      Teuchos::RCP<Thyra::LinearOpBase<double> > pop;
      pop = precObj->getNonconstRightPrecOp();
      if (pop == Teuchos::null)
      	pop = precObj->getNonconstUnspecifiedPrecOp();
      solvePrecOpPtr = 
      	Teuchos::rcp_dynamic_cast<Thyra::EpetraLinearOp>(pop,true)->epetra_op();
    }
    else // no preconditioner
      precObj = Teuchos::null;

  }
  else if (precMatrixSource == UserDefined_) {

    precInterfacePtr->computePreconditioner(x.getEpetraVector(),
					    *precPtr, &p);

    // Wrap the preconditioner so that apply() calls ApplyInverse()
    RCP<const Thyra::LinearOpBase<double> > precOp =
      Thyra::epetraLinearOp(precPtr, 
			    Thyra::NOTRANS, 
			    Thyra::EPETRA_OP_APPLY_APPLY_INVERSE);

    RCP<Thyra::DefaultPreconditioner<double> > precObjDef =
       rcp(new Thyra::DefaultPreconditioner<double>);
    precObjDef->initializeRight(precOp);
    precObj = precObjDef;
    solvePrecOpPtr = precPtr;
  }

  isPrecConstructed = true; 

  // Unscale the linear system
  //if ( !Teuchos::is_null(scaling) )
  //  scaling->unscaleLinearSystem(Problem);

  double endTime = timer.WallTime();
  timeCreatePreconditioner += (endTime - startTime);

  if (precObj != Teuchos::null)  {
    if (utils.isPrintType(Utils::LinearSolverDetails))
      utils.out() << "\n       Time required to create preconditioner : " 
           << (endTime - startTime) << " (sec.)" << std::endl;
  }
  else {
    if (utils.isPrintType(Utils::LinearSolverDetails))
      utils.out() << "\n       No preconditioner requested. " << std::endl;
  }

  return true;
}

//***********************************************************************
bool NOX::Epetra::LinearSystemStratimikos::
recomputePreconditioner(const NOX::Epetra::Vector& x, 
			Teuchos::ParameterList& linearSolverParams) const
{  
cout << " NOX::Epetra::LinearSystemStratimikos::recomputePreconditioner\n"
     << " NOT IMPLEMENTED " << std::endl;
return false;
}

//***********************************************************************
bool NOX::Epetra::LinearSystemStratimikos::destroyPreconditioner() const
{ 
  return true;
}

//***********************************************************************
NOX::Epetra::LinearSystemStratimikos::OperatorType 
NOX::Epetra::LinearSystemStratimikos::getOperatorType(const Epetra_Operator& Op)
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
bool NOX::Epetra::LinearSystemStratimikos::
computeJacobian(const NOX::Epetra::Vector& x)
{
  bool success = jacInterfacePtr->computeJacobian(x.getEpetraVector(), 
						  *jacPtr);
  return success;
}

//***********************************************************************
Teuchos::RCP<NOX::Epetra::Scaling>
NOX::Epetra::LinearSystemStratimikos::getScaling()
{
  return scaling;
}

//***********************************************************************
void NOX::Epetra::LinearSystemStratimikos::
resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& scalingObject)
{
  scaling = scalingObject;
  return;
}

//***********************************************************************
void NOX::Epetra::LinearSystemStratimikos::
throwError(const std::string& functionName, const std::string& errorMsg) const
{
  if (utils.isPrintType(Utils::Error)) {
    utils.out() << "NOX::Epetra::LinearSystemStratimikos::" << functionName 
	 << " - " << errorMsg << std::endl;
  }
  throw "NOX Error";
}

//***********************************************************************
Teuchos::RCP<const NOX::Epetra::Interface::Jacobian> 
NOX::Epetra::LinearSystemStratimikos::getJacobianInterface() const
{
  return jacInterfacePtr;
}

//***********************************************************************
Teuchos::RCP<const NOX::Epetra::Interface::Preconditioner> 
NOX::Epetra::LinearSystemStratimikos::getPrecInterface() const
{
  return precInterfacePtr;
}

//***********************************************************************
bool
NOX::Epetra::LinearSystemStratimikos::hasPreconditioner() const
{
  return precObj != Teuchos::null;
}

//***********************************************************************
bool
NOX::Epetra::LinearSystemStratimikos::isPreconditionerConstructed() const
{
  return isPrecConstructed;
}

//***********************************************************************
Teuchos::RCP<const Epetra_Operator>
NOX::Epetra::LinearSystemStratimikos::getJacobianOperator() const
{
  return jacPtr;
}

//***********************************************************************
Teuchos::RCP<Epetra_Operator>
NOX::Epetra::LinearSystemStratimikos::getJacobianOperator()
{
  return jacPtr;
}

//***********************************************************************
Teuchos::RCP<const Epetra_Operator>
NOX::Epetra::LinearSystemStratimikos::getPrecOperator() const
{
  return solvePrecOpPtr;
}

//***********************************************************************
Teuchos::RCP<const Epetra_Operator> 
NOX::Epetra::LinearSystemStratimikos::getGeneratedPrecOperator() const
{
  return solvePrecOpPtr;
}

//***********************************************************************
Teuchos::RCP<Epetra_Operator>
NOX::Epetra::LinearSystemStratimikos::getGeneratedPrecOperator()
{
  return solvePrecOpPtr;
}

//***********************************************************************
NOX::Epetra::LinearSystem::PreconditionerReusePolicyType 
NOX::Epetra::LinearSystemStratimikos::
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
NOX::Epetra::LinearSystemStratimikos::getTimeCreatePreconditioner() const
{
  return timeCreatePreconditioner;
}

//***********************************************************************
double 
NOX::Epetra::LinearSystemStratimikos::getTimeApplyJacobianInverse() const
{
  return timeApplyJacbianInverse;
}

//***********************************************************************
void
NOX::Epetra::LinearSystemStratimikos::setJacobianOperatorForSolve(
	       const Teuchos::RCP<const Epetra_Operator>& solveJacOp)
{
  jacPtr = Teuchos::rcp_const_cast<Epetra_Operator>(solveJacOp);
  jacType = getOperatorType(*solveJacOp);
}

//***********************************************************************
void
NOX::Epetra::LinearSystemStratimikos::setPrecOperatorForSolve(
	       const Teuchos::RCP<const Epetra_Operator>& solvePrecOp)
{
  solvePrecOpPtr = Teuchos::rcp_const_cast<Epetra_Operator>(solvePrecOp);

  // Wrap the preconditioner so that apply() calls ApplyInverse()
  Teuchos::RCP<const Thyra::LinearOpBase<double> > precOp =
    Thyra::epetraLinearOp(solvePrecOpPtr, 
			  Thyra::NOTRANS, 
			  Thyra::EPETRA_OP_APPLY_APPLY_INVERSE);

  Teuchos::RCP<Thyra::DefaultPreconditioner<double> > precObjDef =
    rcp(new Thyra::DefaultPreconditioner<double>);
  precObjDef->initializeRight(precOp);
  precObj = precObjDef;
  
}

//***********************************************************************
void
NOX::Epetra::LinearSystemStratimikos::setStratimikosPreconditioner() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
     "NOX::Epetra::LinearSystemStratimikos::setStratimikosPreconditioner()\n"
     << " NOT IMPLEMENTED ");
}

//***********************************************************************
void NOX::Epetra::LinearSystemStratimikos::
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

    msg << "NOX::Epetra::LinearSystemStratimikos::" << nox_function << " - The " << prec_type << " preconditioner has returned a nonzero error code of " << error_code << " for the function \"" << prec_function << "\".  Please consult the preconditioner documentation for this error type.";
    
    if (throwErrorOnPrecFailure) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
    }
    else
      utils.out() << msg.str() << std::endl; 
  }
}

//***********************************************************************

int NOX::Epetra::LinearSystemStratimikos::
getNumLinearSolves() {return linearSolveCount;}

int NOX::Epetra::LinearSystemStratimikos::
getLinearItersLastSolve() {return linearSolveIters_last;}
  
int NOX::Epetra::LinearSystemStratimikos::
getLinearItersTotal() {return linearSolveIters_total;}
  
double NOX::Epetra::LinearSystemStratimikos::
getAchievedTol() {return linearSolveAchievedTol;}
//***********************************************************************

#endif //HAVE_NOX_STRATIMIKOS
