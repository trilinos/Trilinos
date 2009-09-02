// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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

// NOX includes
#include "NOX_Epetra_Interface_Required.H"
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
#include "Epetra_RowMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"

// External include files for Stratimikos
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Teuchos_VerboseObject.hpp"

// EpetraExt includes for dumping a matrix
#ifdef HAVE_NOX_DEBUG
#ifdef HAVE_NOX_EPETRAEXT
#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_RowMatrixOut.h"
#endif
#endif

#ifdef HAVE_NOX_ML_EPETRA
#include "Teuchos_ParameterList.hpp"
#endif

#include <typeinfo>

//***********************************************************************
NOX::Epetra::LinearSystemStratimikos::
LinearSystemStratimikos(
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
  initializeStratimikos(linearSolverParams);
  tmpVectorPtr = Teuchos::rcp(new NOX::Epetra::Vector(cloneVector));

  // Create Jac Operator internally, if requested one of 2 ways
  if (jacPtr == Teuchos::null) 
    createJacobianOperator(printParams, linearSolverParams, iReq, cloneVector);
  else if (linearSolverParams.get("Jacobian Operator", "Have Jacobian")
           != "Have Jacobian")
    createJacobianOperator(printParams, linearSolverParams, iReq, cloneVector);
  
  jacType = getOperatorType(*jacPtr);
  precType = jacType;

  reset(linearSolverParams);
}

//***********************************************************************
NOX::Epetra::LinearSystemStratimikos::
LinearSystemStratimikos(
 Teuchos::ParameterList& printParams, 
 Teuchos::ParameterList& linearSolverParams,
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
  precPtr(preconditioner),
  scaling(s),
  conditionNumberEstimate(0.0),
  isPrecConstructed(false),
  precQueryCounter(0),
  maxAgeOfPrec(1),
  timer(cloneVector.getEpetraVector().Comm()),
  timeCreatePreconditioner(0.0),
  timeApplyJacbianInverse(0.0)
{
  // Interface for user-define preconditioning.
  if (precIsAlreadyInverted) precMatrixSource = UserDefined_;
  else                       precMatrixSource = SeparateMatrix;

  initializeStratimikos(linearSolverParams);
  tmpVectorPtr = Teuchos::rcp(new NOX::Epetra::Vector(cloneVector));

  // Both operators are supplied
  jacType = getOperatorType(*jacPtr);
  precType = getOperatorType(*precPtr);

  reset(linearSolverParams);
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
  Teuchos::RCP<Teuchos::FancyOStream> out
     = Teuchos::VerboseObjectBase::getDefaultOStream();

  Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

  linearSolverBuilder.setParameterList(Teuchos::rcp(&stratParams, false));

  // Create a linear solver factory given information read from the
  // parameter list.
  lowsFactory = linearSolverBuilder.createLinearSolveStrategy("");

  // Setup output stream and the verbosity level
  lowsFactory->setOStream(out);
  lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

}
//***********************************************************************
void NOX::Epetra::LinearSystemStratimikos::
reset(Teuchos::ParameterList& linearSolverParams)
{

  // First remove any preconditioner that may still be active
  destroyPreconditioner();
    
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
    string errorMessage = "Option for \"Preconditioner Reuse Policy\" is invalid! \nPossible options are \"Reuse\", \"Rebuild\", and \"Recompute\".";
    throwError("reset()", errorMessage);
  }
  maxAgeOfPrec = linearSolverParams.get("Max Age Of Prec", 1);
  precQueryCounter = 0;

#ifdef HAVE_NOX_DEBUG
#ifdef HAVE_NOX_EPETRAEXT
  linearSolveCount = 0;
#endif
#endif

}

//***********************************************************************
bool NOX::Epetra::LinearSystemStratimikos::createJacobianOperator(
       Teuchos::ParameterList& printParams,
       Teuchos::ParameterList& lsParams,
       const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
       const NOX::Epetra::Vector& cloneVector)
{
  string choice = lsParams.get("Jacobian Operator", "Matrix-Free");

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
  double startTime = timer.WallTime();

  // Need non-const version of the input vector
  // Epetra_LinearProblem requires non-const versions so we can perform
  // scaling of the linear problem.
  NOX::Epetra::Vector& nonConstInput = const_cast<NOX::Epetra::Vector&>(input);
  
  // Zero out the delta X of the linear problem if requested by user.
  if (zeroInitialGuess)
    result.init(0.0);

  // Wrap Thyra objects around Epetra and NOX objects
  Teuchos::RCP<const Thyra::LinearOpBase<double> > linearOp =
    Thyra::epetraLinearOp(jacPtr);

  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > lows;

  // Construct preconditioner if given AND if factory supports it
  if (precMatrixSource == SeparateMatrix &&
      lowsFactory->supportsPreconditionerInputType(
          Thyra::PRECONDITIONER_INPUT_TYPE_AS_MATRIX)) {

    Teuchos::RCP<const Thyra::LinearOpBase<double> > precOp =
       Thyra::epetraLinearOp(precPtr);

    lows = lowsFactory->createOp();
    Thyra::initializeApproxPreconditionedOp<double>(*lowsFactory,linearOp,precOp,&*lows);
   
  }
  else if (precMatrixSource == UserDefined_ &&
      lowsFactory->supportsPreconditionerInputType(
          Thyra::PRECONDITIONER_INPUT_TYPE_AS_OPERATOR)) {
    
    Teuchos::RCP<const Thyra::LinearOpBase<double> > precOp =
       Thyra::epetraLinearOp(precPtr);
    lows = lowsFactory->createOp();

    Thyra::initializePreconditionedOp<double>(
      *lowsFactory,linearOp,Thyra::rightPrec<double>(precOp) ,&*lows);
  }
  else {
    // Default case of using only the Jacobian and no 
    // separate preconditioner object
    lows = Thyra::linearOpWithSolve(*lowsFactory, linearOp);
  }
  

  Teuchos::RCP<Epetra_Vector> resultRCP =
    Teuchos::rcp(&result.getEpetraVector(), false);
  Teuchos::RCP<Epetra_Vector> inputRCP =
    Teuchos::rcp(&nonConstInput.getEpetraVector(), false);

  Teuchos::RCP<Thyra::VectorBase<double> >
    x = Thyra::create_Vector(resultRCP , linearOp->domain() );
  Teuchos::RCP<const Thyra::VectorBase<double> >
    b = Thyra::create_Vector(inputRCP, linearOp->range() );

  // Solve the linear system for x
  lows->solve(Thyra::NONCONJ_ELE, *b, &*x);

/*
  // Make sure preconditioner was constructed if requested
  if (!isPrecConstructed && (precAlgorithm != None_)) {
    throwError("applyJacobianInverse", 
       "Preconditioner is not constructed!  Call createPreconditioner() first.");
  }
*/

  // Set the output parameters in the "Output" sublist
/***
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
***/

  // Dump solution of linear system
#ifdef HAVE_NOX_DEBUG
#ifdef HAVE_NOX_EPETRAEXT
  ++linearSolveCount;
  std::ostringstream iterationNumber;
  iterationNumber << linearSolveCount;
    
  std::string prefixName = p.get("Write Linear System File Prefix", 
				 "NOX_LinSys");
  std::string postfixName = iterationNumber.str();
  postfixName += ".mm";
  if (p.get("Write Linear System", false)) {
    std::string lhsFileName = prefixName + "_LHS_" + postfixName;
    EpetraExt::MultiVectorToMatrixMarketFile(lhsFileName.c_str(), 
					   result.getEpetraVector());
  }
#endif
#endif

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
     << " NOT IMPLEMENTED " << endl;
return false;
/**
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
  else if (precAlgorithm == Stratimikos_) {

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
      TEST_FOR_EXCEPTION(true, std::logic_error, msg);
    }
    else {
      if (utils.isPrintType(NOX::Utils::Warning))
	utils.out() << msg << endl;
    }
    return false;
  }
  
  return true;
**/
}

//***********************************************************************
bool NOX::Epetra::LinearSystemStratimikos::
createPreconditioner(const NOX::Epetra::Vector& x, Teuchos::ParameterList& p, 
		     bool recomputeGraph) const
{
  double startTime = timer.WallTime();  

  if (utils.isPrintType(Utils::LinearSolverDetails))
    utils.out() << "\n       Creating a new preconditioner" << endl;;


  if (precMatrixSource == UseJacobian) {
    // Just set and enforce explicit constuction.
  }
  else if (precMatrixSource == SeparateMatrix) {
    //Epetra_RowMatrix& precMatrix = dynamic_cast<Epetra_RowMatrix&>(*precPtr);
    precInterfacePtr->computePreconditioner(x.getEpetraVector(), 
				      *precPtr, &p);
    //this->setAztecOOJacobian();
    //aztecSolverPtr->SetPrecMatrix(&precMatrix);
    //aztecSolverPtr->ConstructPreconditioner(conditionNumberEstimate);
  }
  else if (precMatrixSource == UserDefined_) {

    precInterfacePtr->computePreconditioner(x.getEpetraVector(),
					    *precPtr, &p);
    //solvePrecOpPtr = precPtr;

  }

  isPrecConstructed = true; 

  // Unscale the linear system
  //if ( !Teuchos::is_null(scaling) )
  //  scaling->unscaleLinearSystem(Problem);

  double endTime = timer.WallTime();
  timeCreatePreconditioner += (endTime - startTime);

  if (utils.isPrintType(Utils::LinearSolverDetails))
    utils.out() << "\n       Time required to create precondtioner : " 
         << (endTime - startTime) << " (sec.)" << endl;;

  return true;
}

//***********************************************************************
bool NOX::Epetra::LinearSystemStratimikos::
recomputePreconditioner(const NOX::Epetra::Vector& x, 
			Teuchos::ParameterList& linearSolverParams) const
{  
cout << " NOX::Epetra::LinearSystemStratimikos::applyRightPreconditioning\n"
     << " NOT IMPLEMENTED " << endl;
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
throwError(const string& functionName, const string& errorMsg) const
{
  if (utils.isPrintType(Utils::Error)) {
    utils.out() << "NOX::Epetra::LinearSystemStratimikos::" << functionName 
	 << " - " << errorMsg << endl;
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
  return true;
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
  return precPtr;
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
		    << precQueryCounter << " / " << maxAgeOfPrec << endl;
    
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
  /*
  solvePrecOpPtr = Teuchos::rcp_const_cast<Epetra_Operator>(solvePrecOp);
  this->setAztecOOPreconditioner();
  */
  TEST_FOR_EXCEPTION(true, std::logic_error,
     "NOX::Epetra::LinearSystemStratimikos::setPrecOperatorForSolve\n"
     << " NOT IMPLEMENTED ");
}

//***********************************************************************
void
NOX::Epetra::LinearSystemStratimikos::setStratimikosPreconditioner() const
{
  TEST_FOR_EXCEPTION(true, std::logic_error,
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
      TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
    }
    else
      utils.out() << msg.str() << endl; 
  }
}

//***********************************************************************

