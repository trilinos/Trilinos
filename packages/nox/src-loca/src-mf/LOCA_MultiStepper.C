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
#include "LOCA_MultiStepper.H"    // class definition
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_Combo.H"
#include "LOCA_StatusTest_Wrapper.H"

// LOCA Includes
#include "NOX_Utils.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_GlobalData.H"
#include "LOCA_Factory.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_MultiPredictor_AbstractStrategy.H"
#include "LOCA_MultiContinuation_AbstractStrategy.H"
#include "LOCA_MultiContinuation_AbstractGroup.H"  
#include "LOCA_MultiContinuation_ExtendedGroup.H"
#include "LOCA_MultiContinuation_ExtendedVector.H"
#include "LOCA_MultiContinuation_ConstrainedGroup.H"
#include "MFLOCA.H"

// Multifario Includes
extern "C" {
#include <MFAtlas.h>
}

LOCA::MultiStepper::MultiStepper(
            const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	   const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& initialGuess,
	   const Teuchos::RefCountPtr< NOX::StatusTest::Generic>& t,
	   const Teuchos::RefCountPtr<Teuchos::ParameterList>& p) :
  globalData(),
  parsedParams(),
  predictor(),
  curGroupPtr(),
  bifGroupPtr(),
  statusTestPtr(),
  paramListPtr(),
  stepperList(),
  solverPtr(),
  paramVec(),
  conParamIDVec(),
  conParamData()
{
  reset(global_data, initialGuess, t, p);
}

LOCA::MultiStepper::~MultiStepper() 
{ 
}

bool 
LOCA::MultiStepper::reset(
		  const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		  const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& initialGuess,
		  const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& t,
		  const Teuchos::RefCountPtr<Teuchos::ParameterList>& p) 
{
  globalData = global_data;
  paramListPtr = p;
  statusTestPtr = t;

  // Parse parameter list
  parsedParams = Teuchos::rcp(new LOCA::Parameter::SublistParser(globalData));
  parsedParams->parseSublists(paramListPtr);

  // Create predictor strategy
  Teuchos::RefCountPtr<Teuchos::ParameterList> predictorParams = 
    parsedParams->getSublist("Predictor");
  predictor = globalData->locaFactory->createPredictorStrategy(
							     parsedParams,
							     predictorParams);

  // Get stepper sublist
  stepperList = parsedParams->getSublist("Stepper");

  // Get continuation parameter vector
  paramVec = initialGuess->getParams();

  // Get continuation parameter info
  getConParamData();

  // Make a copy of the parameter list, change continuation method to
  // natural
  Teuchos::RefCountPtr<Teuchos::ParameterList> firstStepperParams = 
    Teuchos::rcp(new Teuchos::ParameterList(*stepperList));
  firstStepperParams->set("Continuation Method", "Natural");

  // Create constrained group
  Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> constraintsGrp
    = buildConstrainedGroup(initialGuess);

  // Create bifurcation group
  Teuchos::RefCountPtr<Teuchos::ParameterList> bifurcationParams = 
    parsedParams->getSublist("Bifurcation");
  bifGroupPtr = globalData->locaFactory->createBifurcationStrategy(
						       parsedParams,
						       bifurcationParams,
						       constraintsGrp);

  // Create continuation strategy
  curGroupPtr = globalData->locaFactory->createContinuationStrategy(
							parsedParams,
							firstStepperParams,
							bifGroupPtr, predictor,
							conParamIDVec);

  // Set step size			    
  for (unsigned int i=0; i<conParamIDVec.size(); i++)
    curGroupPtr->setStepSize(0.0, i);
  
  // Set previous solution vector in current solution group
  curGroupPtr->setPrevX(curGroupPtr->getX());

  // Create solver using initial conditions
  solverPtr = Teuchos::rcp(new NOX::Solver::Manager(
					   curGroupPtr, 
					   statusTestPtr,
					   parsedParams->getSublist("NOX")));

  printInitializationInfo();

  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters))
    paramListPtr->print(globalData->locaUtils->out());
  
  return true;
}

LOCA::Abstract::Iterator::IteratorStatus
LOCA::MultiStepper::run() {
  NOX::StatusTest::StatusType solverStatus;
  string callingFunction = "LOCA::MultiStepper::run()";

  // Perform solve of initial conditions
  solverStatus = solverPtr->solve();

  // Set up continuation groups
  const LOCA::MultiContinuation::ExtendedGroup& constSolnGrp =
    dynamic_cast<const LOCA::MultiContinuation::ExtendedGroup&>(
       solverPtr->getSolutionGroup());
  Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> underlyingGroup 
    = Teuchos::rcp_const_cast<LOCA::MultiContinuation::AbstractGroup>(constSolnGrp.getUnderlyingGroup());

  // Create continuation strategy
  curGroupPtr = globalData->locaFactory->createContinuationStrategy(
							parsedParams,
							stepperList,
							underlyingGroup, 
							predictor,
							conParamIDVec);

  // If nonlinear solve failed, return (this must be done after continuation 
  // groups are created so MultiStepper::getSolutionGroup() functions 
  // correctly.
  if (solverStatus != NOX::StatusTest::Converged)
    return LOCA::Abstract::Iterator::Failed;
  
  // Save initial solution
  curGroupPtr->printSolution();

  // Create new solver using new continuation groups and combo status test
  solverPtr = Teuchos::rcp(new NOX::Solver::Manager(
					   curGroupPtr, statusTestPtr,
					   parsedParams->getSublist("NOX")));

  MFImplicitMF M;
  MFNRegion Omega;
  MFAtlas A;
  MFNVector u0;
  MFContinuationMethod H;
  LOCAData* data = new LOCAData(globalData, parsedParams, solverPtr, 
				curGroupPtr, paramListPtr, 
				statusTestPtr, 
				Teuchos::rcp(&conParamData,false));
  M=MFIMFCreateLOCA(data);
  Omega=MFNRegionCreateLOCA(data);
  u0=MFCreateLOCANVectorWithData(Teuchos::rcp_dynamic_cast<LMCEV>(curGroupPtr->getX().clone()));

  H=MFCreateHendersonsMethod();

  /* Max distance from TS to M */
  MFHendersonsMethodSetEpsilon(H, stepperList->get("Epsilon", 0.1));

  /* -1 means infinite */
  MFHendersonsMethodSetMaxCharts(H,stepperList->get("Max Charts", -1));

  /* Write info to stdout */
  MFHendersonsMethodSetVerbose(H, stepperList->get("Verbosity", 1));

  /* Page out non-interior polyhedra */
  MFHendersonsMethodSetPage(H, stepperList->get("Page Charts", 1));

  /* Write polyhedra to a plotfile */
  MFHendersonsMethodSetDumpToPlotFile(H, stepperList->get("Dump Polyhedra", true));

  /* Write points to a file */
  MFHendersonsMethodSetDumpToCenterFile(H,stepperList->get("Dump Centers", false)); 

  /* File name to save data */
  const char *fname = 
    stepperList->get("Filename", "MFresults").c_str();
  MFHendersonsMethodSetFilename(H,const_cast<char*>(fname));

  A=MFComputeAtlas(H,M,Omega,u0);
  MFCloseAtlas(H,A);
  printf("\n\tDone computing Atlas\n");fflush(stdout);
  MFFreeAtlas(A);
  MFFreeImplicitMF(M);
  //MFFreeNRegion(Omega);
  MFFreeNVector(u0);
  MFFreeHendersonsMethod(H);

  return LOCA::Abstract::Iterator::Finished;
}

Teuchos::RefCountPtr<const LOCA::MultiContinuation::AbstractGroup>
LOCA::MultiStepper::getSolutionGroup()
{
  return curGroupPtr->getBaseLevelUnderlyingGroup();
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
LOCA::MultiStepper::getList() const
{
  return paramListPtr;
}

void 
LOCA::MultiStepper::printInitializationInfo()
{  
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
    globalData->locaUtils->out() 
      << std::endl 
      << globalData->locaUtils->fill(72, '~') 
      << std::endl;
    globalData->locaUtils->out() 
      << "Beginning Continuation Run \n" 
      << "Method: MultiParmeter Continuation \n"
      << std::endl;
    globalData->locaUtils->out() 
      << globalData->locaUtils->fill(72, '~') 
      << std::endl << std::endl;
  }
}

void 
LOCA::MultiStepper::getConParamData()
{
  string callingFunction = "LOCA::MultiStepper::getConParamInfo()";

  // Get number of continuation parameters
  int numParams = 
    stepperList->get("Number of Continuation Parameters",1);

  // Get data for each continuation parameter
  conParamData.clear();
  conParamIDVec.clear();
  for (int i=1; i<=numParams; i++) {

    // create sublist name
    stringstream sublistStream;
    sublistStream << "Continuation Parameter " << i;
    string sublistName = sublistStream.str();

    // Get sublist for continuation parameter
    if (!stepperList->isSublist(sublistName)) {
      stringstream errorStream;
      errorStream << "No sublist for continuation parameter " << i << "!";
      globalData->locaErrorCheck->throwError(callingFunction, 
					     errorStream.str());
    }
    Teuchos::ParameterList* conSublistPtr = &(stepperList->sublist(sublistName));

    // Create new struct for con param info
    ParamData d;

    // Get continuation parameter name
    if (!conSublistPtr->isParameter("Parameter Name")) {
      stringstream errorStream;
      errorStream << "\"Parameter Name\" for parameter sublist " 
                  << i << " is not set!";
      globalData->locaErrorCheck->throwError(callingFunction, 
					     errorStream.str());
    }
    d.name = conSublistPtr->get("Parameter Name", "None");
    d.ID = paramVec.getIndex(d.name);

    // Get initial value
    if (!conSublistPtr->isParameter("Initial Value")) {
      stringstream errorStream;
      errorStream << "\"Initial Value\" for parameter sublist " 
                  << i << " is not set!";
      globalData->locaErrorCheck->throwError(callingFunction, 
					     errorStream.str());
    }
    d.initialValue = conSublistPtr->get("Initial Value", 0.0);

    // Get max value
    if (!conSublistPtr->isParameter("Max Value")) {
      stringstream errorStream;
      errorStream << "\"Max Value\" for parameter sublist " 
                  << i << " is not set!";
      globalData->locaErrorCheck->throwError(callingFunction, 
					     errorStream.str());
    }
    d.maxValue = conSublistPtr->get("Max Value", 0.0);

    // Get min value
    if (!conSublistPtr->isParameter("Min Value")) {
      stringstream errorStream;
      errorStream << "\"Min Value\" for parameter sublist " 
                  << i << " is not set!";
      globalData->locaErrorCheck->throwError(callingFunction, 
					     errorStream.str());
    }
    d.minValue = conSublistPtr->get("Min Value", 0.0);

    // Get initial step size
    if (!conSublistPtr->isParameter("Initial Step Size")) {
      stringstream errorStream;
      errorStream << "\"Initial Step Size\" for parameter sublist " 
                  << i << " is not set!";
      globalData->locaErrorCheck->throwError(callingFunction, 
					     errorStream.str());
    }
    d.initialStepSize = conSublistPtr->get("Initial Step Size", 0.1);

    // Get initial max size
    if (!conSublistPtr->isParameter("Max Step Size")) {
      stringstream errorStream;
      errorStream << "\"Max Step Size\" for parameter sublist " 
                  << i << " is not set!";
      globalData->locaErrorCheck->throwError(callingFunction, 
					     errorStream.str());
    }
    d.maxStepSize = conSublistPtr->get("Max Step Size", 1.0);

    // Get initial min size
    if (!conSublistPtr->isParameter("Min Step Size")) {
      stringstream errorStream;
      errorStream << "\"Min Step Size\" for parameter sublist " 
                  << i << " is not set!";
      globalData->locaErrorCheck->throwError(callingFunction, 
					     errorStream.str());
    }
    d.minStepSize = conSublistPtr->get("Min Step Size", 1.0e-3);

    // add data struct to list
    conParamData.push_back(d);
    conParamIDVec.push_back(d.ID);
  }
}

Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>
LOCA::MultiStepper::buildConstrainedGroup(
      const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& grp)
{
  // Get constraints sublist
  Teuchos::RefCountPtr<Teuchos::ParameterList> constraintsList =
    parsedParams->getSublist("Constraints");

  // If we don't have a constraint object, return original group
  if (!constraintsList->isParameter("Constraint Object"))
    return grp;

  string methodName = "LOCA::MultiStepper::buildConstrainedGroup()";

  Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface> constraints;
  Teuchos::RefCountPtr< vector<string> > constraintParamNames;

  // Get constraint object
  if ((*constraintsList).INVALID_TEMPLATE_QUALIFIER
      isType< Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface> >("Constraint Object"))
    constraints = (*constraintsList).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface> >("Constraint Object");
  else
    globalData->locaErrorCheck->throwError(methodName,
	  "\"Constraint Object\" parameter is not of type Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface>!");

  // Get parameter names for constraints
  if ((*constraintsList).INVALID_TEMPLATE_QUALIFIER
      isType< Teuchos::RefCountPtr< vector<string> > >("Constraint Parameter Names"))
    constraintParamNames = (*constraintsList).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RefCountPtr< vector<string> > >("Constraint Parameter Names");
  else
    globalData->locaErrorCheck->throwError(methodName,
	  "\"Constraint Parameter Names\" parameter is not of type Teuchos::RefCountPtr< vector<string> >!");

  // Convert names to integer IDs
  vector<int> constraintParamIDs(constraintParamNames->size());
  const LOCA::ParameterVector& pvec = grp->getParams();
  for (unsigned int i=0; i<constraintParamIDs.size(); i++)
    constraintParamIDs[i] = pvec.getIndex((*constraintParamNames)[i]);

  // Create constrained group
  return 
    Teuchos::rcp(new LOCA::MultiContinuation::ConstrainedGroup(
							globalData,
							parsedParams,
							constraintsList,
							grp,
							constraints,
							constraintParamIDs));
}
