// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
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
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
       const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& initialGuess,
       const Teuchos::RCP< NOX::StatusTest::Generic>& t,
       const Teuchos::RCP<Teuchos::ParameterList>& p) :
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
          const Teuchos::RCP<LOCA::GlobalData>& global_data,
          const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& initialGuess,
          const Teuchos::RCP<NOX::StatusTest::Generic>& t,
          const Teuchos::RCP<Teuchos::ParameterList>& p)
{
  globalData = global_data;
  paramListPtr = p;
  statusTestPtr = t;

  // Parse parameter list
  parsedParams = Teuchos::rcp(new LOCA::Parameter::SublistParser(globalData));
  parsedParams->parseSublists(paramListPtr);

  // Create predictor strategy
  Teuchos::RCP<Teuchos::ParameterList> predictorParams =
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
  Teuchos::RCP<Teuchos::ParameterList> firstStepperParams =
    Teuchos::rcp(new Teuchos::ParameterList(*stepperList));
  firstStepperParams->set("Continuation Method", "Natural");

  // Create constrained group
  Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> constraintsGrp
    = buildConstrainedGroup(initialGuess);

  // Create bifurcation group
  Teuchos::RCP<Teuchos::ParameterList> bifurcationParams =
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
  solverPtr = NOX::Solver::buildSolver(curGroupPtr,
                       statusTestPtr,
                       parsedParams->getSublist("NOX"));

  printInitializationInfo();

  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters))
    paramListPtr->print(globalData->locaUtils->out());

  return true;
}

LOCA::Abstract::Iterator::IteratorStatus
LOCA::MultiStepper::run() {
  NOX::StatusTest::StatusType solverStatus;
  std::string callingFunction = "LOCA::MultiStepper::run()";

  // Perform solve of initial conditions
  solverStatus = solverPtr->solve();

  // Set up continuation groups
  const LOCA::MultiContinuation::ExtendedGroup& constSolnGrp =
    dynamic_cast<const LOCA::MultiContinuation::ExtendedGroup&>(
       solverPtr->getSolutionGroup());
  Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> underlyingGroup
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
  solverPtr = NOX::Solver::buildSolver(curGroupPtr, statusTestPtr,
                       parsedParams->getSublist("NOX"));

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
  u0=MFCreateLOCANVectorWithData(Teuchos::rcp_dynamic_cast<LMCEV>(curGroupPtr->getX().clone()), data->mfErrorHandler);

  H=MFCreateMultifariosMethod(data->mfErrorHandler);

  /* Max distance from TS to M */
  MFMultifarioSetRealParameter(H, "epsilon", stepperList->get("Epsilon", 0.1), data->mfErrorHandler);

  /* -1 means infinite */
  MFMultifarioSetIntegerParameter(H, "maxCharts", stepperList->get("Max Charts", -1), data->mfErrorHandler);

  /* Write info to stdout */
  MFMultifarioSetIntegerParameter(H, "verbose", stepperList->get("Verbosity", 1), data->mfErrorHandler);

  /* Page out non-interior polyhedra */
  MFMultifarioSetIntegerParameter(H, "page", stepperList->get("Page Charts", 1), data->mfErrorHandler);

  /* Write polyhedra to a plotfile */
  MFMultifarioSetIntegerParameter(H, "dumpToPlotFile", stepperList->get("Dump Polyhedra", true), data->mfErrorHandler);

  /* Write points to a file */
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile", stepperList->get("Dump Centers", false), data->mfErrorHandler);

  /* File name to save data */
  const char *fname =
    stepperList->get("Filename", "MFresults").c_str();
  MFMultifarioSetFilename(H,const_cast<char*>(fname), data->mfErrorHandler);

  A=MFComputeAtlas(H,M,Omega,u0,data->mfErrorHandler);
  MFCloseAtlas(H,A,data->mfErrorHandler);
  printf("\n\tDone computing Atlas\n");fflush(stdout);
  MFFreeAtlas(A,data->mfErrorHandler);
  MFFreeImplicitMF(M,data->mfErrorHandler);
  //MFFreeNRegion(Omega);
  MFFreeNVector(u0,data->mfErrorHandler);
  MFFreeContinuationMethod(H,data->mfErrorHandler);

  return LOCA::Abstract::Iterator::Finished;
}

Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
LOCA::MultiStepper::getSolutionGroup()
{
  return curGroupPtr->getBaseLevelUnderlyingGroup();
}

Teuchos::RCP<const Teuchos::ParameterList>
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
  std::string callingFunction = "LOCA::MultiStepper::getConParamInfo()";

  // Get number of continuation parameters
  int numParams =
    stepperList->get("Number of Continuation Parameters",1);

  // Get data for each continuation parameter
  conParamData.clear();
  conParamIDVec.clear();
  for (int i=1; i<=numParams; i++) {

    // create sublist name
    std::stringstream sublistStream;
    sublistStream << "Continuation Parameter " << i;
    std::string sublistName = sublistStream.str();

    // Get sublist for continuation parameter
    if (!stepperList->isSublist(sublistName)) {
      std::stringstream errorStream;
      errorStream << "No sublist for continuation parameter " << i << "!";
      globalData->locaErrorCheck->throwError(callingFunction,
                         errorStream.str());
    }
    Teuchos::ParameterList* conSublistPtr = &(stepperList->sublist(sublistName));

    // Create new struct for con param info
    ParamData d;

    // Get continuation parameter name
    if (!conSublistPtr->isParameter("Parameter Name")) {
      std::stringstream errorStream;
      errorStream << "\"Parameter Name\" for parameter sublist "
                  << i << " is not set!";
      globalData->locaErrorCheck->throwError(callingFunction,
                         errorStream.str());
    }
    d.name = conSublistPtr->get("Parameter Name", "None");
    d.ID = paramVec.getIndex(d.name);

    // Get initial value
    if (!conSublistPtr->isParameter("Initial Value")) {
      std::stringstream errorStream;
      errorStream << "\"Initial Value\" for parameter sublist "
                  << i << " is not set!";
      globalData->locaErrorCheck->throwError(callingFunction,
                         errorStream.str());
    }
    d.initialValue = conSublistPtr->get("Initial Value", 0.0);

    // Get max value
    if (!conSublistPtr->isParameter("Max Value")) {
      std::stringstream errorStream;
      errorStream << "\"Max Value\" for parameter sublist "
                  << i << " is not set!";
      globalData->locaErrorCheck->throwError(callingFunction,
                         errorStream.str());
    }
    d.maxValue = conSublistPtr->get("Max Value", 0.0);

    // Get min value
    if (!conSublistPtr->isParameter("Min Value")) {
      std::stringstream errorStream;
      errorStream << "\"Min Value\" for parameter sublist "
                  << i << " is not set!";
      globalData->locaErrorCheck->throwError(callingFunction,
                         errorStream.str());
    }
    d.minValue = conSublistPtr->get("Min Value", 0.0);

    // Get initial step size
    if (!conSublistPtr->isParameter("Initial Step Size")) {
      std::stringstream errorStream;
      errorStream << "\"Initial Step Size\" for parameter sublist "
                  << i << " is not set!";
      globalData->locaErrorCheck->throwError(callingFunction,
                         errorStream.str());
    }
    d.initialStepSize = conSublistPtr->get("Initial Step Size", 0.1);

    // Get initial max size
    if (!conSublistPtr->isParameter("Max Step Size")) {
      std::stringstream errorStream;
      errorStream << "\"Max Step Size\" for parameter sublist "
                  << i << " is not set!";
      globalData->locaErrorCheck->throwError(callingFunction,
                         errorStream.str());
    }
    d.maxStepSize = conSublistPtr->get("Max Step Size", 1.0);

    // Get initial min size
    if (!conSublistPtr->isParameter("Min Step Size")) {
      std::stringstream errorStream;
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

Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
LOCA::MultiStepper::buildConstrainedGroup(
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp)
{
  // Get constraints sublist
  Teuchos::RCP<Teuchos::ParameterList> constraintsList =
    parsedParams->getSublist("Constraints");

  // If we don't have a constraint object, return original group
  if (!constraintsList->isParameter("Constraint Object"))
    return grp;

  std::string methodName = "LOCA::MultiStepper::buildConstrainedGroup()";

  Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface> constraints;
  Teuchos::RCP< std::vector<std::string> > constraintParamNames;

  // Get constraint object
  if ((*constraintsList).INVALID_TEMPLATE_QUALIFIER
      isType< Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface> >("Constraint Object"))
    constraints = (*constraintsList).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface> >("Constraint Object");
  else
    globalData->locaErrorCheck->throwError(methodName,
      "\"Constraint Object\" parameter is not of type Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>!");

  // Get parameter names for constraints
  if ((*constraintsList).INVALID_TEMPLATE_QUALIFIER
      isType< Teuchos::RCP< std::vector<std::string> > >("Constraint Parameter Names"))
    constraintParamNames = (*constraintsList).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RCP< std::vector<std::string> > >("Constraint Parameter Names");
  else
    globalData->locaErrorCheck->throwError(methodName,
      "\"Constraint Parameter Names\" parameter is not of type Teuchos::RCP< std::vector<std::string> >!");

  // Convert names to integer IDs
  std::vector<int> constraintParamIDs(constraintParamNames->size());
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
