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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER
#include "LOCA_MultiStepper.H"    // class definition
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_Combo.H"
#include "LOCA_StatusTest_Wrapper.H"

// LOCA Includes
#include "LOCA_Utils.H"                    // for static function doPrint
#include "LOCA_ErrorCheck.H"                    // for error checking methods
#include "LOCA_MultiContinuation_AbstractGroup.H"   // class data element
#include "LOCA_MultiContinuation_ArcLengthGroup.H"
#include "LOCA_MultiContinuation_NaturalGroup.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"
#include "MFLOCA.H"

// Multifario Includes
extern "C" {
#include <MFAtlas.h>
}

/* Some compilers (in particular the SGI and ASCI Red - TFLOP) 
 * fail to find the max and min function.  Therfore we redefine them 
 * here. 
 */ 
#ifdef max
#undef max
#endif
#define max(a,b) ((a)>(b)) ? (a) : (b);

#ifdef min
#undef min
#endif
#define min(a,b) ((a)<(b)) ? (a) : (b);


LOCA::MultiStepper::MultiStepper(
		       LOCA::MultiContinuation::AbstractGroup& initialGuess, 
		       NOX::StatusTest::Generic& t,
		       NOX::Parameter::List& p) :
  bifGroupManagerPtr(NULL),
  bifGroupPtr(NULL),
  curGroupPtr(NULL),
  statusTestPtr(NULL),
  paramListPtr(NULL),
  solverPtr(NULL),
  paramVec(),
  conParamIDVec(),
  conParamData()
{
  reset(initialGuess, t, p);
}

LOCA::MultiStepper::MultiStepper(const LOCA::MultiStepper& s) :
  bifGroupManagerPtr(NULL),
  bifGroupPtr(NULL),
  curGroupPtr(NULL),
  statusTestPtr(s.statusTestPtr),
  paramListPtr(s.paramListPtr),
  solverPtr(NULL),
  paramVec(s.paramVec),
  conParamIDVec(s.conParamIDVec),
  conParamData(s.conParamData)
{ 
  bifGroupManagerPtr =
    new LOCA::Bifurcation::Manager(*s.bifGroupManagerPtr);
  bifGroupPtr =
    dynamic_cast<LOCA::MultiContinuation::AbstractGroup*>(s.bifGroupPtr->clone());
  curGroupPtr = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedGroup*>(s.curGroupPtr->clone());

  // Right now this doesn't work because we can't copy the solver
}

LOCA::MultiStepper::~MultiStepper() 
{ 
  delete bifGroupManagerPtr;
  delete bifGroupPtr;
  delete curGroupPtr;
  delete solverPtr;
}

bool 
LOCA::MultiStepper::reset(LOCA::MultiContinuation::AbstractGroup& initialGuess,
			  NOX::StatusTest::Generic& t,
			  NOX::Parameter::List& p) 
{
  delete bifGroupPtr;
  delete curGroupPtr;
  delete bifGroupManagerPtr;
  delete solverPtr;

  paramListPtr = &p;
  statusTestPtr = &t;

  // Initialize the utilities
  LOCA::Utils::setUtils(*paramListPtr);

  // Get continuation parameter vector
  paramVec = initialGuess.getParams();

  // Get continuation parameter info
  getConParamData();

  // Reset group, predictor, step-size managers
  bifGroupManagerPtr =
    new LOCA::Bifurcation::Manager(LOCA::Utils::getSublist("Bifurcation"));

  // Create bifurcation group
  bifGroupPtr = dynamic_cast<LOCA::MultiContinuation::AbstractGroup*>
	        (bifGroupManagerPtr->createBifurcationGroup(initialGuess));

  // Create natural continuation group for first step
  curGroupPtr = new LOCA::MultiContinuation::NaturalGroup(
					 *bifGroupPtr, 
					 conParamIDVec, 
				         LOCA::Utils::getSublist("Stepper"));

  // Set step size			    
  for (unsigned int i=0; i<conParamIDVec.size(); i++)
    curGroupPtr->setStepSize(0.0, i);
  
  // Set previous solution vector in current solution group
  curGroupPtr->setPrevX(curGroupPtr->getX());

  // Create solver using initial conditions
  solverPtr = new NOX::Solver::Manager(*curGroupPtr, *statusTestPtr, 
				       LOCA::Utils::getSublist("NOX"));

  printInitializationInfo();

  if (LOCA::Utils::doPrint(LOCA::Utils::Parameters))
    paramListPtr->print(cout);
  
  return true;
}

LOCA::Abstract::Iterator::IteratorStatus
LOCA::MultiStepper::run() {
  NOX::StatusTest::StatusType solverStatus;
  string callingFunction = "LOCA::MultiStepper::run()";

  // Perform solve of initial conditions
  solverStatus = solverPtr->solve();

  // Set up arclength continuation group
  const LOCA::MultiContinuation::ExtendedGroup& constSolnGrp = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedGroup&>(solverPtr->getSolutionGroup());

  LOCA::Continuation::AbstractGroup& solnAbstractGrp = 
    const_cast<LOCA::Continuation::AbstractGroup&>(constSolnGrp.getUnderlyingGroup());
  LOCA::MultiContinuation::AbstractGroup& solnGrp =
    dynamic_cast<LOCA::MultiContinuation::AbstractGroup&>(solnAbstractGrp);
  delete curGroupPtr;
  curGroupPtr = new LOCA::MultiContinuation::ArcLengthGroup(
					solnGrp, 
					conParamIDVec, 
					LOCA::Utils::getSublist("Stepper"));

  // If nonlinear solve failed, return (this must be done after continuation 
  // groups are created so MultiStepper::getSolutionGroup() functions 
  // correctly.
  if (solverStatus != NOX::StatusTest::Converged)
    return LOCA::Abstract::Iterator::Failed;
  
  curGroupPtr->printSolution();

  delete solverPtr;
  solverPtr = new NOX::Solver::Manager(*curGroupPtr, *statusTestPtr, 
				       LOCA::Utils::getSublist("NOX"));

  // Get stepper sublist
  NOX::Parameter::List& stepperList = LOCA::Utils::getSublist("Stepper");

  MFImplicitMF M;
  MFNRegion Omega;
  MFAtlas A;
  MFNVector u0;
  MFContinuationMethod H;
  LOCAData* data = new LOCAData(*solverPtr, *curGroupPtr, *paramListPtr, 
				*statusTestPtr, conParamData);
  M=MFIMFCreateLOCA(data);
  Omega=MFNRegionCreateLOCA(data);
  u0=MFCreateLOCANVectorWithData(dynamic_cast<LMCEV *>(curGroupPtr->getX().clone()));

  H=MFCreateHendersonsMethod();

  /* Max distance from TS to M */
  MFHendersonsMethodSetEpsilon(H, stepperList.getParameter("Epsilon", 0.1));

  /* -1 means infinite */
  MFHendersonsMethodSetMaxCharts(H,stepperList.getParameter("Max Charts", -1));

  /* Write info to stdout */
  MFHendersonsMethodSetVerbose(H, stepperList.getParameter("Verbosity", 1));

  /* Page out non-interior polyhedra */
  MFHendersonsMethodSetPage(H, stepperList.getParameter("Page Charts", 1));

  /* Write polyhedra to a plotfile */
  MFHendersonsMethodSetDumpToPlotFile(H, stepperList.getParameter("Dump Polyhedra", true));

  /* Write points to a file */
  MFHendersonsMethodSetDumpToCenterFile(H,stepperList.getParameter("Dump Centers", false)); 

  /* File name to save data */
  const char *fname = 
    stepperList.getParameter("Filename", "MFresults").c_str();
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

LOCA::MultiContinuation::AbstractGroup& 
LOCA::MultiStepper::getSolutionGroup()
{
  return dynamic_cast<LOCA::MultiContinuation::AbstractGroup&>(curGroupPtr->getUnderlyingGroup());
}

const NOX::Parameter::List& 
LOCA::MultiStepper::getParameterList() const
{
  return *paramListPtr;
}

void 
LOCA::MultiStepper::printInitializationInfo()
{  
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperIteration)) {
    cout << endl << LOCA::Utils::fill(72, '~') << endl;
    cout << "Beginning Continuation Run \n" 
	 << "Method: MultiParmeter Continuation \n"
	 << endl;
    cout << LOCA::Utils::fill(72, '~') << endl << endl;
  }
}

void 
LOCA::MultiStepper::getConParamData()
{
  string callingFunction = "LOCA::MultiStepper::getConParamInfo()";

  // Get stepper sublist
  NOX::Parameter::List& stepperList = LOCA::Utils::getSublist("Stepper");

  // Get number of continuation parameters
  int numParams = 
    stepperList.getParameter("Number of Continuation Parameters",1);

  // Get data for each continuation parameter
  conParamData.clear();
  conParamIDVec.clear();
  for (int i=1; i<=numParams; i++) {

    // create sublist name
    stringstream sublistStream;
    sublistStream << "Continuation Parameter " << i;
    string sublistName = sublistStream.str();

    // Get sublist for continuation parameter
    if (!stepperList.isParameterSublist(sublistName)) {
      stringstream errorStream;
      errorStream << "No sublist for continuation parameter " << i << "!";
      LOCA::ErrorCheck::throwError(callingFunction, errorStream.str());
    }
    NOX::Parameter::List* conSublistPtr = &(stepperList.sublist(sublistName));

    // Create new struct for con param info
    ParamData d;

    // Get continuation parameter name
    if (!conSublistPtr->isParameter("Parameter Name")) {
      stringstream errorStream;
      errorStream << "\"Parameter Name\" for parameter sublist " 
                  << i << " is not set!";
      LOCA::ErrorCheck::throwError(callingFunction, errorStream.str());
    }
    d.name = conSublistPtr->getParameter("Parameter Name", "None");
    d.ID = paramVec.getIndex(d.name);

    // Get initial value
    if (!conSublistPtr->isParameter("Initial Value")) {
      stringstream errorStream;
      errorStream << "\"Initial Value\" for parameter sublist " 
                  << i << " is not set!";
      LOCA::ErrorCheck::throwError(callingFunction, errorStream.str());
    }
    d.initialValue = conSublistPtr->getParameter("Initial Value", 0.0);

    // Get max value
    if (!conSublistPtr->isParameter("Max Value")) {
      stringstream errorStream;
      errorStream << "\"Max Value\" for parameter sublist " 
                  << i << " is not set!";
      LOCA::ErrorCheck::throwError(callingFunction, errorStream.str());
    }
    d.maxValue = conSublistPtr->getParameter("Max Value", 0.0);

    // Get min value
    if (!conSublistPtr->isParameter("Min Value")) {
      stringstream errorStream;
      errorStream << "\"Min Value\" for parameter sublist " 
                  << i << " is not set!";
      LOCA::ErrorCheck::throwError(callingFunction, errorStream.str());
    }
    d.minValue = conSublistPtr->getParameter("Min Value", 0.0);

    // Get initial step size
    if (!conSublistPtr->isParameter("Initial Step Size")) {
      stringstream errorStream;
      errorStream << "\"Initial Step Size\" for parameter sublist " 
                  << i << " is not set!";
      LOCA::ErrorCheck::throwError(callingFunction, errorStream.str());
    }
    d.initialStepSize = conSublistPtr->getParameter("Initial Step Size", 0.1);

    // Get initial max size
    if (!conSublistPtr->isParameter("Max Step Size")) {
      stringstream errorStream;
      errorStream << "\"Max Step Size\" for parameter sublist " 
                  << i << " is not set!";
      LOCA::ErrorCheck::throwError(callingFunction, errorStream.str());
    }
    d.maxStepSize = conSublistPtr->getParameter("Max Step Size", 1.0);

    // Get initial min size
    if (!conSublistPtr->isParameter("Min Step Size")) {
      stringstream errorStream;
      errorStream << "\"Min Step Size\" for parameter sublist " 
                  << i << " is not set!";
      LOCA::ErrorCheck::throwError(callingFunction, errorStream.str());
    }
    d.minStepSize = conSublistPtr->getParameter("Min Step Size", 1.0e-3);

    // add data struct to list
    conParamData.push_back(d);
    conParamIDVec.push_back(d.ID);
  }
}
