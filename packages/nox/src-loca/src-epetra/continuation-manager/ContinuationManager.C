/*
//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
*/

#include "ContinuationManager.H"

// Trilinos headers
#include "LOCA_GlobalData.H"

// DefCont headers
#include "LOCAInterface.H"
//#include "PhaseConstraint.h"
//#include "IOVtkUtils.h"
#include "IOContFileUtils.H"
#ifdef HAVE_NOX_AMESOS
#include "NOX_Epetra_LinearSystem_Amesos.H"
#endif

ContinuationManager::
ContinuationManager( 
    //const Teuchos::RCP< Epetra_MpiComm > & aComm ,
    const Teuchos::RCP< Epetra_Comm > & aComm ,
    const string & taskFileName):
  comm(aComm),
  continuationFileName("continuation.dat"),
  iniStepLabel(0),
  interfaceConstraint(Teuchos::null),
  //isConstrainedProblem(false),
  locaGlobalData(Teuchos::null),
  locaStepper(Teuchos::null),
  locaStepperStatus(LOCA::Abstract::Iterator::NotFinished),
  maxAllowedSteps(100000),
  outputDir("."),
  problem(Teuchos::null),
  solutionFilesExtension("vtk"),
  solutionFilesPrefix("step_"),
  taskList(Teuchos::rcp (new Teuchos::ParameterList())),
  timeCounter(Epetra_Time(*aComm))
{
  // Reading the parameter task from the file
  Teuchos::updateParametersFromXmlFile(taskFileName, taskList.ptr());

  if (comm->MyPID()==0) {
    std::cout << std::endl << "#### Task Parameters from task file \"" << 
      taskFileName << " ####" << std::endl;
    taskList->print(std::cout,2,false,false);
    std::cout << std::endl << "#### End Parameters from "<< taskFileName << 
      " ####" << std::endl << std::endl;
  }

  // Getting the initial step label flag
  iniStepLabel = taskList->sublist("Continuation Manager").
                         sublist("Continuation").
			   get<int>("Label From Step");

//  // If the constraint sublist exists, check whether the constraint is enabled
//  if ( taskList->sublist("Continuation Manager").isSublist("Constraint") )
//    isConstrainedProblem = taskList->sublist("Continuation Manager").
//                                        sublist("Constraint").
//					 get<bool>("Enable Constraint");

  // Getting solution files prefix and extensions
  if ( taskList->sublist("Continuation Manager").sublist("Continuation").
                                                 isParameter("Solution Files Prefix") )
    solutionFilesPrefix = taskList->sublist("Continuation Manager").
                                     sublist("Continuation").
				     get<string>("Solution Files Prefix");

  if ( taskList->sublist("Continuation Manager").sublist("Continuation").
                                                 isParameter("Solution Files Extension") )
    solutionFilesExtension = taskList->sublist("Continuation Manager").
                                     sublist("Continuation").
				     get<string>("Solution Files Extension");

}

ContinuationManager::
~ContinuationManager()
{
  destroyGlobalData(locaGlobalData);
}

bool ContinuationManager::
BuildLOCAStepper()
{

  if (comm->MyPID()==0) std::cout << std::endl << "Building the LOCA stepper..." << std::endl;

  // Make sure the problem has been set
  TEUCHOS_TEST_FOR_EXCEPTION( problem == Teuchos::null, 
      std::logic_error,
      "ContinuationManager has not been given a valid ProblemLOCAPrototype");

  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base 
  // class
  Teuchos::RCP <LOCAInterface> interface = 
    Teuchos::rcp(new LOCAInterface(problem,Teuchos::rcp(&*this,false)));
  Teuchos::RCP <LOCA::Epetra::Interface::Required> interfaceRequired = interface;
  Teuchos::RCP <NOX::Epetra::Interface::Jacobian> interfaceJacobian = interface;

  // Create the Epetra_RowMatrix for the Jacobian/Preconditioner
  Teuchos::RCP <Epetra_RowMatrix> jacobian = problem->GetJacF();

  // Get the initial guess vector
  Teuchos::RCP <Epetra_Vector> initialGuess = 
    problem->GetInitialGuess();

  // NOX and LOCA sublist
  Teuchos::RCP <Teuchos::ParameterList> noxAndLocaList =
    Teuchos::rcp(new Teuchos::ParameterList(taskList->sublist("NOX and LOCA")));

  // Create the lists needed for linear systems
  Teuchos::ParameterList & noxPrintingList = 
    noxAndLocaList->
      sublist("NOX").
        sublist("Printing");
  Teuchos::ParameterList & linearSystemList = 
    noxAndLocaList->
      sublist("NOX").
        sublist("Direction").
          sublist("Newton").
            sublist("Linear Solver");

  // Instntiating the appropriate linear system
  string linearSolverType = linearSystemList.get("Solver","Aztec");
  Teuchos::RCP<NOX::Epetra::LinearSystem> linearSystem = Teuchos::null; 

  if (linearSolverType == "Aztec")
    // Create the linear system
    //Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linearSystem = 
    linearSystem = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(noxPrintingList, 
	    linearSystemList,
	    interfaceRequired,
	    interfaceJacobian,
	    jacobian,
	    *initialGuess));
#ifdef HAVE_NOX_AMESOS
  else if (linearSolverType == "Amesos")
    // Create the linear system
    //Teuchos::RCP<NOX::Epetra::LinearSystemAmesos> linearSystem = 
    linearSystem = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemAmesos(noxPrintingList, 
	    linearSystemList,
	    interfaceRequired,
	    interfaceJacobian,
	    jacobian,
	    *initialGuess));
#endif
  else
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, 
	std::logic_error, "Continuation manager: " + 
			  linearSolverType + " is not a valid linear solver");
  }

  // Get the parameter vector
  LOCA::ParameterVector continuableParams = problem->GetContinuableParams();

  // Create Epetra factory
  Teuchos::RCP <LOCA::Abstract::Factory> epetraFactory =
    Teuchos::rcp(new LOCA::Epetra::Factory);

  // Create the loca vector
  Teuchos::RCP <NOX::Epetra::Vector> locaInitialGuess = 
    Teuchos::rcp (new NOX::Epetra::Vector(*initialGuess));

  // Instantiate the constraint objects
  //if (isConstrainedProblem) {
  if (interfaceConstraint != Teuchos::null) {

//    // Instantiate the constraint
//    Teuchos::RCP <PhaseConstraint> phaseConstraint = 
//      Teuchos::rcp(new PhaseConstraint(problem,*locaInitialGuess));
//
//    // Instantiate the interface
//    Teuchos::RCP <LOCA::MultiContinuation::ConstraintInterface> 
//      interfaceConstraint = phaseConstraint;

    // Instantiate the constraint parameters names
    Teuchos::RCP< vector<string> > constraintParamsNames = 
      Teuchos::rcp(new vector<string>());

    // The user-defined constrained parameters
    Teuchos::ParameterList & constraintParams = 
      taskList->sublist("Continuation Manager").
                   sublist("Constraint").
		    sublist("Constraint Parameters");

    // Getting the parameter names from the user-defined list
    Teuchos::ParameterList::ConstIterator i;
    for (i = constraintParams.begin(); i !=constraintParams.end(); ++i) 
      constraintParamsNames->push_back(
	  constraintParams.get<string>(constraintParams.name(i)));

    // Instantiating the constraint list
    Teuchos::ParameterList & constraintsList =
      noxAndLocaList->sublist("LOCA").sublist("Constraints");
    constraintsList.set("Constraint Object", interfaceConstraint);
    constraintsList.set("Constraint Parameter Names", 
	constraintParamsNames);
    constraintsList.set("Bordered Solver Method", "Householder");

  }

  // Create the global data object
  locaGlobalData = LOCA::createGlobalData(noxAndLocaList, epetraFactory);

  // Create the Group
  Teuchos::RCP<LOCA::Epetra::Group> group = 
    Teuchos::rcp(new LOCA::Epetra::Group(locaGlobalData, noxPrintingList, 
      interfaceRequired, *locaInitialGuess, linearSystem, continuableParams));

  // Create the Solver convergence test
  Teuchos::RCP<NOX::StatusTest::NormF> normTolerance = 
    Teuchos::rcp(new NOX::StatusTest::NormF(
	  taskList->
	    sublist("Continuation Manager").
	      sublist("Continuation").
		get<double>("Nonlinear Step Tolerance"),
	  NOX::StatusTest::NormF::Scaled));
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxIterations = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(
	  noxAndLocaList->
	    sublist("LOCA").
	      sublist("Stepper").
	        get<int>("Max Nonlinear Iterations")));
  Teuchos::RCP<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(normTolerance);
  combo->addStatusTest(maxIterations);

  // Create the stepper
  locaStepper = Teuchos::rcp (new LOCA::Stepper(locaGlobalData, group, combo, noxAndLocaList));

  // Performing some checks
  ValidateLOCAStepper();

  return true;
}


bool ContinuationManager::
BuildLOCAPeriodicStepper(const Teuchos::RCP<EpetraExt::MultiComm> globalComm)
{

  if (comm->MyPID()==0) cout << endl << "Building the LOCA stepper..." << endl;

  // Make sure the problem has been set
  TEUCHOS_TEST_FOR_EXCEPTION( problem == Teuchos::null, 
      std::logic_error,
      "ContinuationManager has not been given a valid ProblemLOCAPrototype");

  // Create the Epetra_RowMatrix for the Jacobian/Preconditioner
  Teuchos::RCP <Epetra_RowMatrix> jacobian = problem->GetJacF();

  // Get the initial guess vector
  Teuchos::RCP <Epetra_Vector> initialGuess = 
    problem->GetInitialGuess();

  // NOX and LOCA sublist
  Teuchos::RCP <Teuchos::ParameterList> noxAndLocaList =
    Teuchos::rcp(new Teuchos::ParameterList(taskList->sublist("NOX and LOCA")));

  // Create the lists needed for linear systems
  Teuchos::ParameterList & noxPrintingList = 
    noxAndLocaList->
      sublist("NOX").
        sublist("Printing");
  Teuchos::ParameterList & linearSystemList = 
    noxAndLocaList->
      sublist("NOX").
        sublist("Direction").
          sublist("Newton").
            sublist("Linear Solver");

  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base 
  // class
  Teuchos::RCP <LOCAInterface> interface = 
    Teuchos::rcp(new LOCAInterface(problem,Teuchos::rcp(&*this,false)));

  Epetra_MultiVector guessMV(initialGuess->Map(), globalComm->NumTimeStepsOnDomain());
  for (int i=0; i<globalComm->NumTimeStepsOnDomain(); i++) *(guessMV(i)) = *initialGuess;


cout << "XXX  num time steps on domain = " <<  globalComm->NumTimeStepsOnDomain() << endl;

  double dt = 1.0;
  Teuchos::RCP <LOCA::Epetra::Interface::xyzt> xyzt_interface = 
   Teuchos::rcp(new LOCA::Epetra::Interface::xyzt(interface, guessMV, jacobian, globalComm, *initialGuess, dt ));


  Teuchos::RCP <LOCA::Epetra::Interface::Required> interfaceRequired = xyzt_interface;
  Teuchos::RCP <NOX::Epetra::Interface::Jacobian> interfaceJacobian = xyzt_interface;

  Teuchos::RCP<Epetra_RowMatrix> Axyzt =
     Teuchos::rcp(&(xyzt_interface->getJacobian()),false);
  Epetra_Vector& solnxyzt = xyzt_interface->getSolution();


  // Create the linear system
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linearSystem = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(noxPrintingList, 
          linearSystemList,
          interfaceRequired,
          interfaceJacobian,
          Axyzt,
          solnxyzt));

  // Get the parameter vector
  LOCA::ParameterVector continuableParams = problem->GetContinuableParams();

  // Create Epetra factory
  Teuchos::RCP <LOCA::Abstract::Factory> epetraFactory =
    Teuchos::rcp(new LOCA::Epetra::Factory);

  // Create the loca vector
  Teuchos::RCP <NOX::Epetra::Vector> locaInitialGuess = 
    Teuchos::rcp (new NOX::Epetra::Vector(solnxyzt));

  // Instantiate the constraint objects
  //if (isConstrainedProblem) 
  if (interfaceConstraint != Teuchos::null) {

//    // Instantiate the constraint
//    Teuchos::RCP <PhaseConstraint> phaseConstraint = 
//      Teuchos::rcp(new PhaseConstraint(problem,*locaInitialGuess));
//
//    // Instantiate the interface
//    Teuchos::RCP <LOCA::MultiContinuation::ConstraintInterface> 
//      interfaceConstraint = phaseConstraint;

    // Instantiate the constraint parameters names
    Teuchos::RCP< vector<string> > constraintParamsNames = 
      Teuchos::rcp(new vector<string>());

    // The user-defined constrained parameters
    Teuchos::ParameterList & constraintParams = 
      taskList->sublist("Continuation Manager").
                   sublist("Constraint").
		    sublist("Constraint Parameters");

    // Getting the parameter names from the user-defined list
    Teuchos::ParameterList::ConstIterator i;
    for (i = constraintParams.begin(); i !=constraintParams.end(); ++i) 
      constraintParamsNames->push_back(
	  constraintParams.get<string>(constraintParams.name(i)));

    // Instantiating the constraint list
    Teuchos::ParameterList & constraintsList =
      noxAndLocaList->sublist("LOCA").sublist("Constraints");
    constraintsList.set("Constraint Object", interfaceConstraint);
    constraintsList.set("Constraint Parameter Names", 
	constraintParamsNames);
    constraintsList.set("Bordered Solver Method", "Householder");

  }

  // Create the global data object
  locaGlobalData = LOCA::createGlobalData(noxAndLocaList, epetraFactory);

  // Create the Group
  Teuchos::RCP<LOCA::Epetra::Group> group = 
    Teuchos::rcp(new LOCA::Epetra::Group(locaGlobalData, noxPrintingList, 
      interfaceRequired, *locaInitialGuess, linearSystem, continuableParams));

  // Create the Solver convergence test
  Teuchos::RCP<NOX::StatusTest::NormF> normTolerance = 
    Teuchos::rcp(new NOX::StatusTest::NormF(
	  taskList->
	    sublist("Continuation Manager").
	      sublist("Continuation").
		get<double>("Nonlinear Step Tolerance"),
	  NOX::StatusTest::NormF::Scaled));
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxIterations = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(
	  noxAndLocaList->
	    sublist("LOCA").
	      sublist("Stepper").
	        get<int>("Max Nonlinear Iterations")));
  Teuchos::RCP<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(normTolerance);
  combo->addStatusTest(maxIterations);

  // Create the stepper
  locaStepper = Teuchos::rcp (new LOCA::Stepper(locaGlobalData, group, combo, noxAndLocaList));

  // Performing some checks
  ValidateLOCAStepper();

  return true;
}


bool ContinuationManager::
ValidateLOCAStepper()
{
  // The number of steps must not exceed maxAllowedSteps
  int numSteps = iniStepLabel + locaStepper->getList()->sublist("LOCA").
                                                     sublist("Stepper").
				                       get<int>("Max Steps");
  ostringstream osstream;
  osstream << maxAllowedSteps;
  string errString = "iniStepLabel + Max number of continuation steps must not exceed " + 
    osstream.str();
  TEUCHOS_TEST_FOR_EXCEPTION( numSteps > maxAllowedSteps, 
      std::logic_error,
      errString);

  return true;
}

bool ContinuationManager::
RunLOCAStepper()
{
  // Touch the continuation file parameter list to get rid of the annoying "[empty]" flag
  TouchContFileParameters( *(problem->GetContinuationFileParameters()) );

  // Initialise continuation file
  if ( (comm->MyPID() == 0) && (iniStepLabel == 0) )
    WriteHeaderToContFile( GetContinuationFileName(),
	*(problem->GetContinuationFileParameters()) );

  // Running
  locaStepperStatus =  locaStepper->run();

  return true;
}

bool ContinuationManager::
PrintLOCAStepperStatistics()
{
  // Ckecking the convergence of a nonlinear step
  if (locaStepperStatus != LOCA::Abstract::Iterator::Finished) 
     if (locaGlobalData->locaUtils->isPrintType(NOX::Utils::Error))
       locaGlobalData->locaUtils->out() << "Stepper failed to converge!" << std::endl;

  // Output the parameter list
  if (locaGlobalData->locaUtils->isPrintType(NOX::Utils::StepperParameters)) {
    locaGlobalData->locaUtils->out() << std::endl << 
      "### Final Parameters ############" << std::endl;
    locaStepper->getList()->print(locaGlobalData->locaUtils->out());
    locaGlobalData->locaUtils->out() << std::endl;
  }

  // The time spent
  if (comm->MyPID() == 0)
    std::cout << std::endl << "#### Statistics ########" << std::endl;

  std::cout << " Time on proc " << comm->MyPID() << " = " 
       << timeCounter.ElapsedTime() << std::endl;

  // Check number of steps
  int numSteps = locaStepper->getStepNumber();
  if (comm->MyPID() == 0) 
    std::cout << " Number of continuation Steps = " << numSteps << std::endl;

  // Check number of failed steps
  int numFailedSteps = locaStepper->getNumFailedSteps();
  if (comm->MyPID() == 0) 
    std::cout << " Number of failed continuation Steps = " << numFailedSteps << std::endl;

  return true;
}

string ContinuationManager:: 
GetSolutionFileName() const
{
  // Number of digits
  int numDigits = 
    static_cast<int>( std::floor( std::log10( static_cast<double>(maxAllowedSteps) ) ) ) + 1;

  // Composing the filename
  ostringstream fileName;
  fileName << outputDir + "/" +
              solutionFilesPrefix +
	      StringifyInt(iniStepLabel + locaStepper->getStepNumber(), numDigits) +
	      "." + solutionFilesExtension;

  return fileName.str();
}

string ContinuationManager:: 
GetContinuationFileName() const
{
  return outputDir + "/" + continuationFileName;
}

int ContinuationManager:: 
GetStepID() const
{
  return (iniStepLabel+locaStepper->getStepNumber());
}

ContinuationManager::SolutionFileAttribute ContinuationManager::
GetSolutionFileAttribute() const
{
  // Default is NoPrint
  ContinuationManager::SolutionFileAttribute attribute = NoPrint;

  // Number of steps
  int numSteps =  locaStepper->getStepNumber();

  // Steps per print
  int stepsPerPrint = 1;

  if ( taskList->sublist("Continuation Manager").
                  sublist("Continuation").
		   isParameter("Steps Per Print"))
    stepsPerPrint = taskList->sublist("Continuation Manager").
                               sublist("Continuation").
			        get<int>("Steps Per Print");

  // Every few steps, do print
  if (  (stepsPerPrint != 0 ) && ( numSteps % stepsPerPrint  == 0 ) ) 
    attribute = Print;

  return attribute;
}

Teuchos::RCP <Teuchos::ParameterList> ContinuationManager::
GetTaskList() const
{
  return taskList;
}

bool ContinuationManager::
SetLOCAProblem( const Teuchos::RCP <ProblemLOCAPrototype> & aProblem )
{
  problem = aProblem;
  return true;
}

bool ContinuationManager::
SetLOCAConstraint( 
    const Teuchos::RCP <LOCA::MultiContinuation::ConstraintInterface> &
    anInterfaceConstraint)
{
  interfaceConstraint = anInterfaceConstraint;
  return true;
}

bool ContinuationManager::
SetTaskList( Teuchos::RCP <Teuchos::ParameterList> aTaskList )
{
  taskList = aTaskList;
  return true;
}

string ContinuationManager::
StringifyInt( const int & intNumber , const int & digits) const 
{
  // The number of Steps
  ostringstream osstream;
  osstream.width(digits);
  osstream.fill('0');
  osstream << right << intNumber;

  return osstream.str();
}
