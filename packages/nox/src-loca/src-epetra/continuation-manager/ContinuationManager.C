// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

#include <sstream>

ContinuationManager::
ContinuationManager(
    //const Teuchos::RCP< Epetra_MpiComm > & aComm ,
    const Teuchos::RCP< Epetra_Comm > & aComm ,
    const std::string & taskFileName):
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
//                     get<bool>("Enable Constraint");

  // Getting solution files prefix and extensions
  if ( taskList->sublist("Continuation Manager").sublist("Continuation").
                                                 isParameter("Solution Files Prefix") )
    solutionFilesPrefix = taskList->sublist("Continuation Manager").
                                     sublist("Continuation").
                     get<std::string>("Solution Files Prefix");

  if ( taskList->sublist("Continuation Manager").sublist("Continuation").
                                                 isParameter("Solution Files Extension") )
    solutionFilesExtension = taskList->sublist("Continuation Manager").
                                     sublist("Continuation").
                     get<std::string>("Solution Files Extension");

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
  std::string linearSolverType = linearSystemList.get("Solver","Aztec");
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
    Teuchos::RCP< std::vector<std::string> > constraintParamsNames =
      Teuchos::rcp(new std::vector<std::string>());

    // The user-defined constrained parameters
    Teuchos::ParameterList & constraintParams =
      taskList->sublist("Continuation Manager").
                   sublist("Constraint").
            sublist("Constraint Parameters");

    // Getting the parameter names from the user-defined list
    Teuchos::ParameterList::ConstIterator i;
    for (i = constraintParams.begin(); i !=constraintParams.end(); ++i)
      constraintParamsNames->push_back(
      constraintParams.get<std::string>(constraintParams.name(i)));

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

  if (comm->MyPID()==0) std::cout << std::endl << "Building the LOCA stepper..." << std::endl;

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


std::cout << "XXX  num time steps on domain = " <<  globalComm->NumTimeStepsOnDomain() << std::endl;

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
    Teuchos::RCP< std::vector<std::string> > constraintParamsNames =
      Teuchos::rcp(new std::vector<std::string>());

    // The user-defined constrained parameters
    Teuchos::ParameterList & constraintParams =
      taskList->sublist("Continuation Manager").
                   sublist("Constraint").
            sublist("Constraint Parameters");

    // Getting the parameter names from the user-defined list
    Teuchos::ParameterList::ConstIterator i;
    for (i = constraintParams.begin(); i !=constraintParams.end(); ++i)
      constraintParamsNames->push_back(
      constraintParams.get<std::string>(constraintParams.name(i)));

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
  std::ostringstream osstream;
  osstream << maxAllowedSteps;
  std::string errString = "iniStepLabel + Max number of continuation steps must not exceed " +
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

std::string ContinuationManager::
GetSolutionFileName() const
{
  // Number of digits
  int numDigits =
    static_cast<int>( std::floor( std::log10( static_cast<double>(maxAllowedSteps) ) ) ) + 1;

  // Composing the filename
  std::ostringstream fileName;
  fileName << outputDir + "/" +
              solutionFilesPrefix +
          StringifyInt(iniStepLabel + locaStepper->getStepNumber(), numDigits) +
          "." + solutionFilesExtension;

  return fileName.str();
}

std::string ContinuationManager::
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

std::string ContinuationManager::
StringifyInt( const int & intNumber , const int & digits) const
{
  // The number of Steps
  std::ostringstream osstream;
  osstream.width(digits);
  osstream.fill('0');
  osstream << std::right << intNumber;

  return osstream.str();
}
