// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// HMX Cook-off Test Problem

/* Solves the nonlinear equation:
 *
 */

// NOX Objects
#include "NOX.H"
#include "NOX_Epetra.H"
#include "NOX_TestCompare.H"
#include "NOX_Epetra_DebugTools.H"

// For parsing command line
#include "Teuchos_CommandLineProcessor.hpp"

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

#ifdef HAVE_NOX_ML_EPETRA
#include "Teuchos_ParameterList.hpp"
#endif

// Headers needed for FD coloring
#include <vector>
#ifdef HAVE_NOX_EPETRAEXT       // Use epetraext package in Trilinos
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#endif

// New coupling library headers
#include "NOX_Multiphysics_Solver_Manager.H"

// User's application specific files
#include "Problem_Manager.H"
#include "Problem_Interface.H"
#include "HMX_PDE.H"

// Added to allow timings
#include "Epetra_Time.h"

using namespace std;

int main(int argc, char *argv[])
{

  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  Teuchos::CommandLineProcessor clp( false );

  // Default run-time options that can be changed from the command line
  bool          verbose         = true          ;
  int           NumGlobalNodes  = 10    ;
  bool          runMF           = true  ;
  bool          useMatlab       = false ;
  bool          doOffBlocks     = true  ;
  bool          libloose        = true  ;
  std::string        outputDir       = "."   ;
  std::string        goldDir         = "."   ;

  clp.setOption( "verbose", "no-verbose", &verbose, "Verbosity on or off." );
  clp.setOption( "n", &NumGlobalNodes, "Number of elements" );
  clp.setOption( "runMF", "loose", &runMF, "Use Matrix-Free strong coupling" );
  clp.setOption( "offblocks", "no-offblocks", &doOffBlocks, "Include off-diagonal blocks in preconditioning matrix" );
  clp.setOption( "matlab", "no-matlab", &useMatlab, "Use Matlab debugging engine" );
  clp.setOption( "noxlib", "no-noxlib", &libloose, "Perform loose coupling using NOX's library (as opposed to hard-coded test driver)." );
  clp.setOption( "outputdir", &outputDir, "Directory to output mesh and results into. Default is \"./\"" );
  clp.setOption( "golddir", &goldDir, "Directory to read gold test from. Default is \"./\"" );

  Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);

  if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL )
    return parse_return;

  outputDir += "/";
  goldDir   += "/";

  // Create and reset the Timer
  Epetra_Time myTimer(Comm);
  double startWallTime = myTimer.WallTime();

  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  NumGlobalNodes++; // convert #elements to #nodes

  // The number of unknowns must be at least equal to the number of processors.
  if (NumGlobalNodes < NumProc)
  {
    std::cout << "numGlobalNodes = " << NumGlobalNodes
     << " cannot be < number of processors = " << NumProc << std::endl;
    exit(1);
  }

  // Begin Nonlinear Solver ************************************

  // NOTE: For now these parameters apply to all problems handled by
  // Problem_Manager.  Each problem could be made to have its own
  // parameter list as wwell as its own convergence test(s).

  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());

  // Set the nonlinear solver method
  nlParams.set("Nonlinear Solver", "Line Search Based");
  //nlParams.set("Nonlinear Solver", "Trust Region Based");

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", MyPID);
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);
  printParams.set("Output Information",
            NOX::Utils::Warning                  +
            NOX::Utils::OuterIteration           +
            NOX::Utils::InnerIteration           +
            NOX::Utils::Parameters               +
            NOX::Utils::Details                  +
            NOX::Utils::OuterIterationStatusTest +
            NOX::Utils::LinearSolverDetails      +
            NOX::Utils::TestDetails               );

  NOX::Utils outputUtils(printParams);

  // Sublist for line search
  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  searchParams.set("Method", "Full Step");

  // Sublist for direction
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  dirParams.set("Method", "Newton");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");

  // Sublist for linear solver for the Newton method
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");
  lsParams.set("Max Iterations", 800);
  lsParams.set("Tolerance", 1e-4);
  lsParams.set("Output Frequency", 50);
  lsParams.set("Preconditioner", "AztecOO");

  // Create the convergence tests
  // Note: as for the parameter list, both (all) problems use the same
  // convergence test(s) for now, but each could have its own.
  Teuchos::RCP<NOX::StatusTest::NormF>          absresid  =
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RCP<NOX::StatusTest::NormUpdate>     update    =
      Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-5));
  Teuchos::RCP<NOX::StatusTest::Combo>          converged =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> finiteValue =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);
  combo->addStatusTest(finiteValue);

  // Create the Problem Manager
  Problem_Manager problemManager(Comm, doOffBlocks, 0, useMatlab);

  // Note that each problem could contain its own nlParams list as well as
  // its own convergence test(s).
  problemManager.registerParameters(nlParamsPtr);
  problemManager.registerStatusTest(combo);

  std::string nameT         = "Temperature";
  double Const_R        = 1.9872 ;
  double Specific_H        = 0.42 ;
  double Density        = 1.90 ;
  double Thermal_K        = 0.8658e-3 ;
  double diffCoef_T        = Thermal_K / (Density * Specific_H);

  double StericCoef_T        = 0.0;
  double PreExp_T         = 0.0;
  double ActEnergy_T    = 0.0 ;
  map<string, double> SrcTermExponent_T; // Leave empty if no volume source
  map<string, double> SrcTermWeight_T;
    SrcTermWeight_T.insert( pair<string, double> ("SpeciesA", -190.0) );
    SrcTermWeight_T.insert( pair<string, double> ("SpeciesB",  570.0) );
    SrcTermWeight_T.insert( pair<string, double> ("SpeciesC", 2280.0) );

  // Create each part of the HMX cook-off problem
  HMX_PDE HMX_TempEq (Comm,
                  diffCoef_T,
                  Const_R,
                  StericCoef_T, // Dummy for Temp Eq.
                  PreExp_T, // Dummy for Temp Eq.
                  ActEnergy_T, // Dummy for Temp Eq.
                  SrcTermExponent_T,
                  SrcTermWeight_T,
                  1*NumGlobalNodes, nameT);

  HMX_TempEq.setTempFieldName(HMX_TempEq.getName());

  // Override default initialization with values we want
  HMX_TempEq.initializeSolution(500.0);

  problemManager.addProblem(HMX_TempEq);


  std::string nameA         = "SpeciesA";
  double diffCoef_A         = 0.0 ;
  //double stericCoef_A        = 0.0 ;
  double stericCoef_A        = 2.0 ; // ROGER: high coupling
  //double stericCoef_A        = 2.0 ; // ROGER: high coupling
  double preExp_A         = exp(48.7) ;
  double actEnergy_A         = 52700.0 ;
  map<string, double> SrcTermExponent_A;
    SrcTermExponent_A.insert( pair<string, double> (nameA, 1.0) );
  map<string, double> SrcTermWeight_A;
    SrcTermWeight_A.insert( pair<string, double> (nameA, -1.0) );

  HMX_PDE HMX_RxnA(Comm,
                diffCoef_A,
                Const_R,
                stericCoef_A,
                preExp_A,
                actEnergy_A,
                SrcTermExponent_A,
                SrcTermWeight_A,
                1*NumGlobalNodes, nameA);

  HMX_RxnA.setTempFieldName(HMX_TempEq.getName());

  // Override default initialization with values we want
  HMX_RxnA.initializeSolution(2.0);

  problemManager.addProblem(HMX_RxnA);


  std::string nameB         = "SpeciesB";
  double diffCoef_B         = 0.0 ;
  double stericCoef_B        = 1.0 ;
  //double stericCoef_B        = -1.0 ;
  double preExp_B         = exp(37.5) ;
  double actEnergy_B         = 44100.0 ;
  map<string, double> SrcTermExponent_B;
    SrcTermExponent_B.insert( pair<string, double> (nameB, 1.0) );
  map<string, double> SrcTermWeight_B;
    SrcTermWeight_B.insert( pair<string, double> (nameA, 1.0) );
    SrcTermWeight_B.insert( pair<string, double> (nameB, -1.0) );

  HMX_PDE HMX_RxnB(Comm,
                diffCoef_B,
                Const_R,
                stericCoef_B,
                preExp_B,
                actEnergy_B,
                SrcTermExponent_B,
                SrcTermWeight_B,
                NumGlobalNodes, nameB);

  HMX_RxnB.setTempFieldName(HMX_TempEq.getName());

  // Override default initialization with values we want
  HMX_RxnB.initializeSolution(1.0);

  problemManager.addProblem(HMX_RxnB);


  std::string nameC         = "SpeciesC";
  double diffCoef_C         = 0.0 ;
  double stericCoef_C        = 0.0 ;
  double preExp_C         = exp(28.1) ;
  double actEnergy_C         = 34100.0 ;
  map<string, double> SrcTermExponent_C;
    SrcTermExponent_C.insert( pair<string, double> (nameC, 2.0) );
  map<string, double> SrcTermWeight_C;
    SrcTermWeight_C.insert( pair<string, double> (nameB, 2.0) );
    SrcTermWeight_C.insert( pair<string, double> (nameC, -2.0) );

  HMX_PDE HMX_RxnC(Comm,
                diffCoef_C,
                Const_R,
                stericCoef_C,
                preExp_C,
                actEnergy_C,
                SrcTermExponent_C,
                SrcTermWeight_C,
                1*NumGlobalNodes, nameC);

  HMX_RxnC.setTempFieldName(HMX_TempEq.getName());

  // Override default initialization with values we want
//    HMX_RxnC.initializeSolution(0.0);

  problemManager.addProblem(HMX_RxnC);

  problemManager.createDependency(HMX_TempEq, HMX_RxnA);
  problemManager.createDependency(HMX_TempEq, HMX_RxnB);
  problemManager.createDependency(HMX_TempEq, HMX_RxnC);

  problemManager.createDependency(HMX_RxnA, HMX_TempEq);

  problemManager.createDependency(HMX_RxnB, HMX_TempEq);
  problemManager.createDependency(HMX_RxnB, HMX_RxnA);

  problemManager.createDependency(HMX_RxnC, HMX_TempEq);
  problemManager.createDependency(HMX_RxnC, HMX_RxnB);

  // Initialize time integration parameters
  int maxTimeSteps = 3;
  int timeStep = 0;
  double time = 0.;
  //double dt = HMX_TempEq.getdt();
  double dt = 10.0 * HMX_TempEq.getdt();

  problemManager.registerComplete();

  problemManager.outputStatus(std::cout);

  // Print initial solution
  if( verbose )
    problemManager.outputSolutions( outputDir);

  // Identify the test problem
  if( outputUtils.isPrintType(NOX::Utils::TestDetails) )
    outputUtils.out() << "Starting epetra/MultiPhysics/example_hmx.exe" << std::endl;

  // Identify processor information
#ifdef HAVE_MPI
  outputUtils.out() << "This test is broken in parallel." << std::endl;
  outputUtils.out() << "Test failed!" << std::endl;
  MPI_Finalize();
  return -1;
#else
  if (outputUtils.isPrintType(NOX::Utils::TestDetails))
    outputUtils.out() << "Serial Run" << std::endl;
#endif

  // Create a TestCompare class
  int status = 0;
  NOX::TestCompare tester( outputUtils.out(), outputUtils );
  double abstol = 1.e-4;
  double reltol = 1.e-4 ;

  // Time integration loop
  while(timeStep < maxTimeSteps) {

    timeStep++;
    time += dt;

    std::cout << "Time Step: " << timeStep << ",\tTime: " << time << std::endl;

    // Solve the coupled problem
    if( runMF )
      problemManager.solveMF(); // Need a status test check here ....
    else if( !libloose )
      problemManager.solve(); // Hard-coded loose coupling
    else // Loose coupling via NOX library
    {
      // Create the loose coupling solver manager
      Teuchos::RCP<vector<Teuchos::RCP<NOX::Solver::Generic> > > solversVec =
        Teuchos::rcp( new std::vector<Teuchos::RCP<NOX::Solver::Generic> > );

      map<int, Teuchos::RCP<NOX::Solver::Generic> >::iterator iter = problemManager.getSolvers().begin(),
                                                                  iter_end = problemManager.getSolvers().end()   ;
      for( ; iter_end != iter; ++iter )
      {
        std::cout << " ........  registered Solver::Manager # " << (*iter).first << std::endl;
        solversVec->push_back( (*iter).second );
      }

      // Package the Problem_Manager as the DataExchange::Intreface
      Teuchos::RCP<NOX::Multiphysics::DataExchange::Interface> dataExInterface =
        Teuchos::rcp( &problemManager, false );

      NOX::Multiphysics::Solver::Manager cplSolv( solversVec, dataExInterface, combo, nlParamsPtr );

      cplSolv.solve();

      // Refresh all problems with solutions from solver
      problemManager.copyAllGroupXtoProblems();

      // Reset all solver groups to force recomputation of residuals
      problemManager.resetAllCurrentGroupX();
    }

  if( verbose )
    problemManager.outputSolutions( outputDir, timeStep );

    map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator iter = problemManager.getProblems().begin(),
                                                                iter_end = problemManager.getProblems().end()   ;
    for( ; iter_end != iter; ++iter )
    {
      GenericEpetraProblem & problem = *(*iter).second;
      std::string msg = "Numerical-to-Gold Solution comparison for problem \"" + problem.getName() + "\"";

      // Get the gold copy to comapre against current solution
      std::string baseFileame = problemManager.createIOname( problem, timeStep );
      std::string goldFileame = goldDir + "gold_hmx/" + baseFileame;

      Epetra_Vector * tmpVec = NULL;
      int ierr = NOX::Epetra::DebugTools::readVector( goldFileame, Comm, tmpVec );
      if( ierr != 0 )
      {
        outputUtils.out() << "ERROR opening gold copy file \"" << goldFileame << "\"." << std::endl;
        status = ierr;
        break;
      }

      Teuchos::RCP<Epetra_Vector> goldSoln = Teuchos::rcp( tmpVec );

      // Need NOX::Epetra::Vectors for tests
      NOX::Epetra::Vector numerical ( problem.getSolution() , NOX::Epetra::Vector::CreateView );
      NOX::Epetra::Vector goldVec   ( goldSoln       , NOX::Epetra::Vector::CreateView );

      status += tester.testVector( numerical, goldVec, reltol, abstol, msg );

      // Quit if test fails
      if( status != 0 )
        break;
    }

    // Quit if test fails
    if( status != 0 )
      break;

    // Reset problems by copying solution into old solution
    problemManager.resetProblems();

  } // end time step while loop

  // Output timing info
  if(MyPID==0)
    std::cout << "\nTimings :\n\tWallTime --> " << myTimer.WallTime() - startWallTime << " sec."
         << "\n\tElapsedTime --> " << myTimer.ElapsedTime() << " sec." << std::endl << std::endl;

  if( 1 ) // this will be turned on later
  {
    // Summarize test results
    if( status == 0 )
      outputUtils.out() << "Test passed!" << std::endl;
    else
      outputUtils.out() << "Test failed!" << std::endl;
  }
  else // force this test to pass for now, but at least warn of failure
  {
    // Summarize test results
    if( status == 0 )
      outputUtils.out() << "Test passed!" << std::endl;
    else
    {
      outputUtils.out() << "This test actually F-A-I-L-E-D." << std::endl;
      outputUtils.out() << "Test passed!" << std::endl;
    }
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return 0 ;
}

