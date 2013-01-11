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
                                                                                
// 1D Finite Element Brusselator Test Problem

/* Solves the nonlinear equation:
 *
 * dT       d2T    
 * --- - D1 --- - alpha + (beta+1)*T - C*T**2 = 0
 * dt       dx2   
 *
 * T(t,0) = T(t,1) = alpha = 0.6
 * T(0,x) = alpha + sinusoidal perturbation
 *
 *
 * dC       d2C    
 * --- - D2 --- - beta*T + C*T**2 = 0
 * dt       dx2   
 *
 * C(t,0) = C(t,1) = beta / alpha = 2.0 / 0.6
 * C(0,x) = beta / alpha + sinusoidal perturbation
 *
 * and
 *      D1 = D2 = 0.025
 *
 * with d representing partial differentiation.
 *
 * This problem is examined with a variety of time integration schemes in:
 * "Studies on the Convergence of Various Time-Integration Schemes for the
 * Radiation-Diffusion Problem," Curtis C. Ober & John N. Shadid, in prep.
 *
 * In this example, only a 1st-order fully implicit (backward Euler)
 * time integration scheme is considered currently.
 *
 * Values for time step size and finite spatial extent are specified in
 * the constructor initialization list in Brusselator.C using
 * variables dt  and xmin,xmax, respectively.
 * The number of time steps to be taken is specified by variable
 * maxTimeSteps below.
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
#include "Equation_A.H"              
#include "Equation_B.H"              
#include "Burgers.H"              

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
  bool          verbose         = true  ;
  int           NumGlobalNodes  = 20    ;
  bool          runMF           = true  ;
  bool          useMatlab       = false ;
  bool          doOffBlocks     = false ;
  bool          convection      = true  ;
  bool          libloose        = true  ;
  std::string        outputDir       = "."   ;
  std::string        goldDir         = "."   ;

  clp.setOption( "verbose", "no-verbose", &verbose, "Verbosity on or off." );
  clp.setOption( "n", &NumGlobalNodes, "Number of elements" );
  clp.setOption( "runMF", "loose", &runMF, "Use Matrix-Free strong coupling" );
  clp.setOption( "offblocks", "no-offblocks", &doOffBlocks, "Include off-diagonal blocks in preconditioning matrix" );
  clp.setOption( "burgers", "no-burgers", &convection, "Include Burgers equation coupling" );
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
  lsParams.set("Output Frequency", 1);    
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
  //converged->addStatusTest(update);
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

  // Create each part of the Brusselator problem class.  
  Equation_A ProblemA(Comm, NumGlobalNodes, "Temperature" );
  Equation_B ProblemB(Comm, NumGlobalNodes, "Species"     );
  Burgers  burgers   (Comm, NumGlobalNodes, "Burgers"     );

  // An interesting note: the order of solving each problem is based on the
  // order of adding.  For this decoupled problem, problem B is linear
  // with respect to its variables, whereas problem A is nonlinear wrt to its
  // variables.  The order of solution appears to strongly affect the rate
  // of convergence of the decoupled Brusselator.  Solving problem A first
  // dramatically reduces the number of total iterations.
  problemManager.addProblem(ProblemA);
  problemManager.addProblem(ProblemB);
  if( convection )
    problemManager.addProblem(burgers);

//  problemManager.createDependency("Temperature", "Species");
  problemManager.createDependency(ProblemA, ProblemB);
  problemManager.createDependency(ProblemB, ProblemA);

  if( convection )
  {
    problemManager.createDependency(ProblemA, burgers);
    problemManager.createDependency(ProblemB, burgers);
    problemManager.createDependency(burgers, ProblemA);
  }

  // Initialize time integration parameters
  int maxTimeSteps = 3;
  int timeStep = 0;
  double time = 0.;
  double dt = 0.100;
  problemManager.setAlldt(dt);
  
  problemManager.registerComplete(); // Trigger setup of groups, solvers, etc.

  problemManager.outputStatus(std::cout);

  // Print initial solution
  if( verbose )
    problemManager.outputSolutions( outputDir );
  
  // Identify the test problem
  if( outputUtils.isPrintType(NOX::Utils::TestDetails) )
    outputUtils.out() << "Starting epetra/MultiPhysics/example_brusselator.exe" << std::endl;

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
  while(timeStep < maxTimeSteps) 
  {
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
      problemManager.outputSolutions( outputDir, timeStep);

    map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator iter = problemManager.getProblems().begin(),
                                                                iter_end = problemManager.getProblems().end()   ;
    for( ; iter_end != iter; ++iter )
    {
      GenericEpetraProblem & problem = *(*iter).second;
      std::string msg = "Numerical-to-Gold Solution comparison for problem \"" + problem.getName() + "\"";

      // Get the gold copy to comapre against current solution
      std::string baseFileame = problemManager.createIOname( problem, timeStep );
      std::string goldFileame = goldDir + "gold_brusselator/" + baseFileame;

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

