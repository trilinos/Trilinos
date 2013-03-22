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
                                                                                
// 1D Finite Element Interfacial Coupling Problem from 
// Yeckel, Pandy & Derby IJNME 2006.

/* Solves the nonlinear equation:
 */

// NOX Objects
#include "NOX.H"
#include "NOX_Epetra.H"
#include "NOX_TestCompare.H"

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
#include "NOX_Epetra_SchurCoupler.H" 
#include "Problem_Interface.H" 
#include "ConvDiff_PDE.H"              

// Added to allow timings
#include "Epetra_Time.h"

// Added temporarily
#include "Epetra_IntVector.h"

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
  CouplingSolveMethod method          = MATRIX_FREE   ;
  bool                verbose         = true          ;
  int                 NumGlobalNodes  = 20            ;
  int                 probSizeRatio   = 1             ;
  bool                useMatlab       = false         ;
  bool                doOffBlocks     = false         ;
  std::string         solvType        = "seidel"      ;
  // Coupling parameters
  double              alpha           = 0.50          ;
  double              beta            = 0.40          ;
  // Physical parameters
  double              radiation       = 5.67          ;
  double              initVal         = 0.995         ;
  std::string              outputDir       = "."           ;
  std::string              goldDir         = "."           ;


  clp.setOption<CouplingSolveMethod>( "solvemethod", &method, 4, SolveMethodValues, SolveMethodNames, "Selects the coupling method to use");
  clp.setOption( "verbose", "no-verbose", &verbose, "Verbosity on or off." );
  clp.setOption( "n", &NumGlobalNodes, "Number of elements" );
  clp.setOption( "nratio", &probSizeRatio, "Ratio of size of problem 2 to problem 1" );
  clp.setOption( "offblocks", "no-offblocks", &doOffBlocks, "Include off-diagonal blocks in preconditioning matrix" );
  clp.setOption( "matlab", "no-matlab", &useMatlab, "Use Matlab debugging engine" );
  clp.setOption( "solvType", &solvType, "Solve Type.  Valid choices are: jacobi, seidel" );
  clp.setOption( "alpha", &alpha, "Interfacial coupling coefficient, alpha" );
  clp.setOption( "beta", &beta, "Interfacial coupling coefficient, beta" );
  clp.setOption( "radiation", &radiation, "Radiation source term coefficient, R" );
  clp.setOption( "initialval", &initVal, "Initial guess for solution values" );
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
  //converged->addStatusTest(update);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(500));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> finiteValue = 
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);
  combo->addStatusTest(finiteValue);

  // Create the Problem Manager
  Problem_Manager problemManager(Comm, false, 0, useMatlab);

  // Note that each problem could contain its own nlParams list as well as
  // its own convergence test(s). 
  problemManager.registerParameters(nlParamsPtr);
  problemManager.registerStatusTest(combo);

  // Domain boundary temperaures
  double Tleft          = 0.98          ;
  double Tright         = 1.0           ;

  // Distinguish certain parameters needed for T1_analytic
  double peclet_1     	= 9.0           ;
  double peclet_2     	= 0.0           ;
  double kappa_1      	= 1.0           ;
  double kappa_2	= 0.1           ;

  double T1_analytic = ConvDiff_PDE::computeAnalyticInterfaceTemp( radiation, Tleft, Tright, kappa_2, peclet_1 );

  // Create Region 1 PDE
  std::string myName         = "Region_1"    ;
  double radiation_reg1 = 0.0           ;
  double xmin  		= 0.0           ;
  double xmax  		= 1.0           ;

  ConvDiff_PDE Reg1_PDE (
                  Comm, 
                  peclet_1,
                  radiation_reg1,
                  kappa_1,
                  alpha,
                  xmin,
                  xmax, 
                  Tleft,
                  T1_analytic,
                  NumGlobalNodes, 
                  myName  );

  // Override default initialization with values we want
  Reg1_PDE.initializeSolution(initVal);

  problemManager.addProblem(Reg1_PDE);


  // Create Region 2 PDE
  myName 		        = "Region_2"    ;
  xmin  		        = 1.0           ;
  xmax  		        = 2.0           ;

  ConvDiff_PDE Reg2_PDE (
                  Comm, 
                  peclet_2,
                  radiation,
                  kappa_2,
                  beta,
                  xmin,
                  xmax, 
                  T1_analytic,
                  Tright,
                  probSizeRatio*NumGlobalNodes, 
                  myName  );

  // For this problem involving interfacial coupling, the problems are given control
  // over whether or not and how to construct off-block contributions to the
  // Jacobian/Preconditioner matrix.  We explicitly told the problem manager to omit
  // off-blocks via the OffBlock_Manager class.  Here, we inform each problem.
  Reg1_PDE.setExpandJacobian( doOffBlocks );
  Reg2_PDE.setExpandJacobian( doOffBlocks );

  // Override default initialization with values we want
  Reg2_PDE.initializeSolution(initVal);

  problemManager.addProblem(Reg2_PDE);

  problemManager.createDependency(Reg1_PDE, Reg2_PDE);
  problemManager.createDependency(Reg2_PDE, Reg1_PDE);

  problemManager.registerComplete();

  // A consistencyy check associated with using BroydenOperator
  if( 0 )
  {
    Epetra_CrsGraph maskGraph(Copy, problemManager.getCompositeSoln()->Map(), 0);

    map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemIter = problemManager.getProblems().begin();
    map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator problemLast = problemManager.getProblems().end();

    // Loop over each problem being managed and ascertain its graph as well
    // as its graph from its dependencies
    for( ; problemIter != problemLast; ++problemIter ) 
    {
      GenericEpetraProblem & problem = *(*problemIter).second;
      int probId = problem.getId();

      // Get the indices map for copying data from this problem into 
      // the composite problem
      map<int, Teuchos::RCP<Epetra_IntVector> > & problemToCmpositeIndices = 
        problemManager.getProblemToCompositeIndices();
      Epetra_IntVector & problemIndices = *(problemToCmpositeIndices[probId]);

      // Get known dependencies on the other problem
      for( unsigned int k = 0; k < problem.getDependentProblems().size(); ++k) 
      {
        // Get the needed objects for the depend problem
        GenericEpetraProblem & dependProblem = *(problemManager.getProblems()[problem.getDependentProblems()[k]]);
        int dependId                         =  dependProblem.getId();
        Epetra_IntVector & dependIndices     = *(problemManager.getProblemToCompositeIndices()[dependId]);

        map<int, std::vector<int> > offBlockIndices;
        problem.getOffBlockIndices( offBlockIndices );

        map<int, std::vector<int> >::iterator indIter     = offBlockIndices.begin(),
                                         indIter_end = offBlockIndices.end()   ;

        for( ; indIter != indIter_end; ++indIter )
        {
          int compositeRow = problemIndices[(*indIter).first];
          std::vector<int> & colIndices = (*indIter).second;

          // Convert column indices to composite values
          for( unsigned int cols = 0; cols < colIndices.size(); ++cols )
            colIndices[cols] = dependIndices[ colIndices[cols] ];

          maskGraph.InsertGlobalIndices( compositeRow, colIndices.size(), &colIndices[0] );
        }
      }
    }
     maskGraph.FillComplete();

     std::cout << maskGraph << std::endl;

     NOX::Epetra::BroydenOperator * broydenOp = dynamic_cast<NOX::Epetra::BroydenOperator*>(
       problemManager.getJacobianOperator().get() );

    broydenOp->removeEntriesFromBroydenUpdate( maskGraph );
#ifdef HAVE_NOX_DEBUG
    broydenOp->outputActiveEntries();
#endif
  }

  problemManager.outputStatus(std::cout);

  std::cout << "\n\tAnalytic solution, T_1 = " << T1_analytic << "\n" << std::endl;

  // Print initial solution
  if( verbose )
    problemManager.outputSolutions( outputDir, 0 );

  // Identify the test problem
  if( outputUtils.isPrintType(NOX::Utils::TestDetails) )
    outputUtils.out() << "Starting epetra/MultiPhysics/example_yeckel.exe" << std::endl;

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

  // Identify the test problem
  if( outputUtils.isPrintType(NOX::Utils::TestDetails) )
    outputUtils.out() << "Starting epetra/MultiPhysics/example_yeckel.exe" << std::endl;

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

  // Solve the coupled problem
  switch( method ) 
  {
    case MATRIX_FREE :
      problemManager.solveMF(); // Need a status test check here ....
      break;

    case SCHUR_BASED :
      problemManager.solveSchurBased();
      break;

    case LOOSE_HARDCODED :
      problemManager.solve(); // Hard-coded loose coupling
      break;

    case LOOSE_LIBRARY :
    default            :
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
      
      Teuchos::RCP<NOX::StatusTest::MaxIters> fixedPt_maxiters = 
        Teuchos::rcp(new NOX::StatusTest::MaxIters(20));

      if( "jacobi" == solvType )
        nlParamsPtr->sublist("Solver Options").set("Fixed Point Iteration Type", "Jacobi");

      NOX::Multiphysics::Solver::Manager cplSolv( solversVec, dataExInterface, fixedPt_maxiters, nlParamsPtr );

      cplSolv.solve();

      // Refresh all problems with solutions from solver
      problemManager.copyAllGroupXtoProblems();

      // Reset all solver groups to force recomputation of residuals
      problemManager.resetAllCurrentGroupX();
    }
  }
  
  // Output timing info
  if( 0 == MyPID )
    std::cout << "\nTimings :\n\tWallTime --> " << myTimer.WallTime() - startWallTime << " sec."
         << "\n\tElapsedTime --> " << myTimer.ElapsedTime() << " sec." << std::endl << std::endl;

  if( verbose )
    problemManager.outputSolutions( outputDir, 1 );

  // Create a TestCompare class
  int status = 0;
  NOX::TestCompare tester( outputUtils.out(), outputUtils );
  double abstol = 1.e-4;
  double reltol = 1.e-4 ;

  map<int, Teuchos::RCP<GenericEpetraProblem> >::iterator iter = problemManager.getProblems().begin(),
                                                              iter_end = problemManager.getProblems().end()   ;
  for( ; iter_end != iter; ++iter )
  {
    ConvDiff_PDE & problem = dynamic_cast<ConvDiff_PDE &>( *(*iter).second );
    std::string msg = "Numerical-to-Exact Solution comparison for problem \"" + problem.getName() + "\"";

    // Need NOX::Epetra::Vectors for tests
    NOX::Epetra::Vector numerical ( problem.getSolution()     , NOX::Epetra::Vector::CreateView );
    NOX::Epetra::Vector analytic  ( problem.getExactSolution(), NOX::Epetra::Vector::CreateView );

    status += tester.testVector( numerical, analytic, reltol, abstol, msg );
  }

  // Summarize test results  
  if( status == 0 )
    outputUtils.out() << "Test passed!" << std::endl;
  else 
    outputUtils.out() << "Test failed!" << std::endl;

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return 0 ;
}

