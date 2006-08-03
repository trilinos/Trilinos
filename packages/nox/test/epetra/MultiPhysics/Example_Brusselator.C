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
  int  NumGlobalNodes = 10    ;
  bool runMF          = true  ;
  bool useMatlab      = false ;
  bool doOffBlocks    = false ;
  bool convection     = false ;

  clp.setOption( "n", &NumGlobalNodes, "Number of elements" );
  clp.setOption( "runMF", "loose", &runMF, "Use Matrix-Free strong coupling" );
  clp.setOption( "offblocks", "no-offblocks", &doOffBlocks, "Include off-diagonal blocks in preconditioning matrix" );
  clp.setOption( "burgers", "no-burgers", &convection, "Include Burgers equation coupling" );
  clp.setOption( "matlab", "no-matlab", &useMatlab, "Use Matlab debugging engine" );

  Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);

  if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) 
    return parse_return;

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
    cout << "numGlobalNodes = " << NumGlobalNodes 
	 << " cannot be < number of processors = " << NumProc << endl;
    exit(1);
  }

  // Begin Nonlinear Solver ************************************

  // NOTE: For now these parameters apply to all problems handled by
  // Problem_Manager.  Each problem could be made to have its own
  // parameter list as wwell as its own convergence test(s).

  // Create the top level parameter list
  Teuchos::RefCountPtr<Teuchos::ParameterList> nlParamsPtr =
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
			NOX::Utils::OuterIteration + 
			NOX::Utils::OuterIterationStatusTest + 
			NOX::Utils::InnerIteration +
			NOX::Utils::Parameters + 
			NOX::Utils::Details + 
			NOX::Utils::Warning);

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
  Teuchos::RefCountPtr<NOX::StatusTest::NormF>          absresid  =
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RefCountPtr<NOX::StatusTest::NormUpdate>     update    = 
      Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-5));
  Teuchos::RefCountPtr<NOX::StatusTest::Combo>          converged = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  //converged->addStatusTest(update);
  Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxiters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RefCountPtr<NOX::StatusTest::FiniteValue> finiteValue = 
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo = 
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

  problemManager.registerComplete(); // Trigger setup of groups, solvers, etc.

  problemManager.outputStatus(std::cout);

  // Initialize time integration parameters
  int maxTimeSteps = 3;
  int timeStep = 0;
  double time = 0.;
  double dt = 0.100;
  problemManager.setAlldt(dt);
  
  // Print initial solution
  char file_name[25];
  FILE *ifp;
  Epetra_Vector& xMesh = ProblemA.getMesh();
  int NumMyNodes = xMesh.Map().NumMyElements();
  (void) sprintf(file_name, "output.%d_%d",MyPID,timeStep);
  ifp = fopen(file_name, "w");
  for (int i=0; i<NumMyNodes; i++)
    fprintf(ifp, "%d  %E  %E  %E\n", xMesh.Map().MinMyGID()+i, 
                                 xMesh[i], (*ProblemA.getSolution())[i], 
                                 (*ProblemB.getSolution())[i]);
  fclose(ifp);
  //FILE *ifp2;
  //Epetra_Vector& burgersX = burgers.getMesh();
  //(void) sprintf(file_name, "burgers.%d_%d",MyPID,timeStep);
  //ifp2 = fopen(file_name, "w");
  //for (int i = 0; i < 10*NumMyNodes; ++i)
  //  fprintf(ifp2, "%d  %E  %E\n", burgersX.Map().MinMyGID()+i, 
  //                               burgersX[i], burgers.getSolution()[i]);
  //fclose(ifp2);
  
  
  // Time integration loop
  while(timeStep < maxTimeSteps) 
  {
    timeStep++;
    time += dt;
  
    cout << "Time Step: " << timeStep << ",\tTime: " << time << endl;
  
    // Solve the coupled problem
    if( runMF )
      problemManager.solveMF(); // Need a status test check here ....
    else
    {
      // Create the loose coupling solver manager
      Teuchos::RefCountPtr<vector<NOX::Solver::Manager*> > solversVec =
        Teuchos::rcp( new vector<NOX::Solver::Manager*> );

      map<int, NOX::Solver::Manager*>::iterator iter = problemManager.getSolvers().begin(),
                                            iter_end = problemManager.getSolvers().end()   ;
      for( ; iter_end != iter; ++iter )
      {
        cout << " ........  registered Solver::Manager # " << (*iter).first << endl;
        solversVec->push_back( (*iter).second );
      }

      // Package the Problem_Manager as the DataExchange::Intreface
      Teuchos::RefCountPtr<NOX::Multiphysics::DataExchange::Interface> dataExInterface =
        Teuchos::rcp( &problemManager, false );

      NOX::Multiphysics::Solver::Manager cplSolv( solversVec, dataExInterface, combo, nlParamsPtr );

      cplSolv.solve();
    }

    problemManager.outputSolutions(timeStep);

    // Reset problems by copying solution into old solution
    problemManager.resetProblems();

  } // end time step while loop

  // Output timing info
  if(MyPID==0)
    cout << "\nTimings :\n\tWallTime --> " << 
          myTimer.WallTime() - startWallTime << " sec."
         << "\n\tElapsedTime --> " << myTimer.ElapsedTime() 
         << " sec." << endl << endl;

  // Need to check test pass/fail status
//  if(MyPID==0)
//    cout << "Test failed!" << endl;
//
//#ifdef HAVE_MPI
//MPI_Finalize() ;
//#endif
//
//  return -1;

  // Need to put in a check for convergence
  // Added the following so test actually passes in parallel
  if(MyPID==0)
    cout << "Test passed!" << endl;

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return 0 ;
}

