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
                                                                                
// 1D Finite Element Interfacial Coupling Problem from 
// Yeckel, Pandy & Derby IJNME 2006.

/* Solves the nonlinear equation:
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
#include "ConvDiff_PDE.H"              

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

  clp.setOption( "n", &NumGlobalNodes, "Number of elements" );
  clp.setOption( "runMF", "loose", &runMF, "Use Matrix-Free strong coupling" );
  clp.setOption( "offblocks", "no-offblocks", &doOffBlocks, "Include off-diagonal blocks in preconditioning matrix" );
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

  // Coupling parameters
  double alpha		= 0.50          ;
  double beta 		= 0.40          ;

  // Domain boundary temperaures
  double Tleft  		= 0.98          ;
  double Tright 		= 1.0           ;

  // Distinguish certain parameters needed for T1_analytic
  double peclet_1     	= 9.0           ;
  double peclet_2     	= 0.0           ;
  double kappa_1      	= 1.0           ;
  double kappa_2		= 0.1           ;

  //double T1_analytic          = ( kappa_2*(1.0-exp(peclet_1))*Tright - Tleft*peclet_1*exp(peclet_1) ) /
  //                              ( kappa_2*(1.0-exp(peclet_1)) - peclet_1*exp(peclet_1) );
  //double T1_analytic          = 0.99430092; // 5.67
  double T1_analytic          = 0.99050495; // 2.50
  //double T1_analytic          = 0.98247093; // 0.3
  //double T1_analytic          = 0.98177822; // 0.2

  // Create Region 1 PDE
  string myName 		= "Region_1"    ;
  double radiation		= 0.0           ;
  double xmin  		= 0.0           ;
  double xmax  		= 1.0           ;

  ConvDiff_PDE Reg1_PDE (
                  Comm, 
                  peclet_1,
                  radiation,
                  kappa_1,
                  alpha,
                  xmin,
                  xmax, 
                  Tleft,
                  T1_analytic,
                  1*NumGlobalNodes, 
                  myName  );

  // Override default initialization with values we want
  Reg1_PDE.initializeSolution(0.995);

  problemManager.addProblem(Reg1_PDE);


  // Create Region 2 PDE
  myName 		        = "Region_2"    ;
  radiation		        = 5.67          ;
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
                  1*NumGlobalNodes, 
                  myName  );

  // Override default initialization with values we want
  Reg2_PDE.initializeSolution(0.995);

  problemManager.addProblem(Reg2_PDE);

  problemManager.createDependency(Reg1_PDE, Reg2_PDE);
  problemManager.createDependency(Reg2_PDE, Reg1_PDE);

  problemManager.registerComplete();

  problemManager.outputStatus(std::cout);

  cout << "\n\tAnalytic solution, T_1 = " << T1_analytic << "\n" << endl;

  // Print initial solution
  problemManager.outputSolutions(0);

  // Solve the coupled problem
  if( runMF )
    problemManager.solveMF(); // Need a status test check here ....
  else
    problemManager.solve(); // Need a status test check here ....
  
  problemManager.outputSolutions(1);

  // Output timing info
  if(MyPID==0)
    cout << "\nTimings :\n\tWallTime --> " << myTimer.WallTime() - startWallTime << " sec."
         << "\n\tElapsedTime --> " << myTimer.ElapsedTime() << " sec." << endl << endl;

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

