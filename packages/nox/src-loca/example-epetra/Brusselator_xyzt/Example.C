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
#include "LOCA.H"
#include "LOCA_Epetra.H"
#include "NOX_Common.H"

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_MapColoring.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

// Added to allow timings
#include "Epetra_Time.h"

// Headers needed for FD coloring 
#include <vector> 
#ifdef HAVE_NOX_EPETRAEXT       // Use epetraext package in Trilinos
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h" 
#endif

// User's application specific files 
#include "Problem_Interface.H" // Interface file to NOX
#include "Brusselator.H"              

#ifdef HAVE_MPI
#ifdef HAVE_NOX_EPETRAEXT
// Comment out following line for usual implicit time stepping on all procs
#define DO_XYZT 1
#endif
#endif

#ifdef DO_XYZT
#include "LOCA_EpetraNew_Interface_xyzt.H"              
#endif

using namespace std;

int main(int argc, char *argv[])
{
  int ierr = 0, i;

  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Get the number of elements from the command line
  if (argc<2) { 
    cout << "Usage: " << argv[0] << " number_of_elements" << endl;
    exit(1);
  }
  int NumGlobalNodes = atoi(argv[1]) + 1;

#ifdef HAVE_MPI

#ifdef DO_XYZT
  // MPI MANIPULATION FOR XYZT PROBLEMS
 
  int ierrmpi, size, rank;
  ierrmpi = MPI_Comm_size(MPI_COMM_WORLD, &size);
  ierrmpi = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int spatialProcs= 1; // default
  if (argc>2) { spatialProcs = atoi(argv[2]);}
  int timeStepsPerProc= 1; // default
  if (argc>3) { timeStepsPerProc = atoi(argv[3]);}

  if (size % spatialProcs != 0) {cout<<"ERROR: num spatial procs "<<spatialProcs
     << " does not divide into num total procs " << size << endl;  exit(-1);  }

  // Create split communicators, the size of spatial decomposition
  MPI_Comm split_MPI_Comm;
  int replica = rank/spatialProcs;
  ierrmpi =  MPI_Comm_split(MPI_COMM_WORLD, replica, rank, &split_MPI_Comm);

  // Construct 2 different epetra communicators
  Epetra_MpiComm Comm(split_MPI_Comm);
  Epetra_MpiComm globalComm(MPI_COMM_WORLD);

#else
  // Create a communicator for Epetra objects
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#endif
#else
  cout << "RUNNING IN SERIAL, NOT MPI " << endl;
  Epetra_SerialComm Comm;
#endif

  // Create and reset the Timer
  Epetra_Time myTimer(Comm);
  double startWallTime = myTimer.WallTime();

  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  // The number of unknowns must be at least equal to the 
  // number of processors.
  if (NumGlobalNodes < NumProc) {
    cout << "numGlobalNodes = " << NumGlobalNodes 
	 << " cannot be < number of processors = " << NumProc << endl;
    exit(1);
  }

  // Create the Brusselator problem class.  This creates all required
  // Epetra objects for the problem and allows calls to the 
  // function (F) and Jacobian evaluation routines.
  Brusselator::OverlapType OType = Brusselator::ELEMENTS;
  Brusselator Problem(NumGlobalNodes, Comm, OType);

  // Get the vector from the Problem
  Epetra_Vector& soln = Problem.getSolution();

  // Begin Nonlinear Solver ************************************

  // Create the top level parameter list

  NOX::Parameter::List paramList;

  // Create LOCA sublist
  NOX::Parameter::List& locaParamsList = paramList.sublist("LOCA");

  // Create the stepper sublist and set the stepper parameters
  NOX::Parameter::List& locaStepperList = locaParamsList.sublist("Stepper");
  locaStepperList.setParameter("Continuation Method", "Natural");
  //locaStepperList.setParameter("Continuation Method", "Arc Length");
  //locaStepperList.setParameter("Continuation Method", "Householder Arc Length");
  locaStepperList.setParameter("Continuation Parameter", "alpha");
  locaStepperList.setParameter("Initial Value", 0.6);
  locaStepperList.setParameter("Max Value", 100.0);
  locaStepperList.setParameter("Min Value", 0.05);
#ifdef DO_XYZT
  locaStepperList.setParameter("Max Steps", 4);
#else
  locaStepperList.setParameter("Max Steps", 0);// must be 0 so just a nonlinear solver
#endif
  locaStepperList.setParameter("Max Nonlinear Iterations", 15);
  locaStepperList.setParameter("Enable Arc Length Scaling", true);
  locaStepperList.setParameter("Goal Arc Length Parameter Contribution", 0.5);
  locaStepperList.setParameter("Max Arc Length Parameter Contribution", 0.7);
  locaStepperList.setParameter("Initial Scale Factor", 1.0);
  locaStepperList.setParameter("Min Scale Factor", 1.0e-8);
  locaStepperList.setParameter("Enable Tangent Factor Step Size Scaling",false);
  locaStepperList.setParameter("Min Tangent Factor", 0.8);
  locaStepperList.setParameter("Tangent Factor Exponent",1.5);

  // Create step size sublist
  NOX::Parameter::List& stepSizeList = locaParamsList.sublist("Step Size");
  stepSizeList.setParameter("Method", "Constant");
 // stepSizeList.setParameter("Method", "Adaptive");
  stepSizeList.setParameter("Initial Step Size", -0.1);
  stepSizeList.setParameter("Min Step Size", 1.0e-3);
  stepSizeList.setParameter("Max Step Size", 2000.0);
  stepSizeList.setParameter("Aggressiveness", 0.1);
  stepSizeList.setParameter("Failed Step Reduction Factor", 0.5);
  stepSizeList.setParameter("Successful Step Increase Factor", 1.00); // for constant

  // Create predictor sublist
  NOX::Parameter::List& predictorList = locaParamsList.sublist("Predictor");
  //predictorList.setParameter("Method", "Constant");
  //predictorList.setParameter("Method", "Tangent");
  predictorList.setParameter("Method", "Secant");

  // Create bifurcation sublist
    NOX::Parameter::List& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.setParameter("Method", "None");

  // Create Anasazi Eigensolver sublist (needs --with-loca-anasazi)
  locaStepperList.setParameter("Compute Eigenvalues",false);
  NOX::Parameter::List& aList = locaStepperList.sublist("Anasazi");
  aList.setParameter("Block Size", 1);
  aList.setParameter("Arnoldi Size", 10);
  aList.setParameter("NEV", 3);
  aList.setParameter("Tol", 2.0e-7);
  aList.setParameter("Convergence Check", 1);
  aList.setParameter("Restarts",2);
  aList.setParameter("Frequency",1);
  aList.setParameter("Debug Level",0);
  
  // Set the LOCA Utilities
  NOX::Parameter::List& locaUtilsList = locaParamsList.sublist("Utilities");
  locaUtilsList.setParameter("MyPID", MyPID);
  locaUtilsList.setParameter("Output Information", 
			     LOCA::Utils::Warning +
			     LOCA::Utils::StepperIteration +
			     LOCA::Utils::StepperDetails +
			     LOCA::Utils::Solver +
			     LOCA::Utils::SolverDetails +
			     LOCA::Utils::Parameters);

  // Create the "Solver" parameters sublist to be used with NOX Solvers
  NOX::Parameter::List& nlParams = paramList.sublist("NOX");
  // Set the nonlinear solver method
  nlParams.setParameter("Nonlinear Solver", "Line Search Based");
  //nlParams.setParameter("Nonlinear Solver", "Trust Region Based");

  // Set the printing parameters in the "Printing" sublist
  NOX::Parameter::List& printParams = nlParams.sublist("Printing");
  printParams.setParameter("MyPID", MyPID); 
  printParams.setParameter("Output Precision", 3);
  printParams.setParameter("Output Processor", 0);
  printParams.setParameter("Output Information", 
			NOX::Utils::OuterIteration + 
			NOX::Utils::OuterIterationStatusTest + 
			NOX::Utils::InnerIteration +
			NOX::Utils::Parameters + 
			NOX::Utils::Details + 
			NOX::Utils::Warning);

  // Sublist for line search 
  NOX::Parameter::List& searchParams = nlParams.sublist("Line Search");
  searchParams.setParameter("Method", "Full Step");
  //searchParams.setParameter("Method", "Interval Halving");
  //searchParams.setParameter("Method", "Polynomial");
  //searchParams.setParameter("Method", "NonlinearCG");
  //searchParams.setParameter("Method", "Quadratic");
  //searchParams.setParameter("Method", "More'-Thuente");

  // Sublist for direction
  NOX::Parameter::List& dirParams = nlParams.sublist("Direction");
  dirParams.setParameter("Method", "Newton");
  NOX::Parameter::List& newtonParams = dirParams.sublist("Newton");
    newtonParams.setParameter("Forcing Term Method", "Constant");
    //newtonParams.setParameter("Forcing Term Method", "Type 1");
    //newtonParams.setParameter("Forcing Term Method", "Type 2");
    //newtonParams.setParameter("Forcing Term Minimum Tolerance", 1.0e-4);
    //newtonParams.setParameter("Forcing Term Maximum Tolerance", 0.1);
  //dirParams.setParameter("Method", "Steepest Descent");
  //NOX::Parameter::List& sdParams = dirParams.sublist("Steepest Descent");
    //sdParams.setParameter("Scaling Type", "None");
    //sdParams.setParameter("Scaling Type", "2-Norm");
    //sdParams.setParameter("Scaling Type", "Quadratic Model Min");
  //dirParams.setParameter("Method", "NonlinearCG");
  //NOX::Parameter::List& nlcgParams = dirParams.sublist("Nonlinear CG");
    //nlcgParams.setParameter("Restart Frequency", 2000);
    //nlcgParams.setParameter("Precondition", "On");
    //nlcgParams.setParameter("Orthogonalize", "Polak-Ribiere");
    //nlcgParams.setParameter("Orthogonalize", "Fletcher-Reeves");

  // Sublist for linear solver for the Newton method
  NOX::Parameter::List& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.setParameter("Aztec Solver", "GMRES");  
  lsParams.setParameter("Max Iterations", 800);  
  lsParams.setParameter("Tolerance", 1e-6);
  lsParams.setParameter("Output Frequency", 50);    
  lsParams.setParameter("Preconditioner", "AztecOO"); 
  //lsParams.setParameter("Aztec Preconditioner", "ilu"); 
  //lsParams.setParameter("Overlap", 2);  
  //lsParams.setParameter("Graph Fill", 2); 
  //lsParams.setParameter("Aztec Preconditioner", "ilut"); 
  //lsParams.setParameter("Overlap", 2);   
  //lsParams.setParameter("Fill Factor", 2);   
  //lsParams.setParameter("Drop Tolerance", 1.0e-12);   
  //lsParams.setParameter("Aztec Preconditioner", "Polynomial"); 
  //lsParams.setParameter("Polynomial Order", 6); 

  // Create the interface between the test problem and the nonlinear solver
  Problem_Interface interface(Problem);

  // Print initial guess
  interface.printSolution(soln, locaStepperList.getParameter("Initial Value", 0.0));

  // Create the Epetra_RowMatrix
  Epetra_RowMatrix& A = Problem.getJacobian();

#ifdef DO_XYZT
  LOCA::EpetraNew::Interface::xyzt ixyzt(interface, interface, interface,
                                 soln, A, A, globalComm, timeStepsPerProc);
  Epetra_RowMatrix& Axyzt = ixyzt.getJacobian();
  Epetra_Vector& solnxyzt = ixyzt.getSolution();

  LOCA::EpetraNew::Interface::Required& iReq = ixyzt;

  // Create the Linear System
  NOX::EpetraNew::Interface::Jacobian& iJac = ixyzt;
  NOX::EpetraNew::LinearSystemAztecOO linSys(printParams, lsParams,
                                          iReq, iJac, Axyzt, solnxyzt);

  NOX::Epetra::Vector initialGuess(solnxyzt, NOX::DeepCopy, true);
#else
  // Use an Epetra Scaling object if desired
  Epetra_Vector scaleVec(soln);
  NOX::Epetra::Scaling scaling;
  scaling.addRowSumScaling(NOX::Epetra::Scaling::Left, scaleVec);
  //grp.setLinearSolveScaling(scaling);

  LOCA::EpetraNew::Interface::Required& iReq = interface;

  // Create the Linear System
  NOX::EpetraNew::Interface::Jacobian& iJac = interface;
  NOX::EpetraNew::LinearSystemAztecOO linSys(printParams, lsParams,
                                          iReq, iJac, A, soln);
                                          //&scaling);

  // Create the Group
  NOX::Epetra::Vector initialGuess(soln, NOX::DeepCopy, true);
#endif

  // Create and initialize the parameter vector
  LOCA::ParameterVector pVector;
  pVector.addParameter("alpha",0.6);
  pVector.addParameter("beta",2.0);

  LOCA::EpetraNew::Group grp(printParams, iReq, initialGuess, linSys, pVector);

  grp.computeF();

  // Create the convergence tests
  NOX::StatusTest::NormF absresid(1.0e-8, NOX::StatusTest::NormF::Unscaled);
  //NOX::StatusTest::NormF relresid(grp, 1.0e-2);
  //NOX::StatusTest::NormUpdate update(1.0e-5);
  //NOX::StatusTest::NormWRMS wrms(1.0e-2, 1.0e-8);
  //NOX::StatusTest::Combo converged(NOX::StatusTest::Combo::AND);
  //converged.addStatusTest(absresid);
  //converged.addStatusTest(relresid);
  //converged.addStatusTest(wrms);
  //converged.addStatusTest(update);
  NOX::StatusTest::MaxIters maxiters(50);
  NOX::StatusTest::Combo combo(NOX::StatusTest::Combo::OR);
  combo.addStatusTest(absresid);
  combo.addStatusTest(maxiters);

  // Create the method
  //NOX::Solver::Manager solver(grp, combo, nlParams);
  LOCA::Stepper stepper(grp, combo, paramList);

  // Initialize time integration parameters
#ifdef DO_XYZT
  int maxTimeSteps = 1; // No longer need a time integration loop
#else
  int maxTimeSteps = 5;
#endif
  int timeStep = 0;
  double time = 0.;
  double dt = Problem.getdt();
  
  // Time integration loop
  while(timeStep < maxTimeSteps) {

    timeStep++;
    time += dt;
  
    cout << "Time Step: " << timeStep << ",\tTime: " << time << endl;
  
//    NOX::StatusTest::StatusType status = solver.solve();
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status != LOCA::Abstract::Iterator::Finished)
       if (MyPID==0) cout << "Stepper failed to converge!" << endl;


    // Get the Epetra_Vector with the final solution from the solver
    const LOCA::EpetraNew::Group& finalGroup = 
      dynamic_cast<const LOCA::EpetraNew::Group&>(stepper.getSolutionGroup());
    const Epetra_Vector& finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

    // End Nonlinear Solver **************************************

#ifndef DO_XYZT
    // Not needed for continuation runs
    interface.printSolution(soln, locaStepperList.getParameter("Initial Value", 0.0));

    Problem.reset(finalSolution);
    grp.setX(finalSolution);
    stepper.reset(grp, combo, paramList);
    grp.computeF();
#endif

  } // end time step while loop

  // Output the parameter list
  NOX::Utils utils(printParams);
  if (utils.isPrintProcessAndType(NOX::Utils::Parameters)) {
    cout << endl << "Final Parameters" << endl
	 << "****************" << endl;
    stepper.getParameterList().print(cout);
    cout << endl;
  }

  // Output timing info
  if(MyPID==0)
    cout << "\nTimings :\n\tWallTime --> " << 
	    myTimer.WallTime() - startWallTime << " sec."
         << "\n\tElapsedTime --> " << myTimer.ElapsedTime() 
         << " sec." << endl << endl;

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

return ierr ;
}
