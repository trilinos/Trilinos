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

using namespace std;

int main(int argc, char *argv[])
{
  int ierr = 0, i;

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

  // Create and reset the Timer
  Epetra_Time myTimer(Comm);
  double startWallTime = myTimer.WallTime();

  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  // Get the number of elements from the command line
  if (argc!=2) { 
    cout << "Usage: " << argv[0] << " number_of_elements" << endl;
    exit(1);
  }
  int NumGlobalNodes = atoi(argv[1]) + 1;

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
  Brusselator::OverlapType OType = Brusselator::NODES;
  Brusselator Problem(NumGlobalNodes, Comm, OType);

  // Get the vector from the Problem
  Epetra_Vector& soln = Problem.getSolution();

  // Begin Nonlinear Solver ************************************

  // Create the top level parameter list
  NOX::Parameter::List nlParams;

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
  lsParams.setParameter("Tolerance", 1e-4);
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

  // Create the Epetra_RowMatrix.  Uncomment one or more of the following:
  // 1. User supplied (Epetra_RowMatrix)
  //Epetra_RowMatrix& A = Problem.getJacobian(); // Default
  //NOX::EpetraNew::Interface::Jacobian& iJac = interface;
  // 2. Matrix-Free (Epetra_Operator)
  //NOX::EpetraNew::MatrixFree A(interface, soln);
  // 3. Finite Difference (Epetra_RowMatrix)
  //NOX::EpetraNew::FiniteDifference FD(interface, soln, Problem.getGraph());
  //  A.setDifferenceMethod(NOX::Epetra::FiniteDifference::Backward);
  // 4. Jacobi Preconditioner
  //NOX::Epetra::JacobiPreconditioner Prec(soln);
  // 5. Finite Difference with Coloring......uncomment the following
// -------------- Uncomment this block to use coloring --------------- //
#ifdef HAVE_NOX_EPETRAEXT 
  // Create a timer for performance
  Epetra_Time fillTime(Comm);

  bool verbose = false; 
  int reordering = 0; 
  bool useParallelColoring = true; 
  EpetraExt::CrsGraph_MapColoring::ColoringAlgorithm algType = 
    EpetraExt::CrsGraph_MapColoring::ALGO_GREEDY;
  EpetraExt::CrsGraph_MapColoring tmpMapColoring(algType, verbose, reordering, useParallelColoring); 
  Epetra_MapColoring* colorMap = &tmpMapColoring(Problem.getGraph());
  EpetraExt::CrsGraph_MapColoringIndex colorMapIndex(*colorMap);
  vector<Epetra_IntVector>* columns = &colorMapIndex(Problem.getGraph());
  NOX::EpetraNew::FiniteDifferenceColoring FDC(interface, soln, 
                                             Problem.getGraph(),
                                             *colorMap, *columns, 
                                             useParallelColoring);
                                             //1.0e-6, 0.e0);
  printf("\n[%d]\tTime to color Jacobian --> %e sec.  for %d colors.\n\n",
              MyPID,fillTime.ElapsedTime(), colorMap->NumColors());

  FDC.setDifferenceMethod(NOX::EpetraNew::FiniteDifference::Centered);

#else 
  cout << "Cannot use Coloring without package epetraext !!!!" << endl;
  exit(0);
#endif 
// -------------- End of block needed to use coloring --------------- */


  // Use an Epetra Scaling object if desired
  Epetra_Vector scaleVec(soln);
  NOX::Epetra::Scaling scaling;
  scaling.addRowSumScaling(NOX::Epetra::Scaling::Left, scaleVec);
  //grp.setLinearSolveScaling(scaling);

  NOX::EpetraNew::Interface::Required& iReq = interface;

  // Create the Linear System
  NOX::EpetraNew::Interface::Jacobian& iJac = FDC;
  NOX::EpetraNew::LinearSystemAztecOO linSys(printParams, lsParams,
                                          iReq, iJac, FDC, soln);
                                          //&scaling);

  // Create the Group
  NOX::Epetra::Vector initialGuess(soln, NOX::DeepCopy, true);
  NOX::EpetraNew::Group grp(printParams, iReq, initialGuess, linSys);
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
  NOX::StatusTest::MaxIters maxiters(25);
  NOX::StatusTest::Combo combo(NOX::StatusTest::Combo::OR);
  combo.addStatusTest(absresid);
  combo.addStatusTest(maxiters);

  // Create the method
  NOX::Solver::Manager solver(grp, combo, nlParams);

  // Initialize time integration parameters
  int maxTimeSteps = 3;
  int timeStep = 0;
  double time = 0.;
  double dt = Problem.getdt();
  
  // Print initial solution
  char file_name[25];
  FILE *ifp;
  Epetra_Vector& xMesh = Problem.getMesh();
  int NumMyNodes = xMesh.Map().NumMyElements();
  (void) sprintf(file_name, "output.%d_%d",MyPID,timeStep);
  ifp = fopen(file_name, "w");
  for (i=0; i<NumMyNodes; i++)
    fprintf(ifp, "%d  %E  %E  %E\n", xMesh.Map().MinMyGID()+i, 
                                 xMesh[i], soln[2*i], soln[2*i+1]);
  fclose(ifp);
  
  // Time integration loop
  while(timeStep < maxTimeSteps) {

    timeStep++;
    time += dt;
  
    cout << "Time Step: " << timeStep << ",\tTime: " << time << endl;
  
    NOX::StatusTest::StatusType status = solver.solve();
  
    if (status != NOX::StatusTest::Converged)
      if (MyPID==0) 
        cout << "Nonlinear solver failed to converge!" << endl;

    // Get the Epetra_Vector with the final solution from the solver
    const NOX::EpetraNew::Group& finalGroup = 
      dynamic_cast<const NOX::EpetraNew::Group&>(solver.getSolutionGroup());
    const Epetra_Vector& finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

    // End Nonlinear Solver **************************************

    // Print solution
    (void) sprintf(file_name, "output.%03d_%05d",MyPID,timeStep);
    ifp = fopen(file_name, "w");
    for (i=0; i<NumMyNodes; i++)
      fprintf(ifp, "%d  %E  %E  %E\n", soln.Map().MinMyGID()+i,
                                   xMesh[i], finalSolution[2*i],
                                   finalSolution[2*i+1]);
    fclose(ifp);

    Problem.reset(finalSolution);
    grp.setX(finalSolution);
    solver.reset(grp, combo, nlParams);
    grp.computeF();

  } // end time step while loop

  // Output the parameter list
  NOX::Utils utils(printParams);
  if (utils.isPrintProcessAndType(NOX::Utils::Parameters)) {
    cout << endl << "Final Parameters" << endl
	 << "****************" << endl;
    solver.getParameterList().print(cout);
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

/* end main
*/
return ierr ;
}
