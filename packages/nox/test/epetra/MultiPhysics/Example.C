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
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

// Headers needed for FD coloring 
#include <vector> 
#ifdef HAVE_NOX_EPETRAEXT       // Use epetraext package in Trilinos
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h" 
#endif

// User's application specific files 
#include "Problem_Manager.H" 
#include "Problem_Interface.H" 
#include "Equation_A.H"              
#include "Equation_B.H"              

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

  // Begin Nonlinear Solver ************************************

  // NOTE: For now these parameters apply to all problems handled by
  // Problem_Manager.  Each problem could be made to have its own
  // parameter list as wwell as its own convergence test(s).

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
  //lsParams.setParameter("Preconditioning", "None");   
//  lsParams.setParameter("Preconditioning", "AztecOO: Jacobian Matrix");   
  lsParams.setParameter("Preconditioner", "AztecOO");
  //lsParams.setParameter("Graph Fill", 2);
  //lsParams.setParameter("Preconditioning", "AztecOO: User RowMatrix"); 
  //lsParams.setParameter("Preconditioning", "User Supplied Preconditioner");
  //lsParams.setParameter("Aztec Preconditioner", "ilu"); 
  //lsParams.setParameter("Overlap", 2);  
  //lsParams.setParameter("Graph Fill", 2); 
  //lsParams.setParameter("Aztec Preconditioner", "ilut"); 
  //lsParams.setParameter("Overlap", 2);   
  //lsParams.setParameter("Fill Factor", 2);   
  //lsParams.setParameter("Drop Tolerance", 1.0e-12);   
  //lsParams.setParameter("Aztec Preconditioner", "Polynomial"); 
  //lsParams.setParameter("Polynomial Order", 6); 

  // Create each part of the Brusselator problem class.  
  Equation_A ProblemA(Comm, NumGlobalNodes, "Temperature");
//  Equation_A ProblemA2(Comm, 11);
//  Equation_A ProblemA3(Comm, 501);
//  Equation_B ProblemB(Comm, 4, "Species");
  Equation_B ProblemB(Comm, NumGlobalNodes, "Species");
//  Equation_B ProblemB2(Comm, 11);
//  Equation_B ProblemB3(Comm, 501);

  // Create the Problem Manager
  Problem_Manager problemManager(Comm, true);

  // An interesting note: the order of solving each problem is based on the
  // order of adding.  For this decoupled problem, problem B is linear
  // with respect to its variables, whereas problem A is nonlinear wrt to its
  // variables.  The order of solution appears to strongly affect the rate
  // of convergence of the decoupled Brusselator.  Solving problem A first
  // dramatically reduces the number of total iterations.
  problemManager.addProblem(ProblemA);
//  problemManager.addProblem(ProblemA2);
//  problemManager.addProblem(ProblemA3);
  problemManager.addProblem(ProblemB);
//  problemManager.addProblem(ProblemB2);
//  problemManager.addProblem(ProblemB3);

//  problemManager.createDependency("Temperature", "Species");
  problemManager.createDependency(ProblemA, ProblemB);
//  problemManager.createDependency(ProblemA2, ProblemB3);
  problemManager.createDependency(ProblemB, ProblemA);
//  problemManager.createDependency(ProblemB2, ProblemA);
//  problemManager.createDependency(ProblemB3, ProblemA2);

  // Create the convergence tests
  // Note: as for the parameter list, both (all) problems use the same 
  // convergence test(s) for now, but each could have its own.
  NOX::StatusTest::NormF absresid(1.0e-8);
  NOX::StatusTest::NormUpdate update(1.0e-5);
  NOX::StatusTest::Combo converged(NOX::StatusTest::Combo::AND);
  converged.addStatusTest(absresid);
  converged.addStatusTest(update);
  NOX::StatusTest::MaxIters maxiters(25);
  NOX::StatusTest::Combo combo(NOX::StatusTest::Combo::OR);
  combo.addStatusTest(converged);
  combo.addStatusTest(maxiters);

  // Note that each problem could contain its own nlParams list as well as
  // its own convergence test(s). 
  problemManager.registerParameters(nlParams);
  problemManager.registerStatusTest(combo);
  problemManager.registerComplete(); // Trigger setup of groups, solvers, etc.

  problemManager.outputStatus();

  // Initialize time integration parameters
  int maxTimeSteps = 1;
  int timeStep = 0;
  double time = 0.;
  double dt = ProblemA.getdt();
  if( dt != ProblemB.getdt() )
    cout << "WARNING: Time steps differ between problems !!" << endl;
  
  // Print initial solution
  char file_name[25];
  FILE *ifp;
  Epetra_Vector& xMesh = ProblemA.getMesh();
  int NumMyNodes = xMesh.Map().NumMyElements();
  (void) sprintf(file_name, "output.%d_%d",MyPID,timeStep);
  ifp = fopen(file_name, "w");
  for (int i=0; i<NumMyNodes; i++)
    fprintf(ifp, "%d  %E  %E  %E\n", xMesh.Map().MinMyGID()+i, 
                                 xMesh[i], ProblemA.getSolution()[i], 
                                 ProblemB.getSolution()[i]);
  fclose(ifp);
  
  // Time integration loop
  while(timeStep < maxTimeSteps) {

    timeStep++;
    time += dt;
  
    cout << "Time Step: " << timeStep << ",\tTime: " << time << endl;
  
    // Solve decoupled
//    problemManager.solve(); // Need a status test check here ....
    // .... OR solve using matrix-free
    problemManager.solveMF(); // Need a status test check here ....
  
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

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}
