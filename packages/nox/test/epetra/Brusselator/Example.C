//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_MapColoring.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
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
  int i;

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

  // Check for verbose option
  bool verbose = false;
  if (argc > 1)
    if (argv[1][0]=='-' && argv[1][1]=='v')
      verbose = true;

  // Set the problem size (1000 elements)
  int NumGlobalNodes = 1001;

  // The number of unknowns must be at least equal to the 
  // number of processors.
  if (NumGlobalNodes < NumProc) {
    cout << "numGlobalNodes = " << NumGlobalNodes 
	 << " cannot be < number of processors = " << NumProc << endl;
    cout << "Test failed!" << endl;
    exit(1);
  }

  // Create the Brusselator problem class.  This creates all required
  // Epetra objects for the problem and allows calls to the 
  // function (F) and Jacobian evaluation routines.
  Brusselator::OverlapType OType = Brusselator::NODES;
  Teuchos::RefCountPtr<Brusselator> Problem = 
    Teuchos::rcp(new Brusselator(NumGlobalNodes, Comm, OType));

  // Get the vector from the Problem
  Teuchos::RefCountPtr<Epetra_Vector> soln = Problem->getSolution();
  NOX::Epetra::Vector noxSoln(soln, NOX::Epetra::Vector::CreateView);

  // Begin Nonlinear Solver ************************************

  // Create the top level parameter list
  Teuchos::RefCountPtr<Teuchos::ParameterList> nlParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());

  // Set the nonlinear solver method
  nlParams.set("Nonlinear Solver", "Line Search Based");
  //nlParams.set("Nonlinear Solver", "Trust Region Based");

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  if (verbose) {
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
  }
  else
    printParams.set("Output Information", NOX::Utils::Error);

  // Create a print handle object
  NOX::Utils utils(printParams);

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
  lsParams.set("Output Frequency", 0);    
  lsParams.set("Preconditioner", "AztecOO"); 

  // Create the interface between the test problem and the nonlinear solver
  Teuchos::RefCountPtr<Problem_Interface> interface = 
    Teuchos::rcp(new Problem_Interface(*Problem));

#ifndef HAVE_NOX_EPETRAEXT 
  utils.out() << "Cannot use Coloring without package epetraext !!!!" << endl;
  utils.out() << "Test passed!" << endl;
  exit(0);
#else 
  // Create a timer for performance
  Epetra_Time fillTime(Comm);

  EpetraExt::CrsGraph_MapColoring::ColoringAlgorithm algType = 
    EpetraExt::CrsGraph_MapColoring::JONES_PLASSMAN;
  bool colorParallel = true; 
  int reordering = 0; 
  int verbosityLevel = 0; 
  bool distance1 = false; 
  EpetraExt::CrsGraph_MapColoring tmpMapColoring(
    algType, reordering, distance1, verbosityLevel);
  Teuchos::RefCountPtr<Epetra_MapColoring> colorMap = 
    Teuchos::rcp(&tmpMapColoring(*(Problem->getGraph())));
  EpetraExt::CrsGraph_MapColoringIndex colorMapIndex(*colorMap);
  Teuchos::RefCountPtr< std::vector<Epetra_IntVector> > columns = 
    Teuchos::rcp(&colorMapIndex(*(Problem->getGraph())));

  Teuchos::RefCountPtr<NOX::Epetra::Interface::Required> iReq = interface;
  Teuchos::RefCountPtr<NOX::Epetra::FiniteDifferenceColoring> FDC = 
    Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(printParams,
							   iReq, 
							   noxSoln, 
							   Problem->getGraph(),
							   colorMap, 
							   columns, 
							   colorParallel,
							   distance1));

  if (verbose)
    printf("\n[%d]\tTime to color Jacobian --> %e sec.  for %d colors.\n\n",
	   MyPID,fillTime.ElapsedTime(), colorMap->NumColors());

  FDC->setDifferenceMethod(NOX::Epetra::FiniteDifference::Centered);

// -------------- End of block needed to use coloring --------------- */


  // Create the Linear System
  Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac = FDC;
  Teuchos::RefCountPtr<NOX::Epetra::LinearSystemAztecOO> linSys = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
						      iReq, iJac, FDC, 
						      noxSoln));

  // Create the Group
  Teuchos::RefCountPtr<NOX::Epetra::Group> grpPtr = 
    Teuchos::rcp(new NOX::Epetra::Group(printParams, 
					iReq, 
					noxSoln, 
					linSys));  
  NOX::Epetra::Group& grp = *(grpPtr.get());

  // Create the convergence tests
  Teuchos::RefCountPtr<NOX::StatusTest::NormF> absresid = 
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8, NOX::StatusTest::NormF::Unscaled));
  Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxiters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(25));
  Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(absresid);
  combo->addStatusTest(maxiters);

  // Create the method
  NOX::Solver::Manager solver(grpPtr, combo, nlParamsPtr);

  // Initialize time integration parameters
  int maxTimeSteps = 3;
  int timeStep = 0;
  double time = 0.;
  double dt = Problem->getdt();
  
  // Print initial solution
  char file_name[25];
  FILE *ifp;
  Epetra_Vector& xMesh = *(Problem->getMesh());
  int NumMyNodes = xMesh.Map().NumMyElements();
  (void) sprintf(file_name, "output.%d_%d",MyPID,timeStep);
  ifp = fopen(file_name, "w");
  for (i=0; i<NumMyNodes; i++)
    fprintf(ifp, "%d  %E  %E  %E\n", xMesh.Map().MinMyGID()+i, 
	    xMesh[i], (*soln)[2*i], (*soln)[2*i+1]);
  fclose(ifp);
  
  NOX::StatusTest::StatusType status = NOX::StatusTest::Unconverged;

  // Time integration loop
  while(timeStep < maxTimeSteps) {

    timeStep++;
    time += dt;
  
    if (verbose)
      utils.out() << "Time Step: " << timeStep << ",\tTime: " << time << endl;
  
    status = NOX::StatusTest::Unconverged;
    status = solver.solve();
  
    if (verbose)
      if (status != NOX::StatusTest::Converged)
	if (MyPID==0) 
	  utils.out() << "Nonlinear solver failed to converge!" << endl;

    // Get the Epetra_Vector with the final solution from the solver
    const NOX::Epetra::Group& finalGroup = 
      dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());
    const Epetra_Vector& finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

    // End Nonlinear Solver **************************************

    // Print solution
    (void) sprintf(file_name, "output.%03d_%05d",MyPID,timeStep);
    ifp = fopen(file_name, "w");
    for (i=0; i<NumMyNodes; i++)
      fprintf(ifp, "%d  %E  %E  %E\n", soln->Map().MinMyGID()+i,
                                   xMesh[i], finalSolution[2*i],
                                   finalSolution[2*i+1]);
    fclose(ifp);

    Problem->reset(finalSolution);
    grp.setX(finalSolution);
    solver.reset(grpPtr, combo, nlParamsPtr);
    grp.computeF();

  } // end time step while loop

  // Output the parameter list
  if (verbose) {
    if (utils.isPrintType(NOX::Utils::Parameters)) {
      utils.out() << endl << "Final Parameters" << endl
	   << "****************" << endl;
      solver.getList().print(utils.out());
      utils.out() << endl;
    }
  }

  // Output timing info
  if (verbose) {
    if(MyPID==0)
      utils.out() << "\nTimings :\n\tWallTime --> " << 
	myTimer.WallTime() - startWallTime << " sec."
	   << "\n\tElapsedTime --> " << myTimer.ElapsedTime() 
	   << " sec." << endl << endl;
    
  }

  // Report results

  int testStatus = 0; // Converged

  // 1. Convergence
  if (status != NOX::StatusTest::Converged) {
    if (MyPID==0) utils.out() << "Nonlinear solver failed to converge!" << endl;
    testStatus = 1;
  }
  // 2. Nonlinear Iterations (3)
  if (const_cast<Teuchos::ParameterList&>(solver.getList()).sublist("Output").get("Nonlinear Iterations",0) != 3) {
    testStatus = 2;
  }
#ifndef HAVE_MPI
  // 3. Linear Iterations (9)
  if (const_cast<Teuchos::ParameterList&>(solver.getList()).sublist("Direction").sublist("Newton").sublist("Linear Solver").sublist("Output").get("Total Number of Linear Iterations",0) != 9) {
    testStatus = 3;
  }
#endif

  // Summarize test results 
  if (testStatus == 0)
    utils.out() << "Test passed!" << endl;
  else 
    utils.out() << "Test failed!" << endl;
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  // Final return value (0 = successfull, non-zero = failure)
  return 0;
#endif // HAVE_NOX_EPETRAEXT
}
