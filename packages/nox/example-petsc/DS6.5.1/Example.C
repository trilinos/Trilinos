//  Simple 2 equation test for quadratic and cubic line searches 
//  from Dennis & Schnabel's book, chp 6.  The test problem is from
//  Example 6.5.1
/*  
 *    U0**2 + U1**2 - 2 = 0
 *    exp(U0-1) + U1**3 -2 = 0
 */

static char help[] = "Solves Dennis & Schnabel example problem in parallel.\n\n";


// Petsc Objects
#include "petscsles.h"
/*
  Include "petscsles.h" so that we can use SLES solvers.  Note that this file
  automatically includes:
     petsc.h       - base PETSc routines   petscvec.h - vectors
     petscsys.h    - system routines       petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
*/

// NOX Library
#include "NOX.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Petsc_Vector.H"
#include "NOX_Petsc_SharedJacobian.H"
#include "NOX_Petsc_Group.H"

// User's application specific files 
#include "Problem_Interface.H" // Interface file to NOX
#include "DennisSchnabel.H"              

using namespace std;

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char *argv[])
{
  int ierr = 0, i;

  // Initialize MPI
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);CHKERRQ(ierr);

  // Get the process ID and the total number of processors
  int MyPID, NumProc;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&MyPID);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&NumProc);CHKERRQ(ierr);

  int NumGlobalElements = 2;  // Hardcoded for D&S Example problem

  // A maximum of 2 procesors is allowed since there are only 2 equations
  if (NumProc >= 3) {
    cout << "ERROR: Maximum number of processors is 2!" << endl;
    exit(1);
  }

  // Create the Problem class.  This creates all required
  // Petsc objects for the problem and allows calls to the 
  // function (RHS) and Jacobian evaluation routines.
  DennisSchnabel Problem(NumGlobalElements);

  // Get the vector from the Problem
  Vec& soln = Problem.getSolution();

  // Initialize Solution. For simplicity, this is done on all (both) procs.
  int globalIndex[2];
  double doubleArray[2];
  globalIndex[0] = 0;
  globalIndex[1] = 1;
  doubleArray[0] = 2.0;
  doubleArray[1] = 0.5;
  ierr = VecSetValues(soln,2,globalIndex,doubleArray,INSERT_VALUES);CHKERRQ(ierr);
//         Assemble vector, using the 2-step process:
//         VecAssemblyBegin(), VecAssemblyEnd()

  ierr = VecAssemblyBegin(soln);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(soln);CHKERRQ(ierr);

  // Begin Nonlinear Solver ************************************

  // Create parameter list
  NOX::Parameter::List nlParams;
  nlParams.setParameter("Output Level", 4);
  nlParams.setParameter("MyPID", MyPID);
  nlParams.setParameter("Nonlinear Solver", "Newton");
  //nlParams.setParameter("Nonlinear Solver", "NonlinearCG");  
  //nlParams.setParameter("Diagonal Precondition", "On");  // default = "Off"
  //nlParams.setParameter("Direction", "Steepest Descent");  // default
  //nlParams.setParameter("Direction", "Richardson");
  //nlParams.setParameter("Max Iterations", 100);
  //nlParams.setParameter("Orthogonalize", "Fletcher-Reeves");  // default
  //nlParams.setParameter("Orthogonalize", "Polak-Ribiere");
  //nlParams.setParameter("Restart Frequency", 5);  // default = 10
  //nlParams.setParameter("Output Frequency", 10);  // default = 1

  // Sublist for line search
  NOX::Parameter::List& searchParams = nlParams.sublist("Line Search");
  //searchParams.setParameter("Method", "Full Step");
  //searchParams.setParameter("Full Step", 0.01);
  //searchParams.setParameter("Method", "Interval Halving");
  searchParams.setParameter("Method", "Polynomial");
  //searchParams.setParameter("Method", "More'-Thuente");
  searchParams.setParameter("Default Step", 1.0000);
  //searchParams.setParameter("Recovery Step", 0.0001);
  //searchParams.setParameter("Minimum Step", 0.0001);

  // Sublist for direction
  NOX::Parameter::List& dirParams = nlParams.sublist("Direction");
  dirParams.setParameter("Method", "Newton");
  //dirParams.setParameter("Method", "Steepest Descent");
  //dirParams.setParameter("Method", "Dogleg Trust Region");
  //dirParams.setParameter("Method", "Broyden");

  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base class:
  // NLS_PetraGroupInterface
  Problem_Interface interface(Problem);

  // Sublist for linear solver
  NOX::Parameter::List& lsParams = dirParams.sublist("Linear Solver");
  lsParams.setParameter("Max Iterations", 800);  
  lsParams.setParameter("Tolerance", 1e-4);
  lsParams.setParameter("Iteration Output Frequency", 50);    
  lsParams.setParameter("Preconditioning Matrix Type", "None"); 

  // Create the Epetra_RowMatrix.  Uncomment one of the following:
  // 1. User supplied
  Mat& A = Problem.getJacobian();
  string jacType = "User Supplied";
  // 2. Matrix-Free
  //NOX::Epetra::MatrixFree A(interface, soln);
  // 3. Finite Difference
  //NOX::Epetra::FiniteDifference A(interface, soln);

  // Create the Group
  NOX::Petsc::Group grp(lsParams, interface, soln, A, jacType);
  grp.computeF();

  // Create the convergence tests
  // Create the convergence tests
  NOX::StatusTest::NormF testNormF(1.0e-6);
  NOX::StatusTest::MaxIters testMaxIters(1000);
  NOX::StatusTest::Combo combo(NOX::StatusTest::Combo::OR, testNormF, testMaxIters);

  // Create the method
  NOX::Solver::Manager solver(grp, combo, nlParams);

//  cout << "Made it here .." << endl;
//  cin.get();

  NOX::StatusTest::StatusType status = solver.solve();

  if (status != NOX::StatusTest::Converged)
    if (MyPID==0) 
      cout << "Nonlinear solver failed to converge!" << endl;

/*
  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());
  const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

  // End Nonlinear Solver **************************************

  // Print solution
  char file_name[25];
  FILE *ifp;
  int NumMyElements = soln.Map().NumMyElements();
  (void) sprintf(file_name, "output.%d",MyPID);
  ifp = fopen(file_name, "w");
  for (i=0; i<NumMyElements; i++)
    fprintf(ifp, "%d  %E\n", soln.Map().MinMyGID()+i, finalSolution[i]);
  fclose(ifp);

*/

  ierr = PetscFinalize();CHKERRQ(ierr);
 
return 0;
} // end main
