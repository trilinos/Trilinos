//  Simple 2 equation test for quadratic and cubic line searches 
//  from Dennis & Schnabel's book, chp 6.  The test problem is from
//  Example 6.5.1
/*  
 *    U0**2 + U1**2 - 2 = 0
 *    exp(U0-1) + U1**3 -2 = 0
 */

static char help[] = 
       "Solves Dennis & Schnabel example problem in parallel.\n\n";


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
  ierr = VecSetValues(soln,2,globalIndex,doubleArray,
                                         INSERT_VALUES);CHKERRQ(ierr);
//         Assemble vector, using the 2-step process:
//         VecAssemblyBegin(), VecAssemblyEnd()

  ierr = VecAssemblyBegin(soln);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(soln);CHKERRQ(ierr);

  // Begin Nonlinear Solver ************************************

  // Create the top level parameter list
  NOX::Parameter::List nlParams;

  // Specify nonlinear solver method
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
  //searchParams.setParameter("Method", "Full Step");
  //searchParams.setParameter("Method", "Interval Halving");
  searchParams.setParameter("Method", "Polynomial");
  //searchParams.setParameter("Method", "Secant");
  //searchParams.setParameter("Method", "Quadratic");
  //searchParams.setParameter("Method", "More'-Thuente");

  // Sublist for direction
  NOX::Parameter::List& dirParams = nlParams.sublist("Direction");
  //
  // Popular choices include the following (others may also exist)
  //
  dirParams.setParameter("Method", "Newton");
  NOX::Parameter::List& newtonParams = dirParams.sublist("Newton");
    newtonParams.setParameter("Forcing Term Method", "Constant");
    //newtonParams.setParameter("Forcing Term Method", "Type 1");
    //newtonParams.setParameter("Forcing Term Method", "Type 2");
    //newtonParams.setParameter("Forcing Term Minimum Tolerance", 1.0e-4);
    //newtonParams.setParameter("Forcing Term Maximum Tolerance", 0.1);
  // OR
  //dirParams.setParameter("Method", "Steepest Descent");
  //NOX::Parameter::List& sdParams = dirParams.sublist("Steepest Descent");
    //sdParams.setParameter("Scaling Type", "None");
    //sdParams.setParameter("Scaling Type", "2-Norm");
    //sdParams.setParameter("Scaling Type", "Quadratic Model Min");
  // OR
  //dirParams.setParameter("Method", "NonlinearCG");
  //NOX::Parameter::List& nlcgParams = dirParams.sublist("Nonlinear CG");
    //nlcgParams.setParameter("Restart Frequency", 2000);
    //nlcgParams.setParameter("Precondition", "On");
    //nlcgParams.setParameter("Orthogonalize", "Polak-Ribiere");
    //nlcgParams.setParameter("Orthogonalize", "Fletcher-Reeves");

  // Sublist for linear solver
  // Note that preconditioning options as well as the following can be
  // specified on the command line or via the file .petscrc
  // See Petsc documentation for more info.
  NOX::Parameter::List& lsParams = dirParams.sublist("Linear Solver");
  lsParams.setParameter("Max Iterations", 800);  
  lsParams.setParameter("Tolerance", 1e-4);
  lsParams.setParameter("Iteration Output Frequency", 50);    
  lsParams.setParameter("Preconditioning Matrix Type", "None"); 

  // Create the interface between the test problem and the nonlinear solver
  // This is created using inheritance of the abstract base class:
  // NOX::Petsc::Interface
  Problem_Interface interface(Problem);

  // Get a reference to the Petsc_RowMatrix created in Problem_Interface.  
  Mat& A = Problem.getJacobian();

  // Create the Group
  NOX::Petsc::Group grp(lsParams, interface, soln, A);
  grp.computeF(); // Needed to establish the initial convergence state

  // Create the convergence tests
  NOX::StatusTest::NormF testNormF(1.0e-6);
  NOX::StatusTest::MaxIters testMaxIters(1000);
  NOX::StatusTest::Combo combo(NOX::StatusTest::Combo::OR, testNormF, testMaxIters);

  // Create the method and solve
  NOX::Solver::Manager solver(grp, combo, nlParams);
  NOX::StatusTest::StatusType status = solver.solve();

  if (status != NOX::StatusTest::Converged)
    if (MyPID==0) 
      cout << "Nonlinear solver failed to converge!" << endl;

  // Get the Petsc_Vector with the final solution from the solver
  const NOX::Petsc::Group& finalGroup = 
      dynamic_cast<const NOX::Petsc::Group&>(solver.getSolutionGroup());
  const Vec& finalSolution = (dynamic_cast<const NOX::Petsc::Vector&>
        (finalGroup.getX())).getPetscVector();

  // End Nonlinear Solver **************************************

  // Print solution
  if(MyPID==0)
    printf("Final solution :\n");
  char file_name[25];
  FILE *ifp;
  (void) sprintf(file_name, "output.%d",MyPID);
  double* finalVals;
  ifp = fopen(file_name, "w");
  VecGetArray( finalSolution, &finalVals );
  for(int i=0; i<3-NumProc; i++) {
    printf("(proc %d)\t%d\t%e\n",MyPID,i,finalVals[i]);
    fprintf(ifp, "%d  %E\n",i, finalVals[i]);
  }
  VecRestoreArray( finalSolution, &finalVals );
  fclose(ifp);

  ierr = PetscFinalize();CHKERRQ(ierr);
 
return 0;

} // end main
