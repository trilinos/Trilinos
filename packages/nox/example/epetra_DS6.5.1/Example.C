//  Simple 2 equation test for quadratic and cubic line searches 
//  from Dennis & Schnabel's book, chp 6.  The test problem is from
//  Example 6.5.1
/*  
 *    U0**2 + U1**2 - 2 = 0
 *    exp(U0-1) + U1**3 -2 = 0
 */

// Trilinos Objects
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"

// NOX Objects
#include "NOX_Epetra_Group.H"
#include "NOX_Solver_Newton.H"
#include "NOX_Status_All_Tests.H"

// User's application specific files 
#include "Problem_Interface.H" // Interface file to NOX
#include "DennisSchnabel.H"              

using namespace std;

int main(int argc, char *argv[])
{
  int ierr = 0, i, j;
  bool debug = false;

  // Initialize MPI
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
#endif

  // Create a communicator for Epetra objects
#ifdef EPETRA_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  int NumGlobalElements = 2;  // Hardcoded for D&S Example problem

  // A maximum of 2 procesors is allowed since there are only 2 equations
  if (NumProc >= 3) {
    cout << "ERROR: Maximum number of processors is 2!" << endl;
    exit(1);
  }

  // Create the Problem class.  This creates all required
  // Epetra objects for the problem and allows calls to the 
  // function (RHS) and Jacobian evaluation routines.
  DennisSchnabel Problem(NumGlobalElements, Comm);

  // Get the linear objects from the Problem
  Epetra_Vector& soln = Problem.getSolution();
  Epetra_RowMatrix& A = Problem.getJacobian();

  // Initialize Solution
  if (MyPID==0) {
    soln[0]=2.0;
    if (NumProc==1) soln[1]=0.5;
  } else soln[0]=0.5;

  cout << soln << endl;
  
  // Begin Nonlinear Solver ************************************

  // Create parameter list
  NOX::Parameter::List nlParams;
  nlParams.setParameter("Output Level", 4);
  nlParams.setParameter("MyPID", MyPID); 

  // Sublist for linear solver
  NOX::Parameter::List& lsParams = nlParams.sublist("Linear Solver");
  lsParams.setParameter("Max Iterations", 800);  
  lsParams.setParameter("Tolerance", 1e-4); 

  // Sublist for line search
  NOX::Parameter::List& searchParams = nlParams.sublist("Line Search");
  // searchParams.setParameter("Method", "Full Step");
  // searchParams.setParameter("Method", "Interval Halving");
  searchParams.setParameter("Method", "Polynomial");
  // searchParams.setParameter("Method", "More'-Thuente");
  searchParams.setParameter("Default Step", 1.0);

  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base class:
  // NLS_PetraGroupInterface
  Problem_Interface interface(Problem);

  // Create the shared Jacobian
  NOX::Epetra::SharedJacobian shareda(A);

  // Create the Groups 
  NOX::Epetra::Group grp(soln, shareda, interface); 
  grp.computeRHS();

  // Create the convergence tests
  NOX::Status::AbsResid absresid(1.0e-6);

  // Create the method
  NOX::Solver::Newton newton(grp, absresid, nlParams);
  NOX::Status::StatusType status = newton.solve();

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(newton.getSolutionGroup());
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

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}
