// 1D Finite Element Test Problem
/* Solves the nonlinear equation:
 *
 * d2u 
 * --- - k * u**2 = 0
 * dx2
 *
 * subject to @ x=0, u=1
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
#include "FiniteElementProblem.H"              

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

  // Get the number of elements from the command line
  if (argc!=2) { 
    cout << "Usage: " << argv[0] << " number_of_elements" << endl;
    exit(1);
  }
  int NumGlobalElements = atoi(argv[1]) + 1;

  // The number of unknowns must be at least equal to the 
  // number of processors.
  if (NumGlobalElements < NumProc) {
    cout << "numGlobalBlocks = " << NumGlobalElements 
	 << " cannot be < number of processors = " << NumProc << endl;
    exit(1);
  }

  // Create the FiniteElementProblem class.  This creates all required
  // Epetra objects for the problem and allows calls to the 
  // function (RHS) and Jacobian evaluation routines.
  FiniteElementProblem Problem(NumGlobalElements, Comm);

  // Get the linear objects from the Problem
  Epetra_Vector& soln = Problem.getSolution();
  Epetra_RowMatrix& A = Problem.getJacobian();

  // Initialize Solution
  soln.PutScalar(1.0);
  
  // Begin Nonlinear Solver ************************************

  // Create parameter list
  NOX::Parameter::List nlParams;
  nlParams.setParameter("Output Level", 4);
  nlParams.setParameter("MyPID", MyPID); 

  // Sublist for linear solver
  NOX::Parameter::List& lsParams = nlParams.sublist("Linear Solver");
  lsParams.setParameter("Max Iterations", 400);  
  lsParams.setParameter("Tolerance", 1e-4); 

  // Sublist for line search
  NOX::Parameter::List& searchParams = nlParams.sublist("Line Search");
  // searchParams.setParameter("Method", "Full Step");
  // searchParams.setParameter("Method", "Interval Halving");
  searchParams.setParameter("Method", "Polynomial");
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
  NOX::Status::Combo combo1(absresid, NOX::Status::Combo::AND);
  NOX::Status::RelResid relresid(grp.getNormRHS(), 1.0e-2);
  combo1.addTest(relresid);
  NOX::Status::MaxResid maxresid(1.0e-10);
  NOX::Status::Combo combo2(maxresid, NOX::Status::Combo::OR);
  combo2.addTest(combo1);

  // Create the method
  NOX::Solver::Newton newton(grp, combo2, nlParams);
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

  //cout << "Final Solution" << (&soln)[0] <<endl;
  //delete &StandardGraph;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}
