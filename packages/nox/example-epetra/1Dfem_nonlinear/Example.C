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
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

// NOX Objects
#include "NOX.H"
#include "NOX_Epetra.H"

// User's application specific files 
#include "Problem_Interface.H" // Interface file to NOX
#include "FiniteElementProblem.H"              

using namespace std;

int main(int argc, char *argv[])
{
  int ierr = 0, i;

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

  // Get the vector from the Problem
  Epetra_Vector& soln = Problem.getSolution();

  // Initialize Solution
  soln.PutScalar(1.0);
  
  // Begin Nonlinear Solver ************************************

  // Create parameter list
  NOX::Parameter::List nlParams;
  nlParams.setParameter("Output Level", 4);
  nlParams.setParameter("MyPID", MyPID); 
  //nlParams.setParameter("Nonlinear Solver", "Newton");
  nlParams.setParameter("Nonlinear Solver", "Line Search");
  //nlParams.setParameter("Nonlinear Solver", "Trust Region"); 
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
  searchParams.setParameter("Method", "Full Step");
  //searchParams.setParameter("Full Step", 0.01);
  //searchParams.setParameter("Method", "Interval Halving");
  //searchParams.setParameter("Method", "Polynomial");
  //searchParams.setParameter("Method", "More'-Thuente");
  searchParams.setParameter("Default Step", 1.0000);
  //searchParams.setParameter("Recovery Step", 0.0001);
  //searchParams.setParameter("Minimum Step", 0.0001);

  // Sublist for direction
  NOX::Parameter::List& dirParams = nlParams.sublist("Direction");
  dirParams.setParameter("Method", "Newton");
  //dirParams.setParameter("Method", "Steepest Descent");

  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base class:
  // NLS_PetraGroupInterface
  Problem_Interface interface(Problem);

  // Sublist for linear solver
  NOX::Parameter::List& lsParams = dirParams.sublist("Linear Solver");
  lsParams.setParameter("Aztec Solver", "CG");  
  lsParams.setParameter("Max Iterations", 800);  
  lsParams.setParameter("Tolerance", 1e-4);
  lsParams.setParameter("Output Frequency", 50);    
  lsParams.setParameter("Preconditioning", "None");   
  lsParams.setParameter("Scaling", "None");          
  //lsParams.setParameter("Preconditioning", "AztecOO: Jacobian Matrix");   
  //lsParams.setParameter("Preconditioning", "AztecOO: User RowMatrix"); 
  //lsParams.setParameter("Preconditioning", "User Supplied Preconditioner");
  //lsParams.setParameter("Aztec Preconditioner", "ilu"); 
  //lsParams.setParameter("Overlap", 2);  
  //lsParams.setParameter("Graph Fill", 2); 
  //lsParams.setParameter("Aztec Preconditioner", "ilut"); 
  //lsParams.setParameter("Overlap", 2);   
  //lsParams.setParameter("Fill Factor", 2.0);   
  //lsParams.setParameter("Drop Tolerance", 1.0e-6);   
  //lsParams.setParameter("Aztec Preconditioner", "Polynomial"); 
  //lsParams.setParameter("Polynomial Order", 6); 

  // Create the Epetra_RowMatrix.  Uncomment one or more of the following:
  // 1. User supplied (Epetra_RowMatrix)
  Epetra_RowMatrix& A = Problem.getJacobian();
  // 2. Matrix-Free (Epetra_Operator)
  //NOX::Epetra::MatrixFree AA(interface, soln);
  // 3. Finite Difference (Epetra_RowMatrix)
  //NOX::Epetra::FiniteDifference AAA(interface, soln);

  // Create the Group
  NOX::Epetra::Group grp(lsParams, interface, soln, A); 
  //NOX::Epetra::Group grp(lsParams, interface, soln, AA, AAA); 
  grp.computeRHS();

  // Create the convergence tests
  NOX::Status::AbsResid absresid(1.0e-6);
  NOX::Status::RelResid relresid(grp.getNormRHS(), 1.0e-2);
  NOX::Status::Combo converged(NOX::Status::Combo::AND);
  converged.addTest(absresid);
  converged.addTest(relresid);
  NOX::Status::MaxResid maxresid(1.0e-10);
  NOX::Status::MaxIters maxiters(2000);
  NOX::Status::Combo combo(NOX::Status::Combo::OR);
  combo.addTest(converged);
  combo.addTest(maxresid);
  combo.addTest(maxiters);

  // Create the method
  NOX::Solver::Manager solver(grp, combo, nlParams);
  NOX::Status::StatusType status = solver.solve();

  if (status != NOX::Status::Converged)
    if (MyPID==0) 
      cout << "Nonlinear solver failed to converge!" << endl;

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

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}
