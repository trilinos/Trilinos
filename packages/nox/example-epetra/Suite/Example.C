// NOX Library
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

// User's application specific files 
#include "Problem_Interface.H" // Interface file to NOX
#include "DennisSchnabel.H"              
#include "FiniteElementProblem.H"              

#include <vector>

using namespace std;

void outputResults(NOX::Solver::Manager& solver, NOX::Parameter::List& print);

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

  int DSnumGlobalElements = 2;  // Hardcoded for D&S Example problem

  // A maximum of 2 procesors is allowed since there are only 2 equations
  if (NumProc >= 3) {
    cout << "ERROR: Cannot run DennisSchnabel problem on more than 2"
         << "processors !!" << endl;
    exit(1);
  }

  // Create the Problem class.  This creates all required
  // Epetra objects for the problem and allows calls to the 
  // function (RHS) and Jacobian evaluation routines.

  // Eventaully, the hope is to pack these into a vector<GenericProblems>
  // and loop (iterate) over them.

  GenericProblem* ProblemPtr;
  ProblemPtr = new FiniteElementProblem(NumGlobalElements, Comm);

  GenericProblem& Problem = *ProblemPtr;

  // Initialize Solution
  Problem.initializeSolution();

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
  //searchParams.setParameter("Method", "Full Step");
  //searchParams.setParameter("Method", "Interval Halving");
  searchParams.setParameter("Method", "Polynomial");
  //searchParams.setParameter("Method", "Secant");
  //searchParams.setParameter("Method", "Quadratic");
  //searchParams.setParameter("Method", "More'-Thuente");

  // Sublist for direction
  NOX::Parameter::List& dirParams = nlParams.sublist("Direction");
  //dirParams.setParameter("Method", "Newton");
  //NOX::Parameter::List& newtonParams = dirParams.sublist("Newton");
    //newtonParams.setParameter("Forcing Term Method", "Constant");
    //newtonParams.setParameter("Forcing Term Method", "Type 1");
    //newtonParams.setParameter("Forcing Term Method", "Type 2");
    //newtonParams.setParameter("Forcing Term Minimum Tolerance", 1.0e-4);
    //newtonParams.setParameter("Forcing Term Maximum Tolerance", 0.1);
    //NOX::Parameter::List& lsParams = newtonParams.sublist("Linear Solver");
  dirParams.setParameter("Method", "Steepest Descent");
  NOX::Parameter::List& sdParams = dirParams.sublist("Steepest Descent");
    NOX::Parameter::List& lsParams = sdParams.sublist("Linear Solver");
    //sdParams.setParameter("Scaling Type", "None");
    //sdParams.setParameter("Scaling Type", "2-Norm");
    //sdParams.setParameter("Scaling Type", "Quadratic Model Min");
  //dirParams.setParameter("Method", "NonlinearCG");
  //NOX::Parameter::List& nlcgParams = dirParams.sublist("Nonlinear CG");
    //nlcgParams.setParameter("Restart Frequency", 2000);
    //nlcgParams.setParameter("Precondition", "On");
    //nlcgParams.setParameter("Orthogonalize", "Polak-Ribiere");
    //nlcgParams.setParameter("Orthogonalize", "Fletcher-Reeves");

  // Sublist for linear solver
  //lsParams.setParameter("Aztec Solver", "GMRES");  
  //lsParams.setParameter("Max Iterations", 800);  
  //lsParams.setParameter("Tolerance", 1e-4);
  //lsParams.setParameter("Output Frequency", 50);    
  //lsParams.setParameter("Scaling", "None");             
  //lsParams.setParameter("Scaling", "Row Sum");          
  //lsParams.setParameter("Preconditioning", "None");   
  //lsParams.setParameter("Preconditioning", "AztecOO: Jacobian Matrix");   
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

  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base class:
  // NLS_PetraGroupInterface
  Problem_Interface interface(Problem);
  
  // Create the Epetra_RowMatrix.  Uncomment one or more of the following:
  // 1. User supplied (Epetra_RowMatrix)
  Epetra_RowMatrix& A = Problem.getJacobian();
  // 2. Matrix-Free (Epetra_Operator)
  //NOX::Epetra::MatrixFree AA(interface, soln);
  // 3. Finite Difference (Epetra_RowMatrix)
  //NOX::Epetra::FiniteDifference AAA(interface, soln);

  // Create the Group
  NOX::Epetra::Group grp(printParams, lsParams, interface, soln, A); 
  grp.computeF();

  // Create the convergence tests
  NOX::StatusTest::NormF testNormF(1.0e-6);
  NOX::StatusTest::MaxIters testMaxIters(1000);
  NOX::StatusTest::Combo combo(NOX::StatusTest::Combo::OR, testNormF, testMaxIters);

  // Create the method
  NOX::Solver::Manager solver(grp, combo, nlParams);
  NOX::StatusTest::StatusType status = solver.solve();

  if (status != NOX::StatusTest::Converged)
    if (MyPID==0) 
      cout << "Nonlinear solver failed to converge!" << endl;
  outputResults(solver, printParams);

  cout << "\n\n\n\n\n" << endl;

// -----------------------------------------------------------
// Do another solve with different solver parameter options
// -----------------------------------------------------------

  dirParams.setParameter("Method", "Residual-Based: Picard");
  NOX::Parameter::List& picardParams = dirParams.sublist("Picard");
    //lsParams = picardParams.sublist("Linear Solver");

  // Reset the solver
  solver.reset(grp, combo, nlParams);
  status = solver.solve();
  if (status != NOX::StatusTest::Converged)
    if (MyPID==0) 
      cout << "Nonlinear solver failed to converge!" << endl;
  outputResults(solver, printParams);

  cout << "\n\n\n\n\n" << endl;

// -----------------------------------------------------------
// Do another solve with different solver parameter options
// -----------------------------------------------------------

  dirParams.setParameter("Method", "Residual-Based: NonlinearCG");
  NOX::Parameter::List& nlcgParams = dirParams.sublist("Nonlinear CG");
    //lsParams = picardParams.sublist("Linear Solver");

  // Reset the solver
  solver.reset(grp, combo, nlParams);
  status = solver.solve();
  if (status != NOX::StatusTest::Converged)
    if (MyPID==0) 
      cout << "Nonlinear solver failed to converge!" << endl;
  outputResults(solver, printParams);

  cout << "\n\n\n\n\n" << endl;

// -----------------------------------------------------------
// Now do a different problem using yet another method
// -----------------------------------------------------------

//  delete ProblemPtr;
  GenericProblem* ProblemPtr2 = new DennisSchnabel(DSnumGlobalElements, Comm);
  GenericProblem& Problem2 = *ProblemPtr2;

  // Initialize Solution
  Problem2.initializeSolution();

  // Get the vector from the Problem
  Epetra_Vector& soln2 = Problem2.getSolution();
  dirParams.setParameter("Method", "Newton");
  NOX::Parameter::List& newtonParams = dirParams.sublist("Newton");
    newtonParams.setParameter("Forcing Term Method", "Constant");
    NOX::Parameter::List& lsParams2 = newtonParams.sublist("Linear Solver");
      lsParams2.setParameter("Aztec Solver", "GMRES");  
      lsParams2.setParameter("Max Iterations", 800);  
      lsParams2.setParameter("Tolerance", 1e-4);
      lsParams2.setParameter("Output Frequency", 50);    
      lsParams2.setParameter("Preconditioning", "AztecOO: Jacobian Matrix");   
  Problem_Interface interface2(Problem2);
  Epetra_RowMatrix& A2 = Problem2.getJacobian();
  NOX::Epetra::Group grp2(printParams, lsParams2, interface2, soln2, A2); 
  grp2.computeF();

  // Reset the solver
//  For some reason, this doesn't work, so that a new solver is needed ??
//  solver.reset(grp2, combo, nlParams);
  NOX::Solver::Manager solver2(grp2, combo, nlParams);
  status = solver2.solve();
  if (status != NOX::StatusTest::Converged)
    if (MyPID==0) 
      cout << "Nonlinear solver failed to converge!" << endl;
  outputResults(solver2, printParams);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
 
return ierr;
} // end main

void outputResults(NOX::Solver::Manager& solver, 
                   NOX::Parameter::List& printParams)
{
  // Output the parameter list
  NOX::Utils utils(printParams);
  if (utils.isPrintProcessAndType(NOX::Utils::Parameters)) {
    cout << endl << "Final Parameters" << endl
	 << "****************" << endl;
    solver.getParameterList().print(cout);
    cout << endl;
  }

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup = 
      dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());
  const Epetra_Vector& finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>
        (finalGroup.getX())).getEpetraVector();

  // Print solution
  char file_name[25];
  FILE *ifp;
  int NumMyElements = finalSolution.Map().NumMyElements();
  (void) sprintf(file_name, "output.%d",finalSolution.Map().Comm().MyPID());
  ifp = fopen(file_name, "w");
  for (int i=0; i<NumMyElements; i++)
    fprintf(ifp,"%d  %E\n",finalSolution.Map().MinMyGID()+i,finalSolution[i]);
  fclose(ifp);
}


