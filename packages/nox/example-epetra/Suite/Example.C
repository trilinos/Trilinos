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

#include <fstream>
#include <vector>

using namespace std;

void setTraceFile(string);
void cleanupTraceFile();
#ifdef NOX_FILEOUTPUT 
streambuf* holdbuf = 0;
#endif
ofstream* outFile = 0;
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

  setTraceFile("1DfemNL_SD");

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
  // Get references to underlying solver parameter sublists
  NOX::Parameter::List& nlParams = Problem.getParameters();
  NOX::Parameter::List& printParams = nlParams.sublist("Printing");
  NOX::Parameter::List& lsParams = Problem.getlsParameters();
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
  Problem.outputResults(solver, printParams);

  cleanupTraceFile();

  cout << "\n\n\n\n\n" << endl;

// -----------------------------------------------------------
// Do another solve with different solver parameter options
// -----------------------------------------------------------

  NOX::Parameter::List& dirParams = nlParams.sublist("Direction");
  dirParams.setParameter("Method", "NonlinearCG");
  NOX::Parameter::List& nlcgParams = dirParams.sublist("Nonlinear CG");
    nlcgParams.setParameter("Restart Frequency", 1);

  // Reset the solution
  Problem.initializeSolution();
  grp.setX(Problem.getSolution()); 

  setTraceFile("1DfemNL_NLCG_noOrth");

  // Reset the solver
  solver.reset(grp, combo, nlParams);
  status = solver.solve();
  if (status != NOX::StatusTest::Converged)
    if (MyPID==0) 
      cout << "Nonlinear solver failed to converge!" << endl;
  outputResults(solver, printParams);

  cleanupTraceFile();

  cout << "\n\n\n\n\n" << endl;

// -----------------------------------------------------------
// Do another solve with different solver parameter options
// -----------------------------------------------------------

  dirParams.setParameter("Method", "NonlinearCG");
    nlcgParams.setParameter("Restart Frequency", 10*NumGlobalElements);
  NOX::Parameter::List& searchParams = nlParams.sublist("Line Search");
  searchParams.setParameter("Method", "NonlinearCG");

  // Reset the solution
  Problem.initializeSolution();
  grp.setX(Problem.getSolution()); 

  setTraceFile("1DfemNL_NLCG");

  // Reset the solver
  solver.reset(grp, combo, nlParams);
  status = solver.solve();
  if (status != NOX::StatusTest::Converged)
    if (MyPID==0) 
      cout << "Nonlinear solver failed to converge!" << endl;
  outputResults(solver, printParams);

  cleanupTraceFile();

  cout << "\n\n\n\n\n" << endl;

  delete ProblemPtr; ProblemPtr = 0;

// -----------------------------------------------------------
// Now do a different problem using yet another method
// -----------------------------------------------------------

  GenericProblem* ProblemPtr2 = new DennisSchnabel(DSnumGlobalElements, Comm);
  GenericProblem& Problem2 = *ProblemPtr2;

  // Initialize Solution
  Problem2.initializeSolution();

  setTraceFile("DSProblem");

  // Get the vector from the Problem
  Epetra_Vector& soln2 = Problem2.getSolution();
  Problem_Interface interface2(Problem2);
  Epetra_RowMatrix& A2 = Problem2.getJacobian();
  NOX::Epetra::Group grp2(Problem2.getParameters().sublist("Printing"),
                          Problem2.getlsParameters(),
                          interface2, soln2, A2); 
  grp2.computeF();

  // Reset the solver
//  For some reason, this doesn't work, so that a new solver is needed ??
//  solver.reset(grp2, combo, nlParams);
  NOX::Solver::Manager solver2(grp2, combo, Problem2.getParameters());
  status = solver2.solve();
  if (status != NOX::StatusTest::Converged)
    if (MyPID==0) 
      cout << "Nonlinear solver failed to converge!" << endl;
  outputResults(solver2, Problem2.getParameters().sublist("Printing"));

  cleanupTraceFile();

  delete ProblemPtr2; ProblemPtr2 = 0;

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

void setTraceFile(string fileBase)
{

#ifdef NOX_FILEOUTPUT
  // To allow Suite to write output to trace files, NOX must be built with
  //  the appropriate define set using --with-cppflags="-DNOX_FILEOUTPUT" 

  string traceFile = fileBase + ".trace";
  outFile = new ofstream(traceFile.c_str());
  // Redirect cout buffer to file
  holdbuf = cout.rdbuf(outFile->rdbuf());
  return;

#else
  return; // Default is to do nothing (ie leave cout alone)
#endif
  
}

void cleanupTraceFile()
{
#ifdef NOX_FILEOUTPUT
  // Restore original buffer to cout
  cout.rdbuf(holdbuf);
  outFile->close();
  delete outFile; outFile = 0;
  return;

#else
  return; // Default is to do nothing (ie leave cout alone)
#endif
}
