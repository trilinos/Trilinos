#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <string>
#ifdef EPETRA_MPI
#include "mpi.h"
#endif
#ifndef __cplusplus
#define __cplusplus
#endif

// Trilinos Objects
#include <Epetra_Comm.h>
#include <Epetra_SerialComm.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Time.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

// NOX Objects
#include <NOX_Epetra_Group.H>
#include <NOX_Solver_Newton.H>
#include <NOX_Status_AbsResid.H>

// User specific files 
#include <Problem_Interface.h>
#include <fill.h>

using namespace std;

int main(int argc, char *argv[])
{
  int ierr = 0, i, j;
  bool debug = false;

  // Initialize MPI
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  int size = 1; // Serial case (not using MPI)
  int rank = 0;
#endif

#ifdef EPETRA_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  // Get the number of local equations from the command line
  if (argc!=2) { 
    cout << "Usage: " << argv[0] << " number_of_elements" << endl;
    exit(1);
  }
  int NumGlobalElements = atoi(argv[1]) + 1;
  int IndexBase = 0;

  if (NumGlobalElements < NumProc) {
    cout << "numGlobalBlocks = " << NumGlobalElements 
	 << " cannot be < number of processors = " << NumProc << endl;
    exit(1);
  }

  // Construct a Source Map that puts approximately the same 
  // Number of equations on each processor in uniform global ordering

  Epetra_Map StandardMap(NumGlobalElements, 0, Comm);
  int NumMyElements = StandardMap.NumMyElements();
  int StandardMyGlobalElements[NumMyElements];
  StandardMap.MyGlobalElements(StandardMyGlobalElements);

  int OverlapNumMyElements;
  int OverlapMinMyGID;

  OverlapNumMyElements = NumMyElements + 2;
  if ((MyPID == 0) || (MyPID == NumProc - 1)) 
    OverlapNumMyElements --;

  if (MyPID==0) 
    OverlapMinMyGID = StandardMap.MinMyGID();
  else 
    OverlapMinMyGID = StandardMap.MinMyGID() - 1;

  int OverlapMyGlobalElements[OverlapNumMyElements];

  for (i = 0; i < OverlapNumMyElements; i ++) 
    OverlapMyGlobalElements[i] = OverlapMinMyGID + i;

  Epetra_Map* tmpMap;
  if (size == 1) 
    tmpMap = new Epetra_Map(StandardMap);  
  else 
    tmpMap = new Epetra_Map(-1, OverlapNumMyElements, 
			    OverlapMyGlobalElements, 0, Comm);
  
  Epetra_Map& OverlapMap(*tmpMap);

  // Construct Linear Objects  
  Epetra_Import Importer(StandardMap, OverlapMap);
  Epetra_Vector soln(StandardMap);
  Epetra_CrsMatrix AA(Copy, StandardMap, 3);

  // Construct the Matrix
  Fill LO;  
  LO.registerFillObjects(StandardMap, OverlapMap, Importer, Comm);

  // Fix things so we can do multiple imports on objects (make graph static)
  // Get Matrix structure
  LO.fillMatrix(&soln, NULL, &AA);
  // Create a graph
  Epetra_CrsMatrix A(Copy, AA.Graph());
  A.TransformToLocal();

  // Re-register so that static matrix A is used
  LO.registerFillObjects(StandardMap, OverlapMap, Importer, Comm);

  // Initialize Solution
  i = soln.PutScalar(1.0);
  
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
  searchParams.setParameter("Method", "Interval Halving");
  searchParams.setParameter("Default Step", 10.0);

  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of abstract base class:
  // NLS_PetraGroupInterface
  Problem_Interface interface;
  interface.registerFill(&LO);

  // Create the shared Jacobian
  NOX::Epetra::SharedJacobian shareda(A);

  // Create the Groups 
  NOX::Epetra::Group grp(soln, shareda, interface);  

  // Create the convergence tests
  NOX::Status::AbsResid test1(1.0e-6);

  // Create the method
  NOX::Solver::Newton newton(grp, test1, nlParams);
  NOX::Status::StatusType status = newton.solve();

  // End Nonlinear Solver **************************************

  // Print solution
  char file_name[25];
  FILE *ifp;
  (void) sprintf(file_name, "output.%d",MyPID);
  ifp = fopen(file_name, "w");
  for (i=0; i<NumMyElements; i++)
    fprintf(ifp, "%d  %E\n", StandardMap.MinMyGID()+i, soln[i]);
  fclose(ifp);

  //cout << "Final Solution" << (&soln)[0] <<endl;

  delete tmpMap;


  //delete &StandardGraph;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}
