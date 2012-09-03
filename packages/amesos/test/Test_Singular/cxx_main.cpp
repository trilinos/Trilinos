#include "Amesos_ConfigDefs.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Util.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Export.h"
#include "Amesos.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Galeri_ReadHB.h"

//============ //
// main driver //
//============ //

int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

// ================= //
  // reading HB matrix //
  // ================= //

  // HB files are for serial matrices. Hence, only
  // process 0 reads this matrix (and if present
  // solution and RHS). Then, this matrix will be redistributed
  // using epetra capabilities.
  // All variables that begin with "read" refer to the
  // HB matrix read by process 0.
  Epetra_Map* readMap;
  Epetra_CrsMatrix* readA;
  Epetra_Vector* readx;
  Epetra_Vector* readb;
  Epetra_Vector* readxexact;

  // Name of input file with a numerical singularity
  std::string matrix_file="bcsstm05_ns.rua";
  try
  {
    Galeri::ReadHB(matrix_file.c_str(), Comm, readMap,
                   readA, readx, readb, readxexact);
  }
  catch(...)
  {
    std::cout << "Caught exception, maybe file name is incorrect" << std::endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#else
    // not to break test harness
    exit(EXIT_SUCCESS);
#endif
  }

  // Create uniform distributed map.
  // Note that linear map are used for simplicity only!
  // Amesos (through Epetra) can support *any* map.

  Epetra_Map* mapPtr = 0;
  if(readMap->GlobalIndicesInt())
    mapPtr = new Epetra_Map((int) readMap->NumGlobalElements(), 0, Comm);
  else if(readMap->GlobalIndicesLongLong())
    mapPtr = new Epetra_Map(readMap->NumGlobalElements(), 0, Comm);
  else
    assert(false);

  Epetra_Map& map = *mapPtr;

  // Create the distributed matrix, based on Map.
  Epetra_CrsMatrix A(Copy, map, 0);

  const Epetra_Map &OriginalMap = readA->RowMatrixRowMap() ;
  assert (OriginalMap.SameAs(*readMap));
  Epetra_Export exporter(OriginalMap, map);

  Epetra_Vector x(map);          // distributed solution
  Epetra_Vector b(map);          // distributed rhs
  Epetra_Vector xexact(map);     // distributed exact solution

  // Exports from data defined on processor 0 to distributed.
  x.Export(*readx, exporter, Add);
  b.Export(*readb, exporter, Add);
  xexact.Export(*readxexact, exporter, Add);
  A.Export(*readA, exporter, Add);
  A.FillComplete();

  // deletes memory
  delete readMap;
  delete readA;
  delete readx;
  delete readb;
  delete readxexact;

  // Creates an epetra linear problem, contaning matrix
  // A, solution x and rhs b.
  Epetra_LinearProblem problem(&A,&x,&b);

  // =========== //
  // AMESOS PART //
  // =========== //

  // Create the factory
  Amesos factory;

  // Create the solver
  std::string solverType = "Klu";
  Teuchos::RCP<Amesos_BaseSolver> solver = Teuchos::rcp( factory.Create(solverType, problem) );

  // Perform symbolic factorization
  int symRet = solver->SymbolicFactorization();
  if (symRet != 0) {
    std::cout << "Processor "<< map.Comm().MyPID() << " : Symbolic factorization did not complete!" << std::endl;
  }

  // Perform numeric factorization
  int numRet = solver->NumericFactorization();
  if (numRet != 0) {
    std::cout << "Processor "<< map.Comm().MyPID() << " : Numeric factorization did not complete!" << std::endl;
  } 

  // Check that all processors returned error -22 (Numerically Singular)!
  int minRet = 0;
  Comm.MinAll( &numRet, &minRet, 1 );

  if (minRet != NumericallySingularMatrixError) {
    if (Comm.MyPID()==0)
      std::cout << std::endl << "End Result: TEST FAILED" << std::endl;
  }
  else {
    if (Comm.MyPID()==0)
      std::cout << std::endl << "End Result: TEST PASSED" << std::endl;
  }

  delete mapPtr;
 
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);
}
