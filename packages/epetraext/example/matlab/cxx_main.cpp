#include "Epetra_ConfigDefs.h"

#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_DataAccess.h"
#include "Epetra_CrsMatrix.h"

#include "EpetraExt_MatlabEngine.h" // contains the EpetraExt_MatlabEngine class

int main(int argc, char *argv[]) {

// standard Epetra MPI/Serial Comm startup	
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  
  int MyPID = comm.MyPID();
  int ierr = 0;
  
  // setup MatlabEngine
  cout << "going to startup a matlab process...\n";
  EpetraExt::EpetraExt_MatlabEngine engine (comm);
  cout << "matlab started\n";
  
  // setup an array of doubles to be used for the examples
  int M = 20;
  int numGlobalElements = M * comm.NumProc();
  int N = 3;
  int numMyEntries = M * N;
  double* A = new double[numMyEntries];
  double* Aptr = A;
  int startValue = numMyEntries * MyPID;

  for(int col=0; col < N; col++) {
	for(int row=0; row < M; row++) {
          *Aptr++ = startValue++;
      }
  }

  // setup an array of ints to be used for the examples
  int* intA = new int[numMyEntries];
  int* intAptr = intA;
  int intStartValue = numMyEntries * MyPID;
  for(int i=0; i < M*N; i++) {
      *intAptr++ = intStartValue++;
  }
  
  // construct a map to be used by distributed objects
  Epetra_Map map (numGlobalElements, 0, comm);
  
  // CrsMatrix example
  // constructs a globally distributed CrsMatrix and then puts it into Matlab
  cout << MyPID << " constructing CrsMatrix...\n";
  Epetra_CrsMatrix crsMatrix (Copy, map, N);
  int* indices = new int[N];
  for (int col=0; col < N; col++) {
    indices[col] = col;	  
  }
  
  double value = startValue;
  int minMyGID = map.minGID();
  for (int row=0; row < M; row++) {
    for (int col=0; col < N; col++) {
      values[col] = value++;
    }
      
    crsMatrix.InsertGlobalValues(minMyGID + row, N, values, indices);
  }
  
  crsMatrix.FillComplete();
  cout << MyPID << " CrsMatrix constructed\n";
  cout << MyPID << " putting CrsMatrix into Matlab as CRSM\n";
  ierr = engine.PutRowMatrix(crsMatrix, "CRSM", false);
  if (ierr) {
    cout << "There was an error in engine.PutRowMatrix(crsMatrix, "CRSM", false): " << ierr << endl;
  }
  
  // MultiVector example
  // constructs a globally distributed MultiVector and then puts it into Matlab
  cout << MyPID << " constructing MultiVector...\n";
  Epetra_MultiVector multiVector (Copy, map, A, M, N);
  cout << MyPID << " MultiVector constructed\n";
  cout << MyPID << " putting MultiVector into Matlab as MV\n";
  ierr = engine.PutMultiVector(multiVector, "MV");
  if (ierr) {
    cout << "There was an error in engine.PutMultiVector(multiVector, "MV"): " << ierr << endl;
  }
  
  // SerialDenseMatrix example
  // constructs a SerialDenseMatrix on every PE
  cout << MyPID << " constructing a SerialDenseMatrix...\n";
  Epetra_SerialDenseMatrix sdMatrix (Copy, A, M, M, N);
  cout << MyPID << " SerialDenseMatrix constructed\n";
  cout << MyPID << " putting SerialDenseMatrix from PE0 into Matlab as SDM_PE0\n";
  // since the third parameter is left out, the SerialDenseMatrix from PE0 is used by default
  ierr = engine.PutSerialDenseMatrix(sdMatrix, "SDM_PE0");
  if (ierr) {
    cout << "There was an error in engine.PutSerialDenseMatrix(sdMatrix, "SDM_PE0"): " << ierr << endl;
  }
  if (comm.NumProc() > 1) {
    cout << MyPID << " putting SerialDenseMatrix from PE1 into Matlab as SDM_PE1\n";
    // specifying 1 as the third parameter will put the SerialDenseMatrix from PE1 into Matlab
    ierr = engine.PutSerialDenseMatrix(sdMatrix, "SDM_PE1", 1);
    if (ierr) {
      cout << "There was an error in engine.PutSerialDenseMatrix(sdMatrix, "SDM_PE1", 1): " << ierr << endl;
    }
  }
  

  // SerialDenseVector example
  // constructs a SerialDenseVector on every PE
  cout << MyPID << " constructing a SerialDenseVector...\n";
  Epetra_SerialDenseVector sdVector (Copy, A, M);
  cout << MyPID << " SerialDenseVector constructed\n";
  // since the third parameter is left out, the SerialDenseMatrix from PE0 is used by default
  cout << MyPID << " putting SerialDenseVector from PE0 into Matlab as SDV_PE0\n";
  ierr = engine.PutSerialDenseMatrix(sdVector, "SDV_PE0");
  if (ierr) {
    cout << "There was an error in engine.PutSerialDenseMatrix(sdVector, "SDV_PE0"): " << ierr << endl;
  }
  if (comm.NumProc() > 1) {
    cout << MyPID << " putting SerialDenseVector from PE1 into Matlab as SDV_PE1\n";
    // specifying 1 as the third parameter will put the SerialDenseVector from PE1 into Matlab
    ierr = engine.PutSerialDenseMatrix(sdMatrix, "SDV_PE1", 1);
    if (ierr) {
      cout << "There was an error in engine.PutSerialDenseMatrix(sdMatrix, "SDV_PE1", 1): " << ierr << endl;
    }
  }

  // IntSerialDenseMatrix example
  // constructs a IntSerialDenseMatrix on every PE
  cout << MyPID << " constructing a IntSerialDenseMatrix...\n";
  Epetra_IntSerialDenseMatrix isdMatrix (Copy, intA, M, M, N);
  cout << MyPID << " IntSerialDenseMatrix constructed\n";
  // since the third parameter is left out, the IntSerialDenseMatrix from PE0 is used by default
  cout << MyPID << " putting IntSerialDenseMatrix from PE0 into Matlab as ISDM_PE0\n";
  ierr = engine.PutIntSerialDenseMatrix(isdMatrix, "ISMD_PE0");
  if (ierr) {
    cout << "There was an error in engine.PutIntSerialDenseMatrix(isdMatrix, "ISMD_PE0"): " << ierr << endl;
  }
  if (comm.NumProc() > 1) {
    cout << MyPID << " putting IntSerialDenseMatrix from PE1 into Matlab as ISDM_PE1\n";
    // specifying 1 as the third parameter will put the IntSerialDenseMatrix from PE1 into Matlab
    ierr = engine.PutSerialDenseMatrix(sdMatrix, "ISDM_PE1", 1);
    if (ierr) {
      cout << "There was an error in engine.PutSerialDenseMatrix(sdMatrix, "ISDM_PE1", 1): " << ierr << endl;
    }
  }


  // IntSerialDenseVector example
  // constructs a IntSerialDenseVector on every PE
  cout << MyPID << " constructing a IntSerialDenseVector...\n";
  Epetra_IntSerialDenseVector isdVector (Copy, intA, M);
  cout << MyPID << " IntSerialDenseVector constructed\n";
  // since the third parameter is left out, the IntSerialDenseVector from PE0 is used by default
  cout << MyPID << " putting IntSerialDenseVector from PE0 into Matlab as ISDV_PE0\n";
  ierr = engine.PutIntSerialDenseMatrix(isdVector, "ISDV_PE0");
  if (ierr) {
    cout << "There was an error in engine.PutIntSerialDenseMatrix(isdVector, "ISDV_PE0"): " << ierr << endl;
  }
  if (comm.NumProc() > 1) {
    cout << MyPID << " putting IntSerialDenseVector from PE1 into Matlab as ISDV_PE1\n";
    // specifying 1 as the third parameter will put the IntSerialDenseVector from PE1 into Matlab
    ierr = engine.PutSerialDenseMatrix(sdMatrix, "ISDV_PE1", 1);
    if (ierr) {
      cout << "There was an error in engine.PutSerialDenseMatrix(sdMatrix, "ISDV_PE1", 1): " << ierr << endl;
    }
  }
  
  /*while(1) {

	// do nothing
	}*/

  
  char s [200] ;
  char matlabBuffer [1024 * 16];
  /*if (comm.NumProc() == 1) {
  int err;
  while(1) {
      // Prompt the user and get a string
      printf(">> ");
      if (fgets(s, BUFSIZE, stdin) == NULL) {
          printf("Bye\n");
          break ;
      }
      printf ("command :%s:\n", s) ;
      
      // Send the command to MATLAB
      // output goes to stdout
      err = engine.EvalString(s, matlabBuffer, MATLABBUF);
      if (err != 0) {
          printf("there was an error: %d", err);
		  err = 0;
      }
      else {
      	  printf("Matlab Output:\n%s", matlabBuffer);
      }
  }
  }*/

  
  //engine.EvalString("size(TEST)", matlabBuffer, MATLABBUF);
  //cout << matlabBuffer << "\n";
  //engine.EvalString("TEST", matlabBuffer, MATLABBUF);
  //cout << matlabBuffer << "\n";

  engine.EvalString("whos", matlabBuffer, MATLABBUF);
  cout << matlabBuffer << "\n";
  
  cout << "\n" << comm.MyPID() << " all done\n";

// standard finalizer for Epetra MPI Comms
#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  // we need to delete engine because the MatlabEngine finalizer shuts down the Matlab process associated with this example
  // if we don't delete the Matlab engine, then this example application will not shut down properly
  delete &engine;
  return(0);
}
