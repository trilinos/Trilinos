#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Flops.h"
#include "Epetra_Export.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "../epetra_test_err.h"
 
int main(int argc, char *argv[])
{
  int ierr = 0, i, j, forierr = 0;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );

#else

  Epetra_SerialComm Comm;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;



  // char tmp; if (Comm.MyPID()==0) { cout << "Press any key to continue..."<< endl;  cin >> tmp;} Comm.Barrier();

  Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
              << " is alive."<<endl;

  bool verbose1 = verbose;

  // Redefine verbose to only print on PE 0
  if (verbose && Comm.MyPID()!=0) verbose = false;

  int NumMyEquations = 10000;

  int NumGlobalEquations = NumMyEquations*NumProc;
  int NumGlobalVariables = 2 * NumGlobalEquations+1;

  // Construct a Map that puts approximately the same Number of equations on each processor

  Epetra_Map RowMap(NumGlobalEquations, 0, Comm);
  Epetra_Map XMap(NumGlobalVariables, 0, Comm);
  Epetra_Map& YMap = RowMap;
  
  int NumMyY = NumMyEquations;
  int IndexBase = 0;
  bool DistributedGlobal = (NumGlobalEquations>NumMyEquations);

  // Get update list and number of local equations from newly created Map
  int * MyGlobalElements = new int[RowMap.NumMyElements()];
  RowMap.MyGlobalElements(MyGlobalElements);

  // Get update list and number of local equations from newly created XMap
  int * XGlobalElements = new int[XMap.NumMyElements()];
  XMap.MyGlobalElements(XGlobalElements);

  // Get update list and number of local variables from newly created YMap
  int * YGlobalElements = new int[YMap.NumMyElements()];
  YMap.MyGlobalElements(YGlobalElements);

  // We need vectors to compute:
  // X = A^T*Y
  // AATY = A*A^T*Y = A*X
  //  and 
  // BY = B*Y

  Epetra_Vector Y(YMap);
  Epetra_Vector X(XMap);
  Epetra_Vector AATY(YMap);
  Epetra_Vector BY(YMap);


  // Fill Y Vector
  Y.Random();
  //Y.PutScalar(1.0);

  // To create A^T explicitly we need an assembly map that is two elements longer than
  // the XMap, because each processor will be making contributions to two rows beyond what
  // it will own.
  int ATAssemblyNumMyElements = 2*MyGlobalElements[NumMyEquations-1] + 2 - 2*MyGlobalElements[0] + 1;
  int * ATAssemblyGlobalElements = new int[ATAssemblyNumMyElements];

  for (i=0; i<ATAssemblyNumMyElements; i++) ATAssemblyGlobalElements[i] = 2*MyGlobalElements[0] + i;
  Epetra_Map ATAssemblyMap(-1, ATAssemblyNumMyElements, ATAssemblyGlobalElements, 0, Comm);

  // Create a Epetra_Matrix with the values of A
  // A is a simple 1D weighted average operator that mimics a restriction operator
  // that might be found in a multigrid code.
  // Also create A^T explicitly

  Epetra_CrsMatrix A(Copy, RowMap, 3);
  Epetra_CrsMatrix ATAssembly(Copy, ATAssemblyMap, 2);
  Epetra_CrsMatrix AT(Copy, XMap, 2);
  
  //cout << "ATAssemblyMap = "<< endl<< ATAssemblyMap << endl
  //     << "XMap = " << endl << XMap << endl
  //     << "RowMap = " << endl << RowMap << endl;
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  double *Values = new double[3];
  int *Indices = new int[3];
  int NumEntries;
  /*
  Values[0] = 0.25;
  Values[1] = 0.5;
  Values[2] = 0.25;
  */
  Values[0] = 0.5;
  Values[1] = 0.25;
  Values[2] = 0.25;
  forierr = 0;
  for (i=0; i<NumMyEquations; i++)
    {
  /*
      Indices[0] = 2*MyGlobalElements[i];
      Indices[1] = 2*MyGlobalElements[i]+1;
      Indices[2] = 2*MyGlobalElements[i]+2;
   */
      Indices[0] = 2*MyGlobalElements[i]+1;
      Indices[1] = 2*MyGlobalElements[i]+2;
      Indices[2] = 2*MyGlobalElements[i];
      NumEntries = 3;
      forierr += !(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
      for (j=0; j<3; j++)
	forierr += !(ATAssembly.InsertGlobalValues(Indices[j],1, &(Values[j]), &(MyGlobalElements[i]))>=0);
    }
  EPETRA_TEST_ERR(forierr,ierr);


  EPETRA_TEST_ERR(!(ATAssembly.TransformToLocal()==0),ierr);
  // Gather AT values from ATAssembly matrix
  Epetra_Export Exporter(ATAssemblyMap, XMap);
  EPETRA_TEST_ERR(!(AT.Export(ATAssembly, Exporter, Add)==0),ierr);

  // Finish up
  EPETRA_TEST_ERR(!(A.TransformToLocal(&XMap, &YMap)==0),ierr);
  EPETRA_TEST_ERR(!(AT.TransformToLocal(&YMap, &XMap)==0),ierr);


  if (verbose1 && NumGlobalEquations<20) { 
    if (verbose) cout << "\n\n Matrix A\n" << endl;
    cout << A << endl;
    if (verbose) cout << " \n\n Matrix A Transpose\n" << endl;
    cout <<  AT << endl;
  }


  // Create a Epetra_Matrix containing B = A*A^T.
  // This matrix will be a square tridiagonal matrix.  We will use it to compare the results
  // of A*(A^T*X) using two methods: (1) with two calls to Multiply using A^T and then A and
  // (2) using B directly.

  Epetra_CrsMatrix B(Copy, RowMap, 3);

  Values[0] = 1.0/16.0;
  Values[1] = 3.0/8.0;
  Values[2] = 1.0/16.0;
  int Valstart;
  forierr = 0;
  for (i=0; i<NumMyEquations; i++)
    {
      if (MyGlobalElements[i] == 0) {
      Indices[0] = MyGlobalElements[i];
      Indices[1] = MyGlobalElements[i]+1;
      NumEntries = 2;
      Valstart = 1;
      }
      else {
	Indices[0] = MyGlobalElements[i]-1;
	Indices[1] = MyGlobalElements[i];
	Indices[2] = MyGlobalElements[i]+1;
	NumEntries = 3;
	Valstart = 0;
      }
      if (MyGlobalElements[i] == NumGlobalEquations-1) NumEntries--;
      forierr += !(B.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values+Valstart, Indices)==0);
    }
  EPETRA_TEST_ERR(forierr,ierr);

  // Finish up
  EPETRA_TEST_ERR(!(B.TransformToLocal()==0),ierr);
  if (verbose && NumGlobalEquations<20) cout << "\n\nMatrix B \n" << endl;
  if (verbose1 && NumGlobalEquations<20) cout << B << endl;


  Epetra_Flops counter;
  A.SetFlopCounter(counter);
  B.SetFlopCounter(A);
  Epetra_Time timer(Comm);
  EPETRA_TEST_ERR(!(B.Multiply(false, Y, BY)==0),ierr); // Compute BY = B*Y
  double elapsed_time = timer.ElapsedTime();
  double total_flops = B.Flops();
  counter.ResetFlops();
  
  double MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for B*Y = " << MFLOPs << endl<< endl;
  if (verbose && NumGlobalEquations<20) cout << "\n\nVector Z = B*Y \n";
  if (verbose1 && NumGlobalEquations<20) cout << BY << endl;
 

  /////////////////////////////////////////////////////////////////////////////////////////////////

  timer.ResetStartTime();
  EPETRA_TEST_ERR(!(A.Multiply(true, Y, X)==0),ierr); // Compute X = A^T*Y
  elapsed_time = timer.ElapsedTime();
  total_flops = A.Flops();
  counter.ResetFlops();
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for A^T*Y using A and trans=true = " << MFLOPs << endl<< endl;
  if (verbose && NumGlobalEquations<20) cout << "\n\nVector Z = AT*Y \n";
  if (verbose1 && NumGlobalEquations<20) cout << X << endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  timer.ResetStartTime();
  EPETRA_TEST_ERR(!(A.Multiply(false, X, AATY)==0),ierr); // Compute AATY = A*X
  elapsed_time = timer.ElapsedTime();
  total_flops = A.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;
  counter.ResetFlops();
  Epetra_Vector resid(YMap);
  resid.Update(1.0, BY, -1.0, AATY, 0.0);
  double residual;
  resid.Norm2(&residual);

  if (verbose) cout << "\n\nTotal MFLOPs for A*X using A and trans=false = " << MFLOPs << endl<< endl;
  if (verbose) cout << "Residual = " << residual << endl<< endl;
  if (verbose && NumGlobalEquations<20) cout << "\n\nVector Z = A*ATY \n";
  if (verbose1 && NumGlobalEquations<20) cout << AATY << endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  AT.SetFlopCounter(counter);
  timer.ResetStartTime();
  EPETRA_TEST_ERR(!(AT.Multiply(false, Y, X)==0),ierr); // Compute X = A^T*Y
  elapsed_time = timer.ElapsedTime();
  total_flops = AT.Flops();
  counter.ResetFlops();
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for A^T*Y using AT and trans=false = " << MFLOPs << endl<< endl;
  if (verbose && NumGlobalEquations<20) cout << "\n\nVector Z = AT*Y \n";
  if (verbose1 && NumGlobalEquations<20) cout << X << endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  timer.ResetStartTime();
  EPETRA_TEST_ERR(!(AT.Multiply(true, X, AATY)==0),ierr); // Compute AATY = A*X
  elapsed_time = timer.ElapsedTime();
  total_flops = AT.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;
  counter.ResetFlops();
  resid.Update(1.0, BY, -1.0, AATY, 0.0);
  resid.Norm2(&residual);

  if (verbose) cout << "\n\nTotal MFLOPs for A*X using AT and trans=true = " << MFLOPs << endl<< endl;
  if (verbose) cout << "Residual = " << residual << endl<< endl;
  if (verbose && NumGlobalEquations<20) cout << "\n\nVector Z = A*ATY \n";
  if (verbose1 && NumGlobalEquations<20) cout <<AATY << endl;

  // Release all objects
  delete [] Values;
  delete [] Indices;
  delete [] MyGlobalElements;
  delete [] XGlobalElements;
  delete [] YGlobalElements;
  delete [] ATAssemblyGlobalElements;

			
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}

