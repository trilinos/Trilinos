#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#ifdef PETRA_MPI
#include "mpi.h"
#endif
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Petra_Comm.h"
#include "Petra_Time.h"
#include "Petra_BLAS_GE_Matrix.h"

#define TYPE double

 
int main(int argc, char *argv[])
{
  int ierr = 0, i, j, k;
  bool debug = false;

#ifdef PETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;



#ifdef PETRA_MPI
  Petra_Comm Comm(MPI_COMM_WORLD);
#else
  Petra_Comm Comm;
#endif


  //  char tmp;
  //  if (rank==0) cout << "Press any key to continue..."<< endl;
  //  if (rank==0) cin >> tmp;
  //  Comm.Barrier();

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  if (verbose) cout << Comm <<endl;

  bool verbose1 = verbose;

  // Redefine verbose to only print on PE 0
  if (verbose && rank!=0) verbose = false;

  Petra_BLAS_GE_Matrix<TYPE> A;
  Petra_BLAS_GE_Matrix<TYPE> X;
  Petra_BLAS_GE_Matrix<TYPE> B;
  Petra_BLAS_GE_Matrix<TYPE> B1;
  

  int N = 100;
  int NRHS = 5;
  if (argc==4) {
    N = atoi(argv[2]);
    NRHS = atoi(argv[3]);
  }

  if (verbose) cout << "Matrix dimension is " << N << "  Number of RHS = " << NRHS << endl;
  A.Shape(N,N);
  B.Shape(N,NRHS);
  B1.Shape(N,NRHS);
  X.Shape(N,NRHS);

  int izero = 0;
  int ione = 1;
  TYPE zero = izero;
  TYPE one = ione;
  for (int i = 0; i<N; i++) {
    TYPE xi = i%5;
    for (int j = 0; j < NRHS; j++) {
      X(i,j) = xi;
      B1(i,j) = zero;
    }
  }

  for (int i = 0; i<N; i++) {
    for (int j = 0; j < N; j++) {
      TYPE fi = i%7;
      TYPE fj = j%6;
      A(i,j) = fi * fj;
    }
  }
  for (int l=0;l<NRHS; l++)
    for (int i = 0; i<N; i++)
      for (int j = 0; j < N; j++)
      B1(i,l) += A(i,j)*X(j,l);

  Petra_Time Timer(Comm);
  double tstart = Timer.ElapsedTime();
  B.Multiply('N', 'N', one, A, X, zero);
  double time = Timer.ElapsedTime() - tstart;

  double FLOPS = B.Flops();
  double MFLOPS = FLOPS/time/1000000.0;
  if (verbose) cout << "MFLOPS for Matrix Multiplication = " << MFLOPS << endl;

  TYPE sum = zero;
  for (int j = 0; j<NRHS; j++)
    for (int i = 0; i<N; i++) sum += (B(i,j) - B1(i,j));

  if (verbose) cout << "Difference between computed and expected RHS = " << sum << endl;

#ifdef PETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}
