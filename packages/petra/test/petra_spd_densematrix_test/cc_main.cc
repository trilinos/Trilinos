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
#include "Petra_Map.h"
#include "Petra_Time.h"
#include "Petra_RDP_SPD_DenseMatrix.h"

// prototypes

int check(Petra_RDP_SPD_DenseMatrix &A, double * A1, int LDA1, 
	  int N1, int NRHS1, double OneNorm1, 
	  double * B1, int LDB1, 
	  double * X1, int LDX1,
	  bool Upper, bool verbose);

void GenerateHilbert(double *A, int LDA, int N);

bool Residual( int N, int NRHS, double * A, int LDA,
	       double * X, int LDX, double * B, int LDB, double * resid);

 
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
  Petra_Comm & Comm = *new Petra_Comm( MPI_COMM_WORLD );
#else
  Petra_Comm & Comm = *new Petra_Comm();
#endif


  //  char tmp;
  //  if (rank==0) cout << "Press any key to continue..."<< endl;
  //  if (rank==0) cin >> tmp;
  //  Comm.Barrier();

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
              << " is alive."<<endl;

  bool verbose1 = verbose;

  // Redefine verbose to only print on PE 0
  if (verbose && rank!=0) verbose = false;

  int N = 20;
  int NRHS = 4;
  double * A = new double[N*N];
  double * A1 = new double[N*N];
  double * X = new double[(N+1)*NRHS];
  double * X1 = new double[(N+1)*NRHS];
  int LDX = N+1;
  int LDX1 = N+1;
  double * B = new double[N*NRHS];
  double * B1 = new double[N*NRHS];
  int LDB = N;
  int LDB1 = N;

  int LDA = N;
  int LDA1 = LDA;
  double OneNorm1;
  bool Upper = false;
  
  Petra_RDP_SPD_DenseMatrix * Matrix;
  for (int kk=0; kk<2; kk++) {
    for (i=1; i<=N; i++) {
      GenerateHilbert(A, LDA, i);
      OneNorm1 = 0.0;
      for (j=1; j<=i; j++) OneNorm1 += 1.0/((double) j); // 1-Norm = 1 + 1/2 + ...+1/n
      
      if (kk==0) {
	Matrix = new Petra_RDP_SPD_DenseMatrix(View, A, LDA, i);
	LDA1 = LDA;
      }
      else {
	Matrix = new Petra_RDP_SPD_DenseMatrix(Copy, A, LDA, i);
	LDA1 = i;
      }
      GenerateHilbert(A1, LDA1, i);
	
      if (kk==1) {
	Matrix->FactorWithEquilibration(true);
	Matrix->FactorUpperTriangle(true);
	Upper = true;
	Matrix->SolveToRefinedSolution(true);
      }
      
      for (k=0; k<NRHS; k++)
	for (j=0; j<i; j++) {
	  B[j+k*LDB] = 1.0/((double) (k+3)*(j+3));
	  B1[j+k*LDB1] = B[j+k*LDB1];
	}
      
      Petra_RDP_DenseMatrix & Petra_B = *new Petra_RDP_DenseMatrix(View, B, LDB, i, NRHS);
      Petra_RDP_DenseMatrix & Petra_X = *new Petra_RDP_DenseMatrix(View, X, LDX, i, NRHS);
      Matrix->SetVectors(Petra_B, Petra_X);
      
      int ierr = check(*Matrix, A1, LDA1,  i, NRHS, OneNorm1, B1, LDB1,  X1, LDX1, Upper, verbose);
      assert (ierr>-1);
      delete Matrix;
      delete &Petra_X;
      delete &Petra_B;
      if (ierr!=0) {
	if (verbose) cout << "Factorization failed due to bad conditioning.  This is normal if RCOND is small." 
			  << endl;
	break;
      }
    }
  }

  delete [] A;
  delete [] A1;
  delete [] X;
  delete [] X1;
  delete [] B;
  delete [] B1;

  // Now test for larger system, both correctness and performance.


  N = 2000;
  NRHS = 5;
  LDA = N;
  LDB = N;
  LDX = N;

  if (verbose) cout << "\n\nComputing factor of an " << N << " x " << N << " SPD matrix...Please wait.\n\n" << endl;

  A = new double[LDA*N];
  B = new double[LDB*NRHS];
  X = new double[LDB*NRHS];
  
  for (i=0; i<N; i++) B[i] = 0.0;
  for (j=0; j<N; j++) {
    for (k=0; k<NRHS; k++) X[j+k*LDX] = 1.0/((double) (j+5+k));
    for (i=0; i<N; i++) {
      if (i==j) A[i+j*LDA] = 100.0 + i;
      else A[i+j*LDA] = -1.0/((double) (i+10)*(j+10));
      for (k=0; k<NRHS; k++) B[i+k*LDB] += A[i+j*LDA]*X[j+k*LDX];
    }
  }
  Matrix = new Petra_RDP_SPD_DenseMatrix(Copy, A, LDA, N);

  Petra_Time & Timer = *new Petra_Time(Comm);
  double tstart = Timer.ElapsedTime();
  assert(Matrix->Factor()==0);
  double time = Timer.ElapsedTime() - tstart;

  double FLOPS = Matrix->Flops();
  double MFLOPS = FLOPS/time/1000000.0;
  if (verbose) cout << "MFLOPS for Factorization = " << MFLOPS << endl;

  Petra_RDP_DenseMatrix &RHS = *new Petra_RDP_DenseMatrix(View, B, LDB, N, NRHS);
  Petra_RDP_DenseMatrix &LHS = *new Petra_RDP_DenseMatrix(View, X, LDX, N, NRHS);
  Matrix->SetVectors(RHS, LHS);

  tstart = Timer.ElapsedTime();
  assert(Matrix->Solve()==0);
  time = Timer.ElapsedTime() - tstart;

  FLOPS = Matrix->Flops();
  MFLOPS = FLOPS/time/1000000.0;
  if (verbose) cout << "MFLOPS for Solve (NRHS = " << NRHS <<") = " << MFLOPS << endl;

  double * resid = new double[NRHS];
  assert(Residual(N, NRHS, A, LDA, Matrix->X(), Matrix->LDX(), Matrix->B(), Matrix->LDB(), resid));

  if (verbose) {
    for (i=0; i<NRHS; i++)
      cout << "Residual[" << i <<"] = "<< resid[i] << endl;
    cout  << endl;
  }
  delete [] resid;
  delete [] A;
  delete [] X;
  delete [] B;
  delete Matrix;
  delete &LHS;
  delete &RHS;


  ///////////////////////////////////////////////////
  // Now test default constructor and index operators
  ///////////////////////////////////////////////////

  N = 5;
  Petra_RDP_SPD_DenseMatrix C; // Implicit call to default constructor, should not need to call destructor
  C.Shape(5); // Make it 5 by 5
  double * C1 = new double[N*N];
  GenerateHilbert(C1, N, N); // Generate Hilbert matrix

  // Fill values of C with Hilbert values
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      C(i,j) = C1[i+j*N];

  // Test if values are correctly written and read
  for (i=0; i<N; i++)
    for (j=0; j<N; j++) {
      assert(C(i,j) == C1[i+j*N]);
      assert(C(i,j) == C[j][i]);
    }

  if (verbose)
    cout << "Default constructor and index operator check OK.  Values of Hilbert matrix = "
      << endl << C << endl
      << "Values should be 1/(i+j+1)" << endl;


  delete [] C1;
  delete &Timer;

  delete &Comm;

#ifdef PETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}

int check(Petra_RDP_SPD_DenseMatrix & A, double * A1, int LDA1, 
	  int N1, int NRHS1, double OneNorm1, 
	  double * B1, int LDB1, 
	  double * X1, int LDX1,
	  bool Upper, bool verbose) {  

  int i, j;
  // Test query functions

  int N= A.N();
  if (verbose) cout << "\n\nNumber of Equations = " << N << endl<< endl;
  assert(N==N1);

  int LDA = A.LDA();
  if (verbose) cout << "\n\nLDA = " << LDA << endl<< endl;
  assert(LDA==LDA1);

  int LDB = A.LDB();
  if (verbose) cout << "\n\nLDB = " << LDB << endl<< endl;
  assert(LDB==LDB1);

  int LDX = A.LDX();
  if (verbose) cout << "\n\nLDX = " << LDX << endl<< endl;
  assert(LDX==LDX1);

  int NRHS = A.NRHS();
  if (verbose) cout << "\n\nNRHS = " << NRHS << endl<< endl;
  assert(NRHS==NRHS1);

  assert(A.ANORM()==-1.0);
  assert(A.RCOND()==-1.0);
  if (!A.A_Equilibrated() && !A.B_Equilibrated()) {
    assert(A.SCOND()==-1.0);
    assert(A.AMAX()==-1.0);
  }


  // Other binary tests

  assert(!A.Factored());
  assert(A.Upper()==Upper);
  assert(!A.SolutionErrorsEstimated());
  assert(!A.Inverted());
  assert(!A.ReciprocalConditionEstimated());
  assert(!A.Solved());
  assert(!A.SolutionRefined());

      
  
  int ierr = A.Factor();
  assert(ierr>-1);
  if (ierr!=0) return(ierr); // Factorization failed due to poor conditioning.
  double rcond;
  assert(A.ReciprocalConditionEstimate(rcond)==0);
  if (verbose) {
    
    double rcond1 = 1.0/exp(3.5*((double)N));
    if (N==1) rcond1 = 1.0;
    cout << "\n\nRCOND = "<< rcond << " should be approx = " 
		    << rcond1 << endl << endl;
  }
  
  ierr = A.Solve();
  assert(ierr>-1);
  if (ierr!=0 && verbose) cout << "LAPACK rules suggest system should be equilibrated." << endl;

  assert(A.Factored());
  assert(A.Upper()==Upper);
  assert(A.ReciprocalConditionEstimated());
  assert(A.Solved());

  if (A.SolutionErrorsEstimated()) {
    if (verbose) {
      cout << "\n\nFERR[0] = "<< A.FERR()[0] << endl;
      cout << "\n\nBERR[0] = "<< A.BERR()[0] << endl<< endl;
    }
  }
  
  A.Unequilibrate_X();
  double * resid = new double[NRHS];
  assert(Residual(N, NRHS, A1, LDA1, A.X(), A.LDX(), B1, LDB1, resid));
  if (verbose) {
    cout << "Residuals using factorization to solve" << endl;
    for (i=0; i<NRHS; i++)
      cout << "Residual[" << i <<"] = "<< resid[i] << endl;
    cout  << endl;
  }


  ierr = A.Invert();
  assert(ierr>-1);

  assert(A.Inverted());
  assert(!A.Factored());
  assert(A.Upper()==Upper);

  Petra_RDP_DenseMatrix &RHS1 = *new Petra_RDP_DenseMatrix(Copy, B1, LDB1, N, NRHS);
  Petra_RDP_DenseMatrix &LHS1 = *new Petra_RDP_DenseMatrix(Copy, X1, LDX1, N, NRHS);
  assert(A.SetVectors(RHS1, LHS1)==0);
  assert(!A.Solved());

  assert(A.Solve()>-1);
	 
  
  A.Unequilibrate_X();

  assert(Residual(N, NRHS, A1, LDA1, A.X(), A.LDX(), B1, LDB1, resid));

  if (verbose) {
    cout << "Residuals using inverse to solve" << endl;
    for (i=0; i<NRHS; i++)
      cout << "Residual[" << i <<"] = "<< resid[i] << endl;
    cout  << endl;
  }
  delete [] resid;
  delete &RHS1;
  delete &LHS1;

  return(0);
}

 void GenerateHilbert(double *A, int LDA, int N) {
   for (int j=0; j<N; j++)
     for (int i=0; i<N; i++)
       A[i+j*LDA] = 1.0/((double)(i+j+1));
   return;
 }

bool Residual( int N, int NRHS, double * A, int LDA,
	       double * X, int LDX, double * B, int LDB, double * resid) {

  Petra_BLAS & Blas = *new Petra_BLAS();
  Blas.GEMM('N', 'N', N, NRHS, N, -1.0, A, LDA,
	    X, LDX, 1.0, B, LDB);
  bool OK = true;
  for (int i=0; i<NRHS; i++) {
    resid[i] = Blas.NRM2(N, B+i*LDB);
    if (resid[i]>1.0E-8) OK = false;
  }

  delete &Blas;
  return(OK);
}
