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
#include "Petra_RDP_DenseMatrix.h"
#include "Petra_RDP_DenseVector.h"

// prototypes

int check(Petra_RDP_DenseMatrix &A, double * A1, int LDA,
	  int N1, int NRHS1, double OneNorm1, 
	  double * B1, int LDB1, 
	  double * X1, int LDX1,
	  bool Transpose, bool verbose);

void GenerateHilbert(double *A, int LDA, int N);

bool Residual( int N, int NRHS, double * A, int LDA, bool Transpose,
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
  bool Transpose = false;
  
  Petra_RDP_DenseMatrix * Matrix;
  for (int kk=0; kk<2; kk++) {
    for (i=1; i<=N; i++) {
      GenerateHilbert(A, LDA, i);
      OneNorm1 = 0.0;
      for (j=1; j<=i; j++) OneNorm1 += 1.0/((double) j); // 1-Norm = 1 + 1/2 + ...+1/n
      
      if (kk==0) {
	Matrix = new Petra_RDP_DenseMatrix(View, A, LDA, i, i);
	LDA1 = LDA;
      }
      else {
	Matrix = new Petra_RDP_DenseMatrix(Copy, A, LDA, i, i);
	LDA1 = i;
      }
      GenerateHilbert(A1, LDA1, i);
	
      if (kk==1) {
	Matrix->FactorWithEquilibration(true);
	Matrix->SolveWithTranspose(true);
	Transpose = true;
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
      
      ierr = check(*Matrix, A1, LDA1,  i, NRHS, OneNorm1, B1, LDB1,  X1, LDX1, Transpose, verbose);
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

  /////////////////////////////////////////////////////////////////////
  // Now test for larger system, both correctness and performance.
  /////////////////////////////////////////////////////////////////////


  N = 2000;
  NRHS = 5;
  LDA = N;
  LDB = N;
  LDX = N;

  if (verbose) cout << "\n\nComputing factor of an " << N << " x " << N << " general matrix...Please wait.\n\n" << endl;

  // Define A and X

  A = new double[LDA*N];
  X = new double[LDB*NRHS];
  
  for (j=0; j<N; j++) {
    for (k=0; k<NRHS; k++) X[j+k*LDX] = 1.0/((double) (j+5+k));
    for (i=0; i<N; i++) {
      if (i==((j+2)%N)) A[i+j*LDA] = 100.0 + i;
      else A[i+j*LDA] = -11.0/((double) (i+5)*(j+2));
    }
  }

  // Define Petra_RDP_DenseMatrix object

  Matrix = new Petra_RDP_DenseMatrix(Copy, A, LDA, N, N);

  // Time factorization

  Petra_Time & Timer = *new Petra_Time(Comm);
  double tstart = Timer.ElapsedTime();
  ierr = Matrix->Factor();
  if (ierr!=0 && verbose) cout << "Error in factorization = "<<ierr<< endl;
  assert(ierr==0);
  double time = Timer.ElapsedTime() - tstart;

  double FLOPS = Matrix->Flops();
  double MFLOPS = FLOPS/time/1000000.0;
  if (verbose) cout << "MFLOPS for Factorization = " << MFLOPS << endl;

  // Define Left hand side and right hand side 
  Petra_RDP_DenseMatrix &LHS = *new Petra_RDP_DenseMatrix(View, X, LDX, N, NRHS);
  Petra_RDP_DenseMatrix RHS;
  RHS.Shape(N,NRHS); // Allocate RHS

  // Compute RHS from A and X

  tstart = Timer.ElapsedTime();
  RHS.Multiply('N', 'N', 1.0, *Matrix, LHS, 0.0);
  time = Timer.ElapsedTime() - tstart;

  FLOPS = RHS.Flops();
  MFLOPS = FLOPS/time/1000000.0;
  if (verbose) cout << "MFLOPS to build RHS (NRHS = " << NRHS <<") = " << MFLOPS << endl;

  // Set LHS and RHS and solve
  Matrix->SetVectors(RHS, LHS);

  tstart = Timer.ElapsedTime();
  ierr = Matrix->Solve();
  if (ierr==1 && verbose) cout << "LAPACK guidelines suggest this matrix might benefit from equilibration." << endl;
  else if (ierr!=0 && verbose) cout << "Error in solve = "<<ierr<< endl;
  assert(ierr>=0);
  time = Timer.ElapsedTime() - tstart;

  FLOPS = Matrix->Flops();
  MFLOPS = FLOPS/time/1000000.0;
  if (verbose) cout << "MFLOPS for Solve (NRHS = " << NRHS <<") = " << MFLOPS << endl;

  double * resid = new double[NRHS];
  bool OK = Residual(N, NRHS, A, LDA, Matrix->Transpose(), Matrix->X(), Matrix->LDX(), Matrix->B(), Matrix->LDB(), resid);

  if (verbose) {
    if (!OK) cout << "************* Residual do not meet tolerance *************" << endl;
    for (i=0; i<NRHS; i++)
      cout << "Residual[" << i <<"] = "<< resid[i] << endl;
    cout  << endl;
  }

  // Solve again using the Petra_RDP_DenseVector class for LHS and RHS

  Petra_RDP_DenseVector X2;
  Petra_RDP_DenseVector B2;
  X2.Size(Matrix->N());
  B2.Size(Matrix->M());
  for (int kk=0; kk<Matrix->N(); kk++) X2[kk] = kk; // Define entries of X2

  tstart = Timer.ElapsedTime();
  B2.Multiply('N', 'N', 1.0, *Matrix, X2, 0.0); // Define B2 = A*X2
  time = Timer.ElapsedTime() - tstart;

  FLOPS = B2.Flops();
  MFLOPS = FLOPS/time/1000000.0;
  if (verbose) cout << "MFLOPS to build single RHS = " << MFLOPS << endl;

  // Set LHS and RHS and solve
  Matrix->SetVectors(B2, X2);

  tstart = Timer.ElapsedTime();
  ierr = Matrix->Solve();
  if (ierr==1 && verbose) cout << "LAPACK guidelines suggest this matrix might benefit from equilibration." << endl;
  else if (ierr!=0 && verbose) cout << "Error in solve = "<<ierr<< endl;
  assert(ierr>=0);
  time = Timer.ElapsedTime() - tstart;

  FLOPS = Matrix->Flops();
  MFLOPS = FLOPS/time/1000000.0;
  if (verbose) cout << "MFLOPS to solve single RHS = " << MFLOPS << endl;

  OK = Residual(N, 1, A, LDA, Matrix->Transpose(), Matrix->X(), Matrix->LDX(), Matrix->B(), Matrix->LDB(), resid);

  if (verbose) {
    if (!OK) cout << "************* Residual do not meet tolerance *************" << endl;
      cout << "Residual = "<< resid[0] << endl;
  }
  delete [] resid;
  delete [] A;
  delete [] X;
  delete Matrix;
  delete &LHS;

  ///////////////////////////////////////////////////
  // Now test default constructor and index operators
  ///////////////////////////////////////////////////

  N = 5;
  Petra_RDP_DenseMatrix C; // Implicit call to default constructor, should not need to call destructor
  C.Shape(5,5); // Make it 5 by 5
  double * C1 = new double[N*N];
  GenerateHilbert(C1, N, N); // Generate Hilber matrix

  C1[1+2*N] = 1000.0;  // Make matrix nonsymmetric

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
	 << "Values should be 1/(i+j+1), except value (1,2) should be 1000" << endl;

  
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

int check(Petra_RDP_DenseMatrix & A, double * A1, int LDA1, 
	  int N1, int NRHS1, double OneNorm1, 
	  double * B1, int LDB1, 
	  double * X1, int LDX1,
	  bool Transpose, bool verbose) {  

  int i, j;
  bool OK;
  // Test query functions

  int M= A.M();
  if (verbose) cout << "\n\nNumber of Rows = " << M << endl<< endl;
  assert(M==N1);

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
    assert(A.ROWCND()==-1.0);
    assert(A.COLCND()==-1.0);
    assert(A.AMAX()==-1.0);
  }


  // Other binary tests

  assert(!A.Factored());
  assert(A.Transpose()==Transpose);
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
  assert(A.Transpose()==Transpose);
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
  OK = Residual(N, NRHS, A1, LDA1, A.Transpose(), A.X(), A.LDX(), B1, LDB1, resid);
  if (verbose) {
    if (!OK) cout << "************* Residual do not meet tolerance *************" << endl;
 /*
    if (A.A_Equilibrated()) {
      double * R = A.R();
      double * C = A.C();
      for (i=0; i<A.M(); i++) 
      cout << "R[" << i <<"] = "<< R[i] << endl;
      for (i=0; i<A.N(); i++) 
      cout << "C[" << i <<"] = "<< C[i] << endl;
    }
 */
    cout << "\n\nResiduals using factorization to solve" << endl;
    for (i=0; i<NRHS; i++)
      cout << "Residual[" << i <<"] = "<< resid[i] << endl;
    cout  << endl;
  }


  ierr = A.Invert();
  assert(ierr>-1);

  assert(A.Inverted());
  assert(!A.Factored());
  assert(A.Transpose()==Transpose);

  
  Petra_RDP_DenseMatrix &RHS1 = *new Petra_RDP_DenseMatrix(Copy, B1, LDB1, N, NRHS);
  Petra_RDP_DenseMatrix &LHS1 = *new Petra_RDP_DenseMatrix(Copy, X1, LDX1, N, NRHS);
  assert(A.SetVectors(RHS1, LHS1)==0);
  assert(!A.Solved());

  assert(A.Solve()>-1);
	 
  
  A.Unequilibrate_X();

  OK = Residual(N, NRHS, A1, LDA1, A.Transpose(), A.X(), A.LDX(), B1, LDB1, resid);

  if (verbose) {
    if (!OK) cout << "************* Residual do not meet tolerance *************" << endl;
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

bool Residual( int N, int NRHS, double * A, int LDA, bool Transpose,
	       double * X, int LDX, double * B, int LDB, double * resid) {

  Petra_BLAS Blas;
  char Transa = 'N';
  if (Transpose) Transa = 'T';
  Blas.GEMM(Transa, 'N', N, NRHS, N, -1.0, A, LDA,
	    X, LDX, 1.0, B, LDB);
  bool OK = true;
  for (int i=0; i<NRHS; i++) {
    resid[i] = Blas.NRM2(N, B+i*LDB);
    if (resid[i]>1.0E-7) OK = false;
  }

  return(OK);
}
