//@HEADER
// ************************************************************************
//
//               Epetra: Linear Algebra Services Package
//                 Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER


#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Ifpack_SerialTriDiMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Ifpack_SerialTriDiSolver.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif
#include "Epetra_SerialComm.h"
#include "../../epetra/test/epetra_test_err.h"
#include "Epetra_Version.h"

// prototypes

int check(Ifpack_SerialTriDiSolver & solver, double * A1, int LDA,
	  int N1, int NRHS1, double OneNorm1,
	  double * B1, int LDB1,
	  double * X1, int LDX1,
	  bool Transpose, bool verbose);

void GenerateHilbert(double *A, int LDA, int N);

bool Residual( int N, int NRHS, double * A, int LDA, bool Transpose,
	       double * X, int LDX, double * B, int LDB, double * resid);

int matrixCpyCtr(bool verbose, bool debug);

void printHeading(const char* heading);
double* getRandArray(int length);
void printMat(const char* name, Ifpack_SerialTriDiMatrix& matrix);
void printArray(double* array, int length);

int main(int argc, char *argv[])
{
  int ierr = 0;

// #ifdef EPETRA_MPI
//   MPI_Init(&argc,&argv);
//   Epetra_MpiComm Comm( MPI_COMM_WORLD );
// #else
  Epetra_SerialComm Comm;
// #endif

  bool verbose = false;

  // Check if we should print results to standard out
  verbose = true;

  if (verbose && Comm.MyPID()==0)
    cout << Epetra_Version() << endl << endl;

  int rank = Comm.MyPID();
  //  char tmp;
  //  if (rank==0) cout << "Press any key to continue..."<< endl;
  //  if (rank==0) cin >> tmp;
  //  Comm.Barrier();

  //  Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  if (verbose) cout << Comm <<endl;

  //  bool verbose1 = verbose;

  // Redefine verbose to only print on PE 0
  if (verbose && rank!=0) verbose = false;
	
  int N = 4;
  int NRHS = 4;
  double * X = new double[NRHS];
  double * ed_X = new double[NRHS];

  double * B = new double[NRHS];
  double * ed_B = new double[NRHS];


  Ifpack_SerialTriDiSolver solver;
  Ifpack_SerialTriDiMatrix * Matrix;

  Epetra_SerialDenseSolver ed_solver;
  Epetra_SerialDenseMatrix * ed_Matrix;

  //  solver.FactorWithEquilibration(false);

  bool Transpose = false;
  bool Refine = false;
  solver.SolveWithTranspose(Transpose);
  solver.SolveToRefinedSolution(Refine);

  ed_solver.SolveWithTranspose(Transpose);
  ed_solver.SolveToRefinedSolution(Refine);

  Matrix = new Ifpack_SerialTriDiMatrix(4,true);
  ed_Matrix = new Epetra_SerialDenseMatrix(4,4);

  // memset(Matrix->A(),0,sizeof(double)*4*(N-1));

  for(int i=0;i<N;++i) {
    B[i] = ed_B[i] =0.04;
    Matrix->D()[i]=-2.0;
    if(i<(N-1)) {
      Matrix->DL()[i]=1.0;
      Matrix->DU()[i]=1.0;
    }
  }

  for(int i=0;i<N;++i) {
	if (i==0) std::cout << "DL  - ";
	else std::cout <<Matrix->DL()[i-1]<<" ";
    }
    std::cout<<std::endl;

    std::cout <<" D ";
    for(int i=0;i<N;++i) {
      std::cout << Matrix->D()[i]<<" ";
    }
    std::cout<<std::endl;

    std::cout <<" DU ";
    for(int i=0;i<N;++i) {
      if (i==0) std::cout << " - ";
      else std::cout << Matrix->DU()[i-1]<<" ";
    }
    std::cout<<std::endl;
    
    std::cout <<" DU2 ";
    for(int i=0;i<N;++i) {
      if (i<2) std::cout << " - ";
      else std::cout << Matrix->DU2()[i-2]<<" ";
    }
    std::cout<<std::endl;
    
  double * ed_a = ed_Matrix->A();
  for(int i=0;i<N;++i)
    for(int j=0;j<N;++j) {
      if(i==j) ed_a[j*N+i] = -2.0;
      else if(abs(i-j) == 1)   ed_a[j*N+i] = 1.0;
      else  ed_a[j*N + i] = 0;
    }


  Epetra_SerialDenseVector LHS(Copy, X, N);
  Epetra_SerialDenseVector RHS(Copy, B, N);

  Epetra_SerialDenseVector ed_LHS(Copy, ed_X, N);
  Epetra_SerialDenseVector ed_RHS(Copy, ed_B, N);

  solver.SetMatrix(*Matrix);
  solver.SetVectors(LHS, RHS);
  
  ed_solver.SetMatrix(*ed_Matrix);
  ed_solver.SetVectors(ed_LHS, ed_RHS);

  ierr = solver.Solve();
  std::cout << " tridi ierr "<<ierr<<std::endl;
  
  ierr = ed_solver.Solve();
  std::cout << " serialdensesolve ierr "<<ierr<<std::endl;

  std::cout << " LHS vals are: "<<std::endl;
  for(int i=0;i<N;++i) 
    std::cout << "["<<i<<"] "<< LHS(i)<<"  "<<ed_LHS(i)<<std::endl;


 //  delete [] A;
//   delete [] A1;
//   delete [] X;
//   delete [] X1;
//   delete [] B;
//   delete [] B1;

//   /////////////////////////////////////////////////////////////////////
//   // Now test norms and scaling functions
//   /////////////////////////////////////////////////////////////////////

//   Ifpack_SerialTriDiMatrix D;
//   double ScalarA = 2.0;

//   int DM = 10;
//   int DN = 8;
//   D.Shape(DM, DN);
//   for (j=0; j<DN; j++)
//     for (i=0; i<DM; i++) D[j][i] = (double) (1+i+j*DM) ;

//   //cout << D << endl;

//   double NormInfD_ref = (double)(DM*(DN*(DN+1))/2);
//   double NormOneD_ref = (double)((DM*DN*(DM*DN+1))/2 - (DM*(DN-1)*(DM*(DN-1)+1))/2 );

//   double NormInfD = D.NormInf();
//   double NormOneD = D.NormOne();

//   if (verbose) {
//     cout << " *** Before scaling *** " << endl
// 	 << " Computed one-norm of test matrix = " << NormOneD << endl
// 	 << " Expected one-norm                = " << NormOneD_ref << endl
// 	 << " Computed inf-norm of test matrix = " << NormInfD << endl
// 	 << " Expected inf-norm                = " << NormInfD_ref << endl;
//   }
//   D.Scale(ScalarA); // Scale entire D matrix by this value
//   NormInfD = D.NormInf();
//   NormOneD = D.NormOne();
//   if (verbose) {
//     cout << " *** After scaling *** " << endl
// 	 << " Computed one-norm of test matrix = " << NormOneD << endl
// 	 << " Expected one-norm                = " << NormOneD_ref*ScalarA << endl
// 	 << " Computed inf-norm of test matrix = " << NormInfD << endl
// 	 << " Expected inf-norm                = " << NormInfD_ref*ScalarA << endl;
//   }


//   /////////////////////////////////////////////////////////////////////
//   // Now test that A.Multiply(false, x, y) produces the same result
//   // as y.Multiply('N','N', 1.0, A, x, 0.0).
//   /////////////////////////////////////////////////////////////////////

//   N = 10;
//   int M = 10;
//   LDA = N;
//   Ifpack_SerialTriDiMatrix smallA(N, M, false);
//   Ifpack_SerialTriDiMatrix x(N, 1, false);
//   Ifpack_SerialTriDiMatrix y1(N, 1, false);
//   Ifpack_SerialTriDiMatrix y2(N, 1, false);

//   for(i=0; i<N; ++i) {
//     for(j=0; j<M; ++j) {
//       smallA(i,j) = 1.0*i+2.0*j+1.0;
//     }
//     x(i,0) = 1.0;
//     y1(i,0) = 0.0;
//     y2(i,0) = 0.0;
//   }

//   //quick check of operator==
//   if (x == y1) {
//     if (verbose) cout << "err in Ifpack_SerialTriDiMatrix::operator==, "
//         << "erroneously returned true." << std::endl;
//     return(-1);
//   }

//   //quick check of operator!=
//   if (x != x) {
//     if (verbose) cout << "err in Ifpack_SerialTriDiMatrix::operator==, "
//         << "erroneously returned true." << std::endl;
//     return(-1);
//   }

//   int err1 = smallA.Multiply(false, x, y1);
//   int err2 = y2.Multiply('N','N', 1.0, smallA, x, 0.0);
//   if (err1 != 0 || err2 != 0) {
//     if (verbose) cout << "err in Ifpack_SerialTriDiMatrix::Multiply"<<endl;
//     return(err1+err2);
//   }

//   for(i=0; i<N; ++i) {
//     if (y1(i,0) != y2(i,0)) {
//       if (verbose) cout << "different versions of Multiply don't match."<<endl;
//       return(-99);
//     }
//   }

//   /////////////////////////////////////////////////////////////////////
//   // Now test for larger system, both correctness and performance.
//   /////////////////////////////////////////////////////////////////////


//   N = 2000;
//   NRHS = 5;
//   LDA = N;
//   LDB = N;
//   LDX = N;

//   if (verbose) cout << "\n\nComputing factor of an " << N << " x " << N << " general matrix...Please wait.\n\n" << endl;

//   // Define A and X

//   A = new double[LDA*N];
//   X = new double[LDB*NRHS];

//   for (j=0; j<N; j++) {
//     for (k=0; k<NRHS; k++) X[j+k*LDX] = 1.0/((double) (j+5+k));
//     for (i=0; i<N; i++) {
//       if (i==((j+2)%N)) A[i+j*LDA] = 100.0 + i;
//       else A[i+j*LDA] = -11.0/((double) (i+5)*(j+2));
//     }
//   }

//   // Define Ifpack_SerialTriDiMatrix object

//   Ifpack_SerialTriDiMatrix BigMatrix(Copy, A, LDA, N, N);
//   Ifpack_SerialTriDiMatrix OrigBigMatrix(View, A, LDA, N, N);

//   Ifpack_SerialTriDiSolver BigSolver;
//   BigSolver.FactorWithEquilibration(true);
//   BigSolver.SetMatrix(BigMatrix);

//   // Time factorization

//   Epetra_Flops counter;
//   BigSolver.SetFlopCounter(counter);
//   Epetra_Time Timer(Comm);
//   double tstart = Timer.ElapsedTime();
//   ierr = BigSolver.Factor();
//   if (ierr!=0 && verbose) cout << "Error in factorization = "<<ierr<< endl;
//   assert(ierr==0);
//   double time = Timer.ElapsedTime() - tstart;

//   double FLOPS = counter.Flops();
//   double MFLOPS = FLOPS/time/1000000.0;
//   if (verbose) cout << "MFLOPS for Factorization = " << MFLOPS << endl;

//   // Define Left hand side and right hand side
//   Ifpack_SerialTriDiMatrix LHS(View, X, LDX, N, NRHS);
//   Ifpack_SerialTriDiMatrix RHS;
//   RHS.Shape(N,NRHS); // Allocate RHS

//   // Compute RHS from A and X

//   Epetra_Flops RHS_counter;
//   RHS.SetFlopCounter(RHS_counter);
//   tstart = Timer.ElapsedTime();
//   RHS.Multiply('N', 'N', 1.0, OrigBigMatrix, LHS, 0.0);
//   time = Timer.ElapsedTime() - tstart;

//   Ifpack_SerialTriDiMatrix OrigRHS = RHS;

//   FLOPS = RHS_counter.Flops();
//   MFLOPS = FLOPS/time/1000000.0;
//   if (verbose) cout << "MFLOPS to build RHS (NRHS = " << NRHS <<") = " << MFLOPS << endl;

//   // Set LHS and RHS and solve
//   BigSolver.SetVectors(LHS, RHS);

//   tstart = Timer.ElapsedTime();
//   ierr = BigSolver.Solve();
//   if (ierr==1 && verbose) cout << "LAPACK guidelines suggest this matrix might benefit from equilibration." << endl;
//   else if (ierr!=0 && verbose) cout << "Error in solve = "<<ierr<< endl;
//   assert(ierr>=0);
//   time = Timer.ElapsedTime() - tstart;

//   FLOPS = BigSolver.Flops();
//   MFLOPS = FLOPS/time/1000000.0;
//   if (verbose) cout << "MFLOPS for Solve (NRHS = " << NRHS <<") = " << MFLOPS << endl;

//   double * resid = new double[NRHS];
//   bool OK = Residual(N, NRHS, A, LDA, BigSolver.Transpose(), BigSolver.X(), BigSolver.LDX(),
// 		     OrigRHS.A(), OrigRHS.LDA(), resid);

//   if (verbose) {
//     if (!OK) cout << "************* Residual do not meet tolerance *************" << endl;
//     for (i=0; i<NRHS; i++)
//       cout << "Residual[" << i <<"] = "<< resid[i] << endl;
//     cout  << endl;
//   }

//   // Solve again using the Epetra_SerialDenseVector class for LHS and RHS

//   Epetra_SerialDenseVector X2;
//   Epetra_SerialDenseVector B2;
//   X2.Size(BigMatrix.N());
//   B2.Size(BigMatrix.M());
//   int length = BigMatrix.N();
//   {for (int kk=0; kk<length; kk++) X2[kk] = ((double ) kk)/ ((double) length);} // Define entries of X2

//   RHS_counter.ResetFlops();
//   B2.SetFlopCounter(RHS_counter);
//   tstart = Timer.ElapsedTime();
//   B2.Multiply('N', 'N', 1.0, OrigBigMatrix, X2, 0.0); // Define B2 = A*X2
//   time = Timer.ElapsedTime() - tstart;

//   Epetra_SerialDenseVector OrigB2 = B2;

//   FLOPS = RHS_counter.Flops();
//   MFLOPS = FLOPS/time/1000000.0;
//   if (verbose) cout << "MFLOPS to build single RHS = " << MFLOPS << endl;

//   // Set LHS and RHS and solve
//   BigSolver.SetVectors(X2, B2);

//   tstart = Timer.ElapsedTime();
//   ierr = BigSolver.Solve();
//   time = Timer.ElapsedTime() - tstart;
//   if (ierr==1 && verbose) cout << "LAPACK guidelines suggest this matrix might benefit from equilibration." << endl;
//   else if (ierr!=0 && verbose) cout << "Error in solve = "<<ierr<< endl;
//   assert(ierr>=0);

//   FLOPS = counter.Flops();
//   MFLOPS = FLOPS/time/1000000.0;
//   if (verbose) cout << "MFLOPS to solve single RHS = " << MFLOPS << endl;

//   OK = Residual(N, 1, A, LDA, BigSolver.Transpose(), BigSolver.X(), BigSolver.LDX(), OrigB2.A(),
// 		OrigB2.LDA(), resid);

//   if (verbose) {
//     if (!OK) cout << "************* Residual do not meet tolerance *************" << endl;
//       cout << "Residual = "<< resid[0] << endl;
//   }
//   delete [] resid;
//   delete [] A;
//   delete [] X;

//   ///////////////////////////////////////////////////
//   // Now test default constructor and index operators
//   ///////////////////////////////////////////////////

//   N = 5;
//   Ifpack_SerialTriDiMatrix C; // Implicit call to default constructor, should not need to call destructor
//   C.Shape(5,5); // Make it 5 by 5
//   double * C1 = new double[N*N];
//   GenerateHilbert(C1, N, N); // Generate Hilber matrix

//   C1[1+2*N] = 1000.0;  // Make matrix nonsymmetric

//   // Fill values of C with Hilbert values
//   for (i=0; i<N; i++)
//     for (j=0; j<N; j++)
//       C(i,j) = C1[i+j*N];

//   // Test if values are correctly written and read
//   for (i=0; i<N; i++)
//     for (j=0; j<N; j++) {
//       assert(C(i,j) == C1[i+j*N]);
//       assert(C(i,j) == C[j][i]);
//     }

//   if (verbose)
//     cout << "Default constructor and index operator check OK.  Values of Hilbert matrix = "
// 	 << endl << C << endl
// 	 << "Values should be 1/(i+j+1), except value (1,2) should be 1000" << endl;

//   delete [] C1;

//   // now test sized/shaped constructor
//   Ifpack_SerialTriDiMatrix shapedMatrix(10, 12);
//   assert(shapedMatrix.M() == 10);
//   assert(shapedMatrix.N() == 12);
//   for(i = 0; i < 10; i++)
//     for(j = 0; j < 12; j++)
//       assert(shapedMatrix(i, j) == 0.0);
//   Epetra_SerialDenseVector sizedVector(20);
//   assert(sizedVector.Length() == 20);
//   for(i = 0; i < 20; i++)
//     assert(sizedVector(i) == 0.0);
//   if (verbose)
//     cout << "Shaped/sized constructors check OK." << endl;

//   // test Copy/View mode in op= and cpy ctr
//   int temperr = 0;
//   temperr = matrixAssignment(verbose, debug);
//   if(verbose && temperr == 0)
//     cout << "Operator = checked OK." << endl;
//   EPETRA_TEST_ERR(temperr, ierr);
//   temperr = matrixCpyCtr(verbose, debug);
//   if(verbose && temperr == 0)
//     cout << "Copy ctr checked OK." << endl;
//   EPETRA_TEST_ERR(temperr, ierr);

//   // Test some vector methods

//   Epetra_SerialDenseVector v1(3);
//   v1[0] = 1.0;
//   v1[1] = 3.0;
//   v1[2] = 2.0;

//   Epetra_SerialDenseVector v2(3);
//   v2[0] = 2.0;
//   v2[1] = 1.0;
//   v2[2] = -2.0;

//   temperr = 0;
//   if (v1.Norm1()!=6.0) temperr++;
//   if (fabs(sqrt(14.0)-v1.Norm2())>1.0e-6) temperr++;
//   if (v1.NormInf()!=3.0) temperr++;
//   if(verbose && temperr == 0)
//     cout << "Vector Norms checked OK." << endl;
//   temperr = 0;
//   if (v1.Dot(v2)!=1.0) temperr++;
//   if(verbose && temperr == 0)
//     cout << "Vector Dot product checked OK." << endl;

// #ifdef EPETRA_MPI
//   MPI_Finalize() ;
// #endif

/* end main
*/
return ierr ;
}


 void GenerateHilbert(double *A, int LDA, int N) {
   for (int j=0; j<N; j++)
     for (int i=0; i<N; i++)
       A[i+j*LDA] = 1.0/((double)(i+j+1));
   return;
 }

bool Residual( int N, int NRHS, double * A, int LDA, bool Transpose,
	       double * X, int LDX, double * B, int LDB, double * resid) {

  Epetra_BLAS Blas;
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


//=========================================================================

//=========================================================================
//=========================================================================
// prints section heading with spacers/formatting
void printHeading(const char* heading) {
	cout << "\n==================================================================\n";
	cout << heading << endl;
	cout << "==================================================================\n";
}

//=========================================================================
// prints SerialTriDiMatrix/Vector with formatting
void printMat(const char* name, Ifpack_SerialTriDiMatrix& matrix) {
	//cout << "--------------------" << endl;
	cout << "*** " << name << " ***" << endl;
	cout << matrix;
	//cout << "--------------------" << endl;
}

//=========================================================================
// prints double* array with formatting
void printArray(double* array, int length) {
	cout << "user array (size " << length << "): ";
	for(int i = 0; i < length; i++)
		cout << array[i] << "  ";
	cout << endl;
}

//=========================================================================
// returns a double* array of a given length, with random values on interval (-1,1).
// this is the same generator used in SerialTriDiMatrix
double* getRandArray(int length) {
  const double a = 16807.0;
	const double BigInt = 2147483647.0;
	const double DbleOne = 1.0;
	const double DbleTwo = 2.0;
	double seed = rand();

	double* array = new double[length];

	for(int i = 0; i < length; i++) {
		seed = fmod(a * seed, BigInt);
		array[i] = DbleTwo * (seed / BigInt) - DbleOne;
	}

	return(array);
}

