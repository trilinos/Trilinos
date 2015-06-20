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
#include "Epetra_Version.h"

// prototypes

int check(Ifpack_SerialTriDiSolver & solver, double * A1, int LDA,
	  int N1, int NRHS1, double OneNorm1,
	  double * B1, int LDB1,
	  double * X1, int LDX1,
	  bool Transpose, bool verbose);



bool Residual( int N, int NRHS, double * A, int LDA, bool Transpose,
	       double * X, int LDX, double * B, int LDB, double * resid);

int matrixCpyCtr(bool verbose, bool debug);

void printHeading(const char* heading);
void printMat(const char* name, Ifpack_SerialTriDiMatrix& matrix);
void printArray(double* array, int length);

using namespace std;

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  bool verbose = false;

  // Check if we should print results to standard out
  verbose = true;

  if (verbose && Comm.MyPID()==0)
    cout << Epetra_Version() << endl << endl;

  int rank = Comm.MyPID();

  if (verbose) cout << Comm <<endl;

  // Redefine verbose to only print on PE 0
  if (verbose && rank!=0) verbose = false;
	
  int N = 5;
  int NRHS = 5;
  double * X = new double[NRHS];
  double * ed_X = new double[NRHS];

  double * B = new double[NRHS];
  double * ed_B = new double[NRHS];

  Ifpack_SerialTriDiSolver solver;
  Ifpack_SerialTriDiMatrix * Matrix;

  Epetra_SerialDenseSolver ed_solver;
  Epetra_SerialDenseMatrix * ed_Matrix;

  bool Transpose = false;
  bool Refine = false;
  solver.SolveWithTranspose(Transpose);
  solver.SolveToRefinedSolution(Refine);

  ed_solver.SolveWithTranspose(Transpose);
  ed_solver.SolveToRefinedSolution(Refine);

  Matrix = new Ifpack_SerialTriDiMatrix(5,true);
  ed_Matrix = new Epetra_SerialDenseMatrix(5,5);

  for(int i=0;i<N;++i) {
    B[i] = ed_B[i] =2;
    Matrix->D()[i]=2.0;
    if(i<(N-1)) {
      Matrix->DL()[i]=-1.0;
      if(i!=2) Matrix->DU()[i]=-1.0;
    }
  }

  Matrix->Print(std::cout);

  double * ed_a = ed_Matrix->A();
  for(int i=0;i<N;++i)
    for(int j=0;j<N;++j) {
      if(i==j) ed_a[j*N+i] = 2.0;
      else if(abs(i-j) == 1)   ed_a[j*N+i] = -1.0;
      else  ed_a[j*N + i] = 0;
      if(i==2&&j==3) ed_a[j*N+i] = 0.0;
    }


  Epetra_SerialDenseVector LHS(Copy, X, N);
  Epetra_SerialDenseVector RHS(Copy, B, N);

  Epetra_SerialDenseVector ed_LHS(Copy, ed_X, N);
  Epetra_SerialDenseVector ed_RHS(Copy, ed_B, N);

  solver.SetMatrix(*Matrix);
  solver.SetVectors(LHS, RHS);
  
  ed_solver.SetMatrix(*ed_Matrix);
  ed_solver.SetVectors(ed_LHS, ed_RHS);

  solver.Solve();  
  ed_solver.Solve();

  std::cout << " LHS vals are: "<<std::endl;
  bool TestPassed=true;
  for(int i=0;i<N;++i) { 
    std::cout << "["<<i<<"] "<< LHS(i)<<"  "<<ed_LHS(i)<<" delta "<<LHS(i)-ed_LHS(i)<<std::endl;
    if( fabs( (LHS(i)- ed_LHS(i))/(LHS(i)+ ed_LHS(i)) ) > 1.0e-12 ) {
       TestPassed = false;
       std::cout << " not equal for "<<i<<" delta "<< LHS(i)- ed_LHS(i)<<std::endl;
    }
  }

  Ifpack_SerialTriDiMatrix * tdfac = solver.FactoredMatrix();
  Epetra_SerialDenseMatrix * sdfac = ed_solver.FactoredMatrix();

  std::cout << " Tri Di Factored "<<std::endl;
  tdfac->Print(std::cout);
  std::cout << " Dense Factored "<<std::endl;
  sdfac->Print(std::cout);

  delete Matrix;
  delete ed_Matrix;
  delete [] X;
  delete [] ed_X;
  delete [] B;
  delete [] ed_B;


  if (!TestPassed) {
    cout << "Test `TestRelaxation.exe' failed!" << endl;
    exit(EXIT_FAILURE);
  }
  
#ifdef HAVE_MPI
  MPI_Finalize(); 
#endif

  cout << endl;
  cout << "Test `TestRelaxation.exe' passed!" << endl;
  cout << endl;
  return(EXIT_SUCCESS);
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

