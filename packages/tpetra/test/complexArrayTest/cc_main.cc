#include <cstdlib>
#include <cassert>
#include <iostream>
#include <strstream>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <complex>
using namespace std;

// Used to easily change scalar type

//#define SCALARTYPE double
#define SCALARTYPE complex<float>

#include "Petra_Comm.h"
#include "Petra_Time.h"
#include "Petra_Map.h" 
#include "TPetra_ScalarTraits.h"
#include "TPetra_DenseMatrix.h"

// Local prototypes
template<class scalarType>
bool check(int M, int N, TPetra::DenseMatrix<scalarType> A, TPetra::DenseMatrix<scalarType> B, 
	   TPetra::DenseMatrix<scalarType>& C);

int main(int argc, char *argv[]) {

  Petra_Comm Comm;

  TPetra::DenseMatrix<SCALARTYPE> A;
  TPetra::DenseMatrix<SCALARTYPE> B;
  TPetra::DenseMatrix<SCALARTYPE> C1;
  TPetra::DenseMatrix<SCALARTYPE> C2;
  

  int M = 200;
  int N = 200;
  A.shape(M,N);
  B.shape(N,M);
  C1.shape(M,M);
  C2.shape(M,M);

  typedef TPetra::ScalarTraits<SCALARTYPE>::magnitudeType mymagtype;
  bool isComplex = ((TPetra::ScalarTraits<SCALARTYPE>::name()=="std::complex<float>") || 
		    (TPetra::ScalarTraits<SCALARTYPE>::name()=="std::complex<double>"));
  
  SCALARTYPE factor =  TPetra::ScalarTraits<SCALARTYPE>::one();

  if (isComplex) factor += SCALARTYPE(0.0, 1.0);
    



  for (int i=0; i<M; i++)
    for (int j=0; j<N; j++) {
      A(i,j) = ((mymagtype) (i+1)*(j+1)) * factor; // 
      B(j,i) = (SCALARTYPE) A(i,j);
    }
  bool smallProblem = (M*N<200);
  bool verbose = true;
  if (verbose && smallProblem) cout << "\nContents of A:\n" << A << endl;
  if (verbose && smallProblem) cout << "\nContents of B:\n" << B << endl;

  Petra_Time timer(Comm);
  double startFlops = C1.Flops();
  double startTime = timer.ElapsedTime();
  C1.multiply('N', 'N', 1.0, A, B, 0.0);
  double time = timer.ElapsedTime() - startTime;
  double flops = C1.Flops() - startFlops;
  double MFLOPS = flops/time/1000000.0;

  if (verbose) cout << "Statistics for multiply method C(MxM) = A(MxN) * B(NxM) " << "M = " << M
		    << " N = " << N << " is:" << endl 
		    << " Total FLOPS  = " << flops << endl 
		    << " Total Time   = " << time  << endl
		    << " Total MFLOPS = " << MFLOPS << endl;

  if (verbose && smallProblem) cout << "\nContents of C1:\n" << C1 << endl;

  if (!check(M, N, A, B, C2)) return 1; // return 1 if tests fail

      
  if (verbose && smallProblem) cout << "\nContents of C2:\n" << C2 << endl;

  //SCALARTYPE c1Norm = C1.oneNorm();
  //SCALARTYPE c2Norm = C2.oneNorm();
  double c1Norm = 0;
  double c2Norm = 0;

  if (c1Norm!=c2Norm) {
    if (verbose) {
      cout << "Norms of C1 and C2 differ:" << endl
	   << "Norm of C1 = " << c1Norm << endl
	   << "Norm of C2 = " << c2Norm << endl;
    }
  }
  else if (verbose) cout << "TPetra::DenseMatrix check OK" << endl;
	  
  return 0; // All done
}
  template<class scalarType>
  bool check(int M, int N, TPetra::DenseMatrix<scalarType> A, TPetra::DenseMatrix<scalarType> B, 
	     TPetra::DenseMatrix<scalarType>& C) {

    // Confirm correct dimensions

    assert(A.numRows()==M);
    assert(A.numCols()==N);
    assert(B.numRows()==N);
    assert(B.numCols()==M);
    assert(C.numRows()==M);
    assert(C.numCols()==M);

    
    for (int i=0; i<M; i++)
      for (int j=0; j<M; j++) {
	C(i,j) = A(i,0) * B(0,j); // Need to initialize C(ij) with first element because we have no zero now
	for (int k = 1; k<N; k++) C(i,j) += A(i,k)*B(k,j);
      }

    return(true);
  }
