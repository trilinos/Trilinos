#include <cstdlib>
#include <cassert>
#include <iostream>
#include <strstream>
#include <cstring>
#include <cmath>
#include <iomanip>
using namespace std;

#include "Petra_Comm.h"
#include "Petra_Map.h"

#include "TPetra_DenseMatrix.h"
#include "WhizBangMatrixMultiplier.h"

#define FLOAT double
int main(int argc, char *argv[]) {

  Petra_Comm Comm;
  Petra_TSF_DenseMatrix<FLOAT> A;
  Petra_TSF_DenseMatrix<FLOAT> B;
  Petra_TSF_DenseMatrix<FLOAT> C1;
  Petra_TSF_DenseMatrix<FLOAT> C2;
  

  int M = 5;
  int N = 6;
  A.shape(M,N);
  B.shape(N,M);
  C1.shape(M,M);
  C2.shape(M,M);

  for (int i=0; i<M; i++)
    for (int j=0; j<N; j++) {
      A(i,j) = (FLOAT) i*j;
      B(i,j) = (FLOAT) i*j;


  cout << "\nContents of A:\n" << A << endl;
  cout << "\nContents of B:\n" << B << endl;

  
  C1.Multiply('N', 'T', 1.0, A, B, 0.0);

  cout << "\nContents of C1:\n" << C1 << endl;

  WhizBangMatrixMultiplier<FLOAT> multiplier(&A, &B, &C2);

  multiplier.multiply(); // Compute C2 = A*B

  cout << "\nContents of C2:\n" << C2 << endl;

  return 0; // All done
}
