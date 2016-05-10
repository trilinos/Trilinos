/*
//@HEADER
// ************************************************************************
// 
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
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
// Questions Contact  H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

// EXERCISE 1 Goal:
// Use Kokkos to parallelize the outer loop of <y,Ax> using Kokkos::parallel_reduce

#include "mtk_kokkos.h"

TEST_F(MTK_Kokkos, TUTORIAL_01) {

  int N = std::pow(2, 12);
  int M = std::pow(2, 10);
  int nrepeat = 100;

  // Allocate y, x vectors and Matrix A:
  double * const y = new double[ N ];
  double * const x = new double[ M ];
  double * const A = new double[ N * M ];

  // Initialize y vector
  // EXERCISE: convert outer loop to Kokkos::parallel_for
  for ( int i = 0 ; i < N ; ++i ) {
    y[i] = 1;
  }

  // Initialize x vector
  // EXERCISE: convert outer loop to Kokkos::parallel_for
  for ( int i = 0 ; i < M ; ++i ) {
    x[i] = 1;
  }

  // EXERCISE: convert outer loop to Kokkos::parallel_for
  for ( int j = 0 ; j < N ; ++j ) {
    for ( int i = 0 ; i < M ; ++i ) {
      A[ j * M + i ] = 1; 
    }
  }

  // Timer products
  struct timeval begin,end;

  gettimeofday(&begin,NULL);

  for ( int repeat = 0; repeat < nrepeat; repeat++) {

    //Application: <y,Ax> = y^T*A*x
    // EXERCISE: convert outer loop to Kokkos::parallel_reduce
    double result = 0;
    for ( int j = 0 ; j < N ; ++j ) {
      double temp2 = 0 ;
      for ( int i = 0 ; i < M ; ++i ) {
        temp2 += A[ j * M + i ] * x[ i ];
      }
      result += y[ j ] * temp2;
    }

    //Output result
    if ( repeat == (nrepeat - 1) )
      printf("  Computed result for %d x %d is %lf\n", N, M, result);
    const double solution = (double)N *(double)M;

    if ( result != solution ) {
      printf("  Error: result( %lf ) != solution( %lf )\n",result,solution);
    }
  }

  gettimeofday(&end,NULL);

  // Calculate time
  double time = 1.0*(end.tv_sec-begin.tv_sec) +
                1.0e-6*(end.tv_usec-begin.tv_usec);

  // Calculate bandwidth.
  // Each matrix A row (each of length M) is read once.
  // The x vector (of length M) is read N times.
  // The y vector (of length N) is read once.
  // double Gbytes = 1.0e-9 * double(sizeof(double) * ( 2 * M * N + N ));
  double Gbytes = 1.0e-9 * double(sizeof(double) * ( M + M * N + N )) ;

  // Print results (problem size, time and bandwidth in GB/s)
  printf("  M( %d ) N( %d ) nrepeat ( %d ) problem( %g MB ) time( %g s ) bandwidth( %g GB/s )\n",
         M , N, nrepeat, Gbytes * 1000, time, Gbytes * nrepeat / time );

  delete[] A ;
  delete[] y ;
  delete[] x ;

}
