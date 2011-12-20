/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
*/

#include <iostream>
#include <cstdlib>

namespace Test{
  void test_Host (const int beg, const int end, const int runs, const int num_iters);
  void test_TPI (const int beg, const int end, const int runs, const int num_iters, const int threads);
  void test_Pthread (const int beg, const int end, const int runs, const int num_iters, const int threads);
  void test_TBB (const int beg, const int end, const int runs, const int num_iters, const int threads);
  void test_Cuda (const int beg, const int end, const int runs, const int num_iters);
}

int 
main (int argc, char ** argv)
{
  int beg = 10;
  int end = 15;
  int runs = 3;
  int num_iters = 30;
  int threads = 0; // Let the library guess

  //
  // Command-line arguments:
  // beg: whatever Kurtis thinks that is, I have no idea
  // end: ditto
  // runs: Number of trials (entire runs, from assembly to solve)
  // num_iters: Number of (GMRES) solver iterations
  // threads: Number of threads (only applies to Pthread, TPI, and TBB devices)
  if (argc > 1) {
    beg = atoi (argv[1]);
  }    
  if (argc > 2) {
    end = atoi (argv[2]);
  }    
  if (argc > 3) {
    runs = atoi (argv[3]);
  }    
  if (argc > 4) {
    num_iters = atoi (argv[4]);
  }    
  if (argc > 5) {
    threads = atoi (argv[5]);
  }    

#ifdef TEST_KOKKOS_HOST
  Test::test_Host (beg, end, runs, num_iters);
#endif
#ifdef TEST_KOKKOS_PTHREAD
  Test::test_Pthread (beg, end, runs, num_iters, threads);
#endif
#ifdef TEST_KOKKOS_TPI
  Test::test_TPI (beg, end, runs, num_iters, threads);
#endif
#ifdef TEST_KOKKOS_TBB
  Test::test_TBB (beg, end, runs, num_iters, threads);
#endif
#ifdef TEST_KOKKOS_CUDA
  Test::test_Cuda (beg , end, runs, num_iters);
#endif

  return 0;
}
