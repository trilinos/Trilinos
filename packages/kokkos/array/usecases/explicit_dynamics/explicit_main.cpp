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
  void test_Host(int beg, int end, int r, int threads);
  void test_Cuda(int beg, int end, int r);
}

int main(int argc, char ** argv)
{
  int beg = 4 ;
  int end = 12 ;
  int runs = 3 ;
  int host_threads = 4 ;

  if ( 1 < argc ) {
    host_threads = atoi(argv[1]);
  }
  if ( argc == 5) {
    beg = atoi(argv[2]);
    end = atoi(argv[3]);
    runs = atoi(argv[4]);
  }

std::cout << "\" " << argv[0]
                   << " host_threads begin end runs \"" << std::endl ;
std::cout << "\" " << argv[0]
          << " " << host_threads
          << " " << beg
          << " " << end
          << " " << runs
          << " \"" << std::endl ;


  if ( 0 <= host_threads ) {
    Test::test_Host(beg, end, runs, host_threads);
  }

#ifdef TEST_KOKKOSARRAY_CUDA
  Test::test_Cuda(beg , end, runs);
#endif

  return 0;
}
