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

#include "mtk_kokkos.h"
#include <stk_topology/topology.hpp>
  
void test_topo(int N, int M, int nrepeat) {

    // typedef Kokkos::Threads  ExecSpace ;

#ifdef KOKKOS_ENABLE_OPENMP
  typedef Kokkos::OpenMP   ExecSpace ;
#elif KOKKOS_ENABLE_CUDA
  typedef Kokkos::Cuda     ExecSpace ;
#else
  typedef Kokkos::Serial   ExecSpace ;
#endif

#ifdef KOKKOS_ENABLE_OPENMP
   typedef Kokkos::OpenMP       MemSpace;
#elif KOKKOS_ENABLE_CUDA
   typedef Kokkos::CudaSpace    MemSpace;
#else
   typedef Kokkos::HostSpace    MemSpace;
#endif

  // typedef Kokkos::CudaUVMSpace MemSpace; 

  typedef Kokkos::LayoutLeft   Layout ;
  // typedef Kokkos::LayoutRight  Layout ;

  typedef Kokkos::RangePolicy<ExecSpace> range_policy ;

  // Timer products
  struct timeval begin,end;

  gettimeofday(&begin,NULL);

  double result = 0;
  Kokkos::parallel_reduce(N, KOKKOS_LAMBDA(int i, double& update) {
    stk::topology hex8 = stk::topology::HEX_8;
    update += hex8.num_nodes();
  }, result);

  stk::topology hex8 = stk::topology::HEX_8;
  EXPECT_EQ(hex8.num_nodes()*N, result);

  gettimeofday(&end,NULL);

  // Calculate time
  double time = 1.0*(end.tv_sec-begin.tv_sec) +
                1.0e-6*(end.tv_usec-begin.tv_usec);

  // Print results (problem size, time and bandwidth in GB/s)
  printf("  M( %d ) N( %d ) nrepeat ( %d ) result( %g ) time(%g)\n", M , N, nrepeat, result, time);
}


TEST_F(MTK_Kokkos, topo_test) {
  int N = 5;
  int M = 10;
  int nrepeat = 1;
  test_topo(N, M, nrepeat);
}
