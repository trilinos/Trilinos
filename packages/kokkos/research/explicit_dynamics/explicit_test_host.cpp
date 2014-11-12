/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <iostream>
#include <iomanip>
#include <sys/time.h>

#include <Kokkos_Value.hpp>
#include <Kokkos_MultiVector.hpp>
#include <Kokkos_MDArray.hpp>

#include <Kokkos_Core.hpp>

#include <explicit_dynamics_app.hpp>

namespace Test{

void test_Host( int beg, int end, int runs, int threads)
{
  const unsigned numa_count = Kokkos::hwloc::get_available_numa_count();
  const unsigned core_numa = Kokkos::hwloc::get_available_cores_per_numa();
  const unsigned thread_core = Kokkos::hwloc::get_available_threads_per_core();

  if ( 0 < threads ) {
    const size_t node_thread_count = ( threads + numa_count - 1 ) / numa_count ;

    Kokkos::Threads::initialize( numa_count * node_thread_count , numa_count );

    std::cout << std::endl << "\"Threads with manually set threads = \" , "
              << numa_count * node_thread_count << std::endl ;
  }
  else {
    Kokkos::Threads::initialize( numa_count * core_numa * thread_core , numa_count );

    std::cout << std::endl << "\"Threads with detected sequential threads = \" , "
              << numa_count * node_thread_count << std::endl ;
  }

  explicit_dynamics::driver<float,Kokkos::Threads>("Threads-float", beg, end, runs);
  explicit_dynamics::driver<double,Kokkos::Threads>("Threads-double", beg, end, runs);

  Kokkos::Threads::finalize();
}//test_host

}// namespace


