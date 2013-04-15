/*
//@HEADER
// ************************************************************************
//
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
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

#include <KokkosArray_Value.hpp>
#include <KokkosArray_MultiVector.hpp>
#include <KokkosArray_MDArray.hpp>

#include <KokkosArray_Host.hpp>

#include <impl/KokkosArray_Host_macros.hpp>
#include <explicit_dynamics_app.hpp>
#include <impl/KokkosArray_Clear_macros.hpp>

namespace Test{

void test_Host( int beg, int end, int runs, int threads)
{
  const std::pair<unsigned,unsigned> core_topo = KokkosArray::hwloc::get_core_topology();
  const unsigned                     core_cap  = KokkosArray::hwloc::get_core_capacity();

  if ( 0 < threads ) {
    const size_t node_thread_count = ( threads + core_topo.first - 1 ) / core_topo.first ;

    KokkosArray::Host::initialize( core_topo.first , node_thread_count );

    std::cout << std::endl << "\"Host with manually set threads = \" , "
              << core_topo.first * node_thread_count << std::endl ;
  }
  else {
    KokkosArray::Host::initialize( core_topo.first , core_topo.second * core_cap );

    std::cout << std::endl << "\"Host with detected sequential threads = \" , "
              << core_topo.first * node_thread_count << std::endl ;
  }

  explicit_dynamics::driver<float,KokkosArray::Host>("Host-float", beg, end, runs);
  explicit_dynamics::driver<double,KokkosArray::Host>("Host-double", beg, end, runs);

  KokkosArray::Host::finalize();
}//test_host

}// namespace


