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

#include "KokkosArray_hwloc.hpp"

#include <TestStochastic.hpp>

#include <Host/KokkosArray_Host_ProductTensor.hpp>
#include <Host/KokkosArray_Host_StochasticProductTensor.hpp>
#include <Host/KokkosArray_Host_SymmetricDiagonalSpec.hpp>
#include <Host/KokkosArray_Host_BlockCrsMatrix.hpp>
#include <Host/KokkosArray_Host_CrsMatrix.hpp>

namespace unit_test {

template<typename Scalar>
struct performance_test_driver<Scalar,KokkosArray::Host> {

  static void run(bool test_flat, bool test_orig, bool test_block) {
    typedef KokkosArray::Host Device;

    int nGrid;
    int nIter; 
    
    // All methods compared against flat-original
    if (test_flat) {
      nGrid = 12 ;
      nIter = 1 ; 
      performance_test_driver_all<Scalar,Device>( 
	3 , 1 ,  9 , nGrid , nIter , test_block );
      performance_test_driver_all<Scalar,Device>( 
	5 , 1 ,  5 , nGrid , nIter , test_block );
    }
    
    // Just polynomial methods compared against original
    if (test_orig) {
      nGrid = 64 ;
      nIter = 1 ; 
      performance_test_driver_poly<Scalar,Device>( 
	3 , 1 , 12 , nGrid , nIter , test_block );
      performance_test_driver_poly<Scalar,Device>( 
	5 , 1 ,  6 , nGrid , nIter , test_block );
    }
    
    //------------------------------
  }

};

}

template <typename Scalar>
int mainHost(bool test_flat, bool test_orig, bool test_block)
{
  const std::pair<unsigned,unsigned> core_topo =
    KokkosArray::hwloc::get_core_topology();
  const size_t core_capacity = KokkosArray::hwloc::get_core_capacity();

  const size_t gang_count = core_topo.first ;
  const size_t gang_worker_count = core_topo.second * core_capacity;

  KokkosArray::Host::initialize( gang_count , gang_worker_count );

  std::cout << std::endl << "\"Host Performance with "
            << gang_count * gang_worker_count << " threads\"" << std::endl ;
  unit_test::performance_test_driver<Scalar,KokkosArray::Host>::run(
    test_flat, test_orig, test_block);

  KokkosArray::Host::finalize();

  return 0 ;
}

template int mainHost<float>(bool, bool, bool);
template int mainHost<double>(bool, bool, bool);
