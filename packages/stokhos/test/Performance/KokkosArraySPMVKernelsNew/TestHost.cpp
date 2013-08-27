// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include <iostream>

#include "Kokkos_hwloc.hpp"

#include "Stokhos_ConfigDefs.h"
#if defined(HAVE_STOKHOS_OPENMP) && defined(HAVE_STOKHOS_MKL)
#include <omp.h>
#endif

#include "TestStochastic.hpp"

#include "Stokhos_Threads_CrsMatrix.hpp"
#include "Stokhos_Threads_BlockCrsMatrix.hpp"
#include "Stokhos_Threads_StochasticProductTensor.hpp"
#include "Stokhos_Threads_SymmetricDiagonalSpec.hpp"
#include "Stokhos_Threads_CrsProductTensor.hpp"
#include "Stokhos_Threads_TiledCrsProductTensor.hpp"
#include "Stokhos_Threads_CooProductTensor.hpp"
#include "Stokhos_Threads_FlatSparse3Tensor.hpp"
#include "Stokhos_Threads_FlatSparse3Tensor_kji.hpp"
#include "Stokhos_Threads_LexicographicBlockSparse3Tensor.hpp"
#include "Stokhos_Threads_LinearSparse3Tensor.hpp"

namespace unit_test {

template<typename Scalar>
struct performance_test_driver<Scalar,Kokkos::Threads> {

  static void run(bool test_flat, bool test_orig, bool test_deg, bool test_lin,
                  bool test_block, bool symmetric, bool mkl) {
    typedef Kokkos::Threads Device;

    int nGrid;
    int nIter;

    // All methods compared against flat-original
    if (test_flat) {
      nGrid = 12 ;
      nIter = 1 ;
      performance_test_driver_all<Scalar,Device>(
        3 , 1 ,  9 , nGrid , nIter , test_block , symmetric );
      performance_test_driver_all<Scalar,Device>(
        5 , 1 ,  5 , nGrid , nIter , test_block , symmetric );
    }

    // Just polynomial methods compared against original
    if (test_orig) {
      nGrid = 64 ;
      nIter = 1 ;
      if (mkl) {
#ifdef HAVE_STOKHOS_MKL
        performance_test_driver_poly<Scalar,Device,Stokhos::MKLSparseMatOps>(
          3 , 1 , 12 , nGrid , nIter , test_block , symmetric );
        performance_test_driver_poly<Scalar,Device,Stokhos::MKLSparseMatOps>(
          5 , 1 ,  6 , nGrid , nIter , test_block , symmetric );
#else
        std::cout << "MKL support not enabled!" << std::endl;
#endif
      }
      else {
        performance_test_driver_poly<Scalar,Device,Stokhos::DefaultSparseMatOps>(
          3 , 1 , 12 , nGrid , nIter , test_block , symmetric );
        performance_test_driver_poly<Scalar,Device,Stokhos::DefaultSparseMatOps>(
          5 , 1 ,  6 , nGrid , nIter , test_block , symmetric );
      }
    }

    // Just polynomial methods compared against original
    if (test_deg) {
      nGrid = 64 ;
      nIter = 1 ;
      if (mkl) {
#ifdef HAVE_STOKHOS_MKL
        performance_test_driver_poly_deg<Scalar,Device,Stokhos::MKLSparseMatOps>(
          3 , 1 , 12 , nGrid , nIter , test_block , symmetric );
#else
        std::cout << "MKL support not enabled!" << std::endl;
#endif
      }
      else {
        performance_test_driver_poly_deg<Scalar,Device,Stokhos::DefaultSparseMatOps>(
          3 , 1 , 12 , nGrid , nIter , test_block , symmetric );
      }
    }

    // Just polynomial methods compared against original
    if (test_lin) {
      nGrid = 64 ;
      nIter = 10 ;
      performance_test_driver_linear<Scalar,Device,Stokhos::DefaultSparseMatOps>(
        31 ,  255 , 32 , nGrid , nIter , test_block , symmetric );
    }

    //------------------------------
  }

};

}

template <typename Scalar>
int mainHost(bool test_flat, bool test_orig, bool test_deg, bool test_lin,
             bool test_block, bool symmetric, bool mkl)
{
  const std::pair<unsigned,unsigned> core_topo =
    Kokkos::hwloc::get_core_topology();
  const size_t core_capacity = Kokkos::hwloc::get_core_capacity();
  //const size_t core_capacity = 1;

  const size_t gang_count = core_topo.first ;
  const size_t gang_worker_count = core_topo.second * core_capacity;

#if defined(HAVE_STOKHOS_OPENMP) && defined(HAVE_STOKHOS_MKL)
  // Call a little OpenMP parallel region so that MKL will get the right
  // number of threads.  This isn't perfect in that the thread binding
  // doesn't seem right, and only works at all when using GNU threads with MKL.
#pragma omp parallel
  {
    int numThreads = omp_get_num_threads();
#pragma omp single
    std::cout << " num_omp_threads = " << numThreads << std::endl;
  }
#endif

  Kokkos::Threads::initialize( std::make_pair(gang_count , gang_worker_count),
                               core_topo );
  Kokkos::Threads::print_configuration( std::cout );

  std::cout << std::endl << "\"Host Performance with "
            << gang_count * gang_worker_count << " threads\"" << std::endl ;

  unit_test::performance_test_driver<Scalar,Kokkos::Threads>::run(
    test_flat, test_orig, test_deg, test_lin, test_block, symmetric, mkl);

  Kokkos::Threads::finalize();

  return 0 ;
}

template int mainHost<float>(bool, bool, bool, bool, bool, bool, bool);
template int mainHost<double>(bool, bool, bool, bool, bool, bool, bool);
