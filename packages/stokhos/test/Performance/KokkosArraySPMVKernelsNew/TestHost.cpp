// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include <iostream>

#include "Stokhos_ConfigDefs.h"
#if defined(HAVE_STOKHOS_OPENMP) && defined(HAVE_STOKHOS_MKL)
#include <omp.h>
#endif

#include "TestStochastic.hpp"

#include "Stokhos_Host_CrsMatrix.hpp"
#include "Stokhos_Host_BlockCrsMatrix.hpp"
#include "Stokhos_Host_StochasticProductTensor.hpp"
#include "Stokhos_Host_CrsProductTensor.hpp"
#include "Stokhos_Host_FlatSparse3Tensor.hpp"
#include "Stokhos_Host_FlatSparse3Tensor_kji.hpp"
#include "Stokhos_Host_LexicographicBlockSparse3Tensor.hpp"

namespace unit_test {

template<typename Scalar>
struct performance_test_driver<Scalar,KokkosArray::Host> {

  static void run(bool test_flat, bool test_orig, bool test_block,
                  bool symmetric, bool mkl) {
    typedef KokkosArray::Host Device;

    int nGrid;
    int nIter;

    // All methods compared against flat-original
    // if (test_flat) {
    //   nGrid = 12 ;
    //   nIter = 1 ;
    //   performance_test_driver_all<Scalar,Device>(
    //  3 , 1 ,  9 , nGrid , nIter , test_block );
    //   performance_test_driver_all<Scalar,Device>(
    //  5 , 1 ,  5 , nGrid , nIter , test_block );
    // }

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
          5 , 4 ,  4 , nGrid , nIter , test_block , symmetric );
      }
    }

    //------------------------------
  }

};

}

template <typename Scalar>
int mainHost(bool test_flat, bool test_orig, bool test_block, bool symmetric,
             bool mkl)
{
  const size_t gang_count = KokkosArray::Host::detect_gang_capacity();
  const size_t gang_worker_count = KokkosArray::Host::detect_gang_worker_capacity() ;

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

  KokkosArray::Host::initialize( gang_count , gang_worker_count );

  std::cout << std::endl << "\"Host Performance with "
            << gang_count * gang_worker_count << " threads\"" << std::endl ;

  unit_test::performance_test_driver<Scalar,KokkosArray::Host>::run(
    test_flat, test_orig, test_block, symmetric, mkl);

  KokkosArray::Host::finalize();

  return 0 ;
}

template int mainHost<float>(bool, bool, bool, bool, bool);
template int mainHost<double>(bool, bool, bool, bool, bool);
