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

#include "Kokkos_Core.hpp"

#ifdef KOKKOS_HAVE_PTHREAD
#include "Stokhos_Threads_CrsProductTensor.hpp"
#endif
#ifdef KOKKOS_HAVE_OPENMP
#include "Stokhos_OpenMP_CrsProductTensor.hpp"
#endif

#include "TestStochastic.hpp"

namespace unit_test {

template<typename Scalar, typename Device>
struct performance_test_driver {

  static void run(bool test_flat, bool test_orig, bool test_deg, bool test_lin,
                  bool test_block, bool symmetric, bool mkl) {

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
#ifdef __MIC__
      nGrid = 32 ;
#else
      nGrid = 32 ;
#endif
      nIter = 1 ;
      // Something funny happens when we go to larger problem sizes on the
      // MIC where it appears to slow down subsequent calculations (i.e.,
      // the degree 5 cases will run slower).  Maybe it is getting too hot?
#ifdef __MIC__
      performance_test_driver_poly<Scalar,Device,Stokhos::DefaultMultiply>(
        3 , 1 , 9 , nGrid , nIter , test_block , symmetric );
#else
      performance_test_driver_poly<Scalar,Device,Stokhos::DefaultMultiply>(
        3 , 1 , 12 , nGrid , nIter , test_block , symmetric );
#endif
      performance_test_driver_poly<Scalar,Device,Stokhos::DefaultMultiply>(
        5 , 1,  6 , nGrid , nIter , test_block , symmetric );
    }

    // Just polynomial methods compared against original
    if (test_deg) {
 #ifdef __MIC__
      nGrid = 32 ;
#else
      nGrid = 64 ;
#endif
      nIter = 1 ;
      performance_test_driver_poly_deg<Scalar,Device,Stokhos::DefaultMultiply>(
        3 , 1 , 12 , nGrid , nIter , test_block , symmetric );
    }

    // Just polynomial methods compared against original
    if (test_lin) {
#ifdef __MIC__
      nGrid = 32 ;
#else
      nGrid = 64 ;
#endif
      nIter = 10 ;
      performance_test_driver_linear<Scalar,Device,Stokhos::DefaultMultiply>(
        31 ,  255 , 32 , nGrid , nIter , test_block , symmetric );
    }

    //------------------------------
  }

};

}

template <typename Scalar, typename Device>
int mainHost(bool test_flat, bool test_orig, bool test_deg, bool test_lin,
             bool test_block, bool symmetric, bool mkl)
{
  const size_t team_count =
    Kokkos::hwloc::get_available_numa_count() *
    Kokkos::hwloc::get_available_cores_per_numa();
  const size_t threads_per_team =
    Kokkos::hwloc::get_available_threads_per_core();

  Device::initialize( team_count * threads_per_team );

  std::string name = "Host";
#ifdef KOKKOS_HAVE_PTHREAD
  Kokkos::Threads::print_configuration( std::cout );
  if (Kokkos::Impl::is_same<Device,Kokkos::Threads>::value)
    name = "Threads";
#endif
#ifdef KOKKOS_HAVE_OPENMP
  Kokkos::OpenMP::print_configuration( std::cout );
  if (Kokkos::Impl::is_same<Device,Kokkos::OpenMP>::value)
    name = "OpenMP";
#endif
  std::cout << std::endl << "\"" << name << " Performance with "
            << team_count * threads_per_team << " threads\"" << std::endl ;

  unit_test::performance_test_driver<Scalar,Device>::run(
    test_flat, test_orig, test_deg, test_lin, test_block, symmetric, mkl);

  Device::finalize();

  return 0 ;
}

#ifdef KOKKOS_HAVE_PTHREAD
template int mainHost<float,Kokkos::Threads>(bool, bool, bool, bool, bool, bool, bool);
template int mainHost<double,Kokkos::Threads>(bool, bool, bool, bool, bool, bool, bool);
#endif

#ifdef KOKKOS_HAVE_OPENMP
template int mainHost<float,Kokkos::OpenMP>(bool, bool, bool, bool, bool, bool, bool);
template int mainHost<double,Kokkos::OpenMP>(bool, bool, bool, bool, bool, bool, bool);
#endif
