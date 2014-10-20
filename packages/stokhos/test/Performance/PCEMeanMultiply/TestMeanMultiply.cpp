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

// Devices
#include "Kokkos_Core.hpp"

// Utilities
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#ifdef KOKKOS_HAVE_CUDA
#include "cuda_runtime_api.h"
#endif

template <typename Scalar, typename Ordinal, typename Device>
void performance_test_driver( const Ordinal nGrid,
                              const Ordinal nIter,
                              const Ordinal order,
                              const Ordinal min_var,
                              const Ordinal max_var );

int main(int argc, char *argv[])
{
  bool success = true;
  bool verbose = false;
  try {

    const size_t num_sockets = Kokkos::hwloc::get_available_numa_count();
    const size_t num_cores_per_socket =
      Kokkos::hwloc::get_available_cores_per_numa();
    // const size_t num_threads_per_core =
    //   Kokkos::hwloc::get_available_threads_per_core();
    // const size_t num_threads =
    //   num_sockets * num_cores_per_socket * num_threads_per_core;

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This test performance of mean-based UQ::PCE multiply routines.\n");
    int nGrid = 32;
    CLP.setOption("n", &nGrid, "Number of mesh points in the each direction");
    int nIter = 10;
    CLP.setOption("ni", &nIter, "Number of multiply iterations");
    int order = 3;
    CLP.setOption("order", &order, "Polynomial order");
    int dim_min = 1;
    CLP.setOption("dmin", &dim_min, "Starting stochastic dimension");
    int dim_max = 12;
    CLP.setOption("dmax", &dim_max, "Stopping stochastic dimension");
    int numa = num_sockets;
    CLP.setOption("numa", &numa,  "Number of numa nodes");
    int cores = num_cores_per_socket;
    CLP.setOption("cores", &cores, "Cores per numa node");
#ifdef KOKKOS_HAVE_PTHREAD
    int threads = 0;
    CLP.setOption("threads", &threads, "Number of threads for Threads device");
#endif
#ifdef KOKKOS_HAVE_OPENMP
    int openmp = 0;
    CLP.setOption("openmp", &openmp, "Number of threads for OpenMP device");
#endif
#ifdef KOKKOS_HAVE_CUDA
    bool cuda = false;
    CLP.setOption("cuda", "no-cuda", &cuda, "Enable Cuda device");
    int device_id = 0;
    CLP.setOption("device", &device_id, "CUDA device ID");
#endif
    CLP.parse( argc, argv );

    typedef int Ordinal;
    typedef double Scalar;

#ifdef KOKKOS_HAVE_PTHREAD
    if (threads > 0) {
      typedef Kokkos::Threads Device;

      Kokkos::Threads::initialize(threads, numa, cores);

      std::cout << std::endl
                << "Threads performance with " << threads
                << " threads, " << numa << " numas, " << cores
                << " cores/numa:" << std::endl;

      performance_test_driver<Scalar,Ordinal,Device>(
        nGrid, nIter, order, dim_min, dim_max);

      Kokkos::Threads::finalize();
    }
#endif

#ifdef KOKKOS_HAVE_OPENMP
    if (openmp > 0) {
      typedef Kokkos::OpenMP Device;

      Kokkos::OpenMP::initialize(openmp, numa, cores);

      std::cout << std::endl
                << "OpenMP performance with " << openmp
                << " threads, " << numa << " numas, " << cores
                << " cores/numa:" << std::endl;

      performance_test_driver<Scalar,Ordinal,Device>(
        nGrid, nIter, order, dim_min, dim_max);

      Kokkos::OpenMP::finalize();
    }
#endif

#ifdef KOKKOS_HAVE_CUDA
    if (cuda) {
      typedef Kokkos::Cuda Device;

      Kokkos::HostSpace::execution_space::initialize();
      Kokkos::Cuda::initialize(Kokkos::Cuda::SelectDevice(device_id));

      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, device_id);
      std::cout << std::endl
                << "CUDA performance for device " << device_id << " ("
                << deviceProp.name << "):"
                << std::endl;

      performance_test_driver<Scalar,Ordinal,Device>(
        nGrid, nIter, order, dim_min, dim_max);

      Kokkos::HostSpace::execution_space::finalize();
      Kokkos::Cuda::finalize();
    }
#endif

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  if (success)
    return 0;
  return -1;
}
