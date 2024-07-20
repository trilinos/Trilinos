// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

// Devices
#include "Kokkos_Core.hpp"

// Utilities
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#ifdef KOKKOS_ENABLE_CUDA
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
#ifdef KOKKOS_ENABLE_THREADS
    int threads = 0;
    CLP.setOption("threads", &threads, "Number of threads for Threads device");
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    int openmp = 0;
    CLP.setOption("openmp", &openmp, "Number of threads for OpenMP device");
#endif
#ifdef KOKKOS_ENABLE_CUDA
    bool cuda = false;
    CLP.setOption("cuda", "no-cuda", &cuda, "Enable Cuda device");
    int device_id = 0;
    CLP.setOption("device", &device_id, "CUDA device ID");
#endif
    CLP.parse( argc, argv );

    typedef int Ordinal;
    typedef double Scalar;

#ifdef KOKKOS_ENABLE_THREADS
    if (threads > 0) {
      typedef Kokkos::Threads Device;

      Kokkos::InitializationSettings init_args;
      init_args.set_num_threads(threads);
      Kokkos::initialize( init_args );

      std::cout << std::endl
                << "Threads performance with " << threads
                << " threads, " << numa << " numas, " << cores
                << " cores/numa:" << std::endl;

      performance_test_driver<Scalar,Ordinal,Device>(
        nGrid, nIter, order, dim_min, dim_max);

      Kokkos::finalize();
    }
#endif

#ifdef KOKKOS_ENABLE_OPENMP
    if (openmp > 0) {
      typedef Kokkos::OpenMP Device;

      Kokkos::InitializationSettings init_args;
      init_args.set_num_threads(openmp);
      Kokkos::initialize( init_args );

      std::cout << std::endl
                << "OpenMP performance with " << openmp
                << " threads, " << numa << " numas, " << cores
                << " cores/numa:" << std::endl;

      performance_test_driver<Scalar,Ordinal,Device>(
        nGrid, nIter, order, dim_min, dim_max);

      Kokkos::finalize();
    }
#endif

#ifdef KOKKOS_ENABLE_CUDA
    if (cuda) {
      typedef Kokkos::Cuda Device;

      Kokkos::InitializationSettings init_args;
      init_args.set_device_id(device_id);
      Kokkos::initialize( init_args );

      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, device_id);
      std::cout << std::endl
                << "CUDA performance for device " << device_id << " ("
                << deviceProp.name << "):"
                << std::endl;

      performance_test_driver<Scalar,Ordinal,Device>(
        nGrid, nIter, order, dim_min, dim_max);

      Kokkos::finalize();
    }
#endif

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  if (success)
    return 0;
  return -1;
}
