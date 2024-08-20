// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

// Tests
#include "TestSpMv.hpp"

// Devices
#include "Kokkos_Core.hpp"

// Utilities
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#ifdef KOKKOS_ENABLE_CUDA
#include "cuda_runtime_api.h"
#endif

template <typename Storage>
void mainHost(int nGrid, int nIter, KokkosSparse::DeviceConfig dev_config);
template <typename Storage>
void mainCuda(int nGrid, int nIter, KokkosSparse::DeviceConfig dev_config);

int main(int argc, char *argv[])
{
  bool success = true;
  bool verbose = false;
  try {

    const size_t num_sockets = Kokkos::hwloc::get_available_numa_count();
    const size_t num_cores_per_socket =
      Kokkos::hwloc::get_available_cores_per_numa();
    const size_t num_threads_per_core =
      Kokkos::hwloc::get_available_threads_per_core();

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This test performance of MP::Vector multiply routines.\n");
    int nGrid = 32;
    CLP.setOption("n", &nGrid, "Number of mesh points in the each direction");
    int nIter = 10;
    CLP.setOption("ni", &nIter, "Number of multiply iterations");
    int num_cores = num_cores_per_socket * num_sockets;
    CLP.setOption("cores", &num_cores,
                  "Number of CPU cores to use (defaults to all)");
    int num_hyper_threads = num_threads_per_core;
    CLP.setOption("hyperthreads", &num_hyper_threads,
                  "Number of hyper threads per core to use (defaults to all)");
    int threads_per_vector = 1;
    CLP.setOption("threads_per_vector", &threads_per_vector,
                  "Number of threads to use within each vector");
#ifdef KOKKOS_ENABLE_THREADS
    bool threads = true;
    CLP.setOption("threads", "no-threads", &threads, "Enable Threads device");
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    bool openmp = true;
    CLP.setOption("openmp", "no-openmp", &openmp, "Enable OpenMP device");
#endif
#ifdef KOKKOS_ENABLE_CUDA
    bool cuda = true;
    CLP.setOption("cuda", "no-cuda", &cuda, "Enable Cuda device");
    int cuda_threads_per_vector = 16;
    CLP.setOption("cuda_threads_per_vector", &cuda_threads_per_vector,
                  "Number of Cuda threads to use within each vector");
    int cuda_block_size = 0;
    CLP.setOption("cuda_block_size", &cuda_block_size,
                  "Cuda block size (0 implies the default choice)");
    int num_cuda_blocks = 0;
    CLP.setOption("num_cuda_blocks", &num_cuda_blocks,
                  "Number of Cuda blocks (0 implies the default choice)");
    int device_id = 0;
    CLP.setOption("device", &device_id, "CUDA device ID");
#endif
    CLP.parse( argc, argv );

    typedef int Ordinal;
    typedef double Scalar;

#ifdef KOKKOS_ENABLE_THREADS
    if (threads) {
      typedef Kokkos::Threads Device;
      typedef Stokhos::StaticFixedStorage<Ordinal,Scalar,1,Device> Storage;

      Kokkos::InitializationSettings init_args;
      init_args.set_num_threads(num_cores*num_hyper_threads);
      Kokkos::initialize( init_args );

      std::cout << std::endl
                << "Threads performance with " << num_cores*num_hyper_threads
                << " threads:" << std::endl;

      KokkosSparse::DeviceConfig dev_config(num_cores,
                                       threads_per_vector,
                                       num_hyper_threads / threads_per_vector);

      mainHost<Storage>(nGrid, nIter, dev_config);

      Kokkos::finalize();
    }
#endif

#ifdef KOKKOS_ENABLE_OPENMP
    if (openmp) {
      typedef Kokkos::OpenMP Device;
      typedef Stokhos::StaticFixedStorage<Ordinal,Scalar,1,Device> Storage;

      Kokkos::InitializationSettings init_args;
      init_args.set_num_threads(num_cores*num_hyper_threads);
      Kokkos::initialize( init_args );

      std::cout << std::endl
                << "OpenMP performance with " << num_cores*num_hyper_threads
                << " threads:" << std::endl;

      KokkosSparse::DeviceConfig dev_config(num_cores,
                                       threads_per_vector,
                                       num_hyper_threads / threads_per_vector);

      mainHost<Storage>(nGrid, nIter, dev_config);

      Kokkos::finalize();
    }
#endif

#ifdef KOKKOS_ENABLE_CUDA
    if (cuda) {
      typedef Kokkos::Cuda Device;
      typedef Stokhos::StaticFixedStorage<Ordinal,Scalar,1,Device> Storage;

      Kokkos::InitializationSettings init_args;
      init_args.set_device_id(device_id);
      Kokkos::initialize( init_args );

      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, device_id);
      std::cout << std::endl
                << "CUDA performance for device " << device_id << " ("
                << deviceProp.name << "):"
                << std::endl;

      KokkosSparse::DeviceConfig dev_config(
        num_cuda_blocks,
        cuda_threads_per_vector,
        cuda_threads_per_vector == 0 ? 0 : cuda_block_size / cuda_threads_per_vector);

      mainCuda<Storage>(nGrid,nIter,dev_config);

      Kokkos::finalize();
    }
#endif

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  if (success)
    return 0;
  return -1;
}
