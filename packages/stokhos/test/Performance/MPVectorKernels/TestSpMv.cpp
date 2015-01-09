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

// Tests
#include "TestSpMv.hpp"

// Devices
#include "Kokkos_Core.hpp"

// Utilities
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#ifdef KOKKOS_HAVE_CUDA
#include "cuda_runtime_api.h"
#endif

template <typename Storage>
void mainHost(int nGrid, int nIter, Kokkos::DeviceConfig dev_config);
template <typename Storage>
void mainCuda(int nGrid, int nIter, Kokkos::DeviceConfig dev_config);

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
#ifdef KOKKOS_HAVE_PTHREAD
    bool threads = true;
    CLP.setOption("threads", "no-threads", &threads, "Enable Threads device");
#endif
#ifdef KOKKOS_HAVE_OPENMP
    bool openmp = true;
    CLP.setOption("openmp", "no-openmp", &openmp, "Enable OpenMP device");
#endif
#ifdef KOKKOS_HAVE_CUDA
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

#ifdef KOKKOS_HAVE_PTHREAD
    if (threads) {
      typedef Kokkos::Threads Device;
      typedef Stokhos::StaticFixedStorage<Ordinal,Scalar,1,Device> Storage;

      Kokkos::Threads::initialize(num_cores*num_hyper_threads);

      std::cout << std::endl
                << "Threads performance with " << num_cores*num_hyper_threads
                << " threads:" << std::endl;

      Kokkos::DeviceConfig dev_config(num_cores,
                                       threads_per_vector,
                                       num_hyper_threads / threads_per_vector);

      mainHost<Storage>(nGrid, nIter, dev_config);

      Kokkos::Threads::finalize();
    }
#endif

#ifdef KOKKOS_HAVE_OPENMP
    if (openmp) {
      typedef Kokkos::OpenMP Device;
      typedef Stokhos::StaticFixedStorage<Ordinal,Scalar,1,Device> Storage;

      Kokkos::OpenMP::initialize(num_cores*num_hyper_threads);

      std::cout << std::endl
                << "OpenMP performance with " << num_cores*num_hyper_threads
                << " threads:" << std::endl;

      Kokkos::DeviceConfig dev_config(num_cores,
                                       threads_per_vector,
                                       num_hyper_threads / threads_per_vector);

      mainHost<Storage>(nGrid, nIter, dev_config);

      Kokkos::OpenMP::finalize();
    }
#endif

#ifdef KOKKOS_HAVE_CUDA
    if (cuda) {
      typedef Kokkos::Cuda Device;
      typedef Stokhos::StaticFixedStorage<Ordinal,Scalar,1,Device> Storage;

      Kokkos::HostSpace::execution_space::initialize();
      Kokkos::Cuda::initialize(Kokkos::Cuda::SelectDevice(device_id));

      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, device_id);
      std::cout << std::endl
                << "CUDA performance for device " << device_id << " ("
                << deviceProp.name << "):"
                << std::endl;

      Kokkos::DeviceConfig dev_config(
        num_cuda_blocks,
        cuda_threads_per_vector,
        cuda_threads_per_vector == 0 ? 0 : cuda_block_size / cuda_threads_per_vector);

      mainCuda<Storage>(nGrid,nIter,dev_config);

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
