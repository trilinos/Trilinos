//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#include "Kokkos_ThrustGPUNode.hpp"

// mfh 08 Jan 2015: Don't enable the contents of this file unless the
// appropriate CMake option is enabled.  This avoids deprecation
// warnings once we deprecate this Node type.
#ifdef HAVE_TPETRACLASSIC_THRUST

#include <Teuchos_Assert.hpp>
#include <iostream>
#include <cuda_runtime.h>

#ifdef HAVE_KOKKOSCLASSIC_KOKKOSCORE
#include "Kokkos_Core.hpp"
#endif

namespace KokkosClassic {

  ThrustGPUNode::ThrustGPUNode () {
    using std::cout;
    using std::cerr;
    using std::endl;

    Teuchos::ParameterList params = getDefaultParameters();
    int device = params.get<int>("Device Number");
    const int verbose = params.get<int>("Verbose");

    // set device
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    TEUCHOS_TEST_FOR_EXCEPTION(
      deviceCount == 0, std::runtime_error,
      "KokkosClassic::ThrustGPUNode constructor: system has no CUDA devices.");

    if (device < 0 || device >= deviceCount) {
      cerr << "KokkosClassic::ThrustGPUNode constructor: The specified device "
           << "number " << device << " is not valid.  Using Device 0 instead."
           << endl;
      device = 0;
    }
    cudaDeviceProp deviceProp;
    cudaSetDevice (device);
    cudaGetDeviceProperties (&deviceProp, device);
    // as of CUDA 2.1, device prop contains the following fields
    // char name[256];
    // size_t totalGlobalMem, sharedMemPerBlock;
    // int regsPerBlock, warpSize;
    // size_t memPitch;
    // int maxThreadsPerBlock, maxThreadsDim[3], maxGridSize[3];
    // size_t totalConstMem;
    // int major, minor;
    // int clockRate;
    // size_t textureAlignment;
    // int deviceOverlap;
    // int multiProcessorCount;
    // int kernelExecTimeoutEnabled;
    if (verbose) {
      cout << "KokkosClassic::ThrustGPUNode attached to device #" << device << " \""
           << deviceProp.name << "\", of compute capability "
           << deviceProp.major << "." << deviceProp.minor << endl;
    }
    totalMem_ = deviceProp.totalGlobalMem;

#if defined(HAVE_KOKKOSCLASSIC_KOKKOSCORE) && defined(KOKKOS_HAVE_CUDA)
    if (! Kokkos::HostSpace::execution_space::is_initialized ())
      Kokkos::HostSpace::execution_space::initialize ();
    if (! Kokkos::Cuda::is_initialized ())
      Kokkos::Cuda::initialize (Kokkos::Cuda::SelectDevice (device));
#endif
  }

  ThrustGPUNode::ThrustGPUNode (Teuchos::ParameterList &pl)
  {
    using std::cout;
    using std::cerr;
    using std::endl;

    // get node parameters
    Teuchos::ParameterList params = getDefaultParameters ();
    params.setParameters (pl);
    int device = params.get<int>("Device Number");
    const int verbose = params.get<int>("Verbose");

    // set device
    int deviceCount; cudaGetDeviceCount(&deviceCount);
    TEUCHOS_TEST_FOR_EXCEPTION(
      deviceCount == 0, std::runtime_error,
      "Kokkos::ThrustGPUNode constructor: system has no CUDA devices."
    );
    if (device < 0 || device >= deviceCount) {
      cerr << "Kokkos::ThrustGPUNode constructor: The specified device number "
           << device << " is not valid. Using Device 0." << endl;
      device = 0;
    }
    cudaDeviceProp deviceProp;
    cudaSetDevice (device);
    cudaGetDeviceProperties (&deviceProp, device);
    // as of CUDA 2.1, device prop contains the following fields
    // char name[256];
    // size_t totalGlobalMem, sharedMemPerBlock;
    // int regsPerBlock, warpSize;
    // size_t memPitch;
    // int maxThreadsPerBlock, maxThreadsDim[3], maxGridSize[3];
    // size_t totalConstMem;
    // int major, minor;
    // int clockRate;
    // size_t textureAlignment;
    // int deviceOverlap;
    // int multiProcessorCount;
    // int kernelExecTimeoutEnabled;
    if (verbose) {
      cout << "Kokkos::ThrustGPUNode attached to device #" << device << " \""
           << deviceProp.name << "\", of compute capability "
           << deviceProp.major << "." << deviceProp.minor << endl;
    }
    totalMem_ = deviceProp.totalGlobalMem;

#if defined(HAVE_KOKKOSCLASSIC_KOKKOSCORE) && defined(KOKKOS_HAVE_CUDA)
    if (! Kokkos::HostSpace::execution_space::is_initialized ()) {
      Kokkos::HostSpace::execution_space::initialize ();
    }
    if (! Kokkos::Cuda::is_initialized ()) {
      Kokkos::Cuda::initialize (Kokkos::Cuda::SelectDevice (device));
    }
#endif
  }

  ThrustGPUNode::~ThrustGPUNode() {}

  Teuchos::ParameterList ThrustGPUNode::getDefaultParameters () {
    Teuchos::ParameterList params;
    params.set ("Verbose",       0);
    params.set ("Device Number", 0);
    return params;
  }

  void ThrustGPUNode::sync () const {
    cudaError err = cudaThreadSynchronize();
    TEUCHOS_TEST_FOR_EXCEPTION(
      cudaSuccess != err, std::runtime_error,
      "KokkosClassic::ThrustGPUNode::sync: cudaThreadSynchronize() returned "
      "error " << err);
  }

  std::string ThrustGPUNode::name () {
    return "ThrustGPU";
  }

} // namespace KokkosClassic

#endif // HAVE_TPETRACLASSIC_THRUST
