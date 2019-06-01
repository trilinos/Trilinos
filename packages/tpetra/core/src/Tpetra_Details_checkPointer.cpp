// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
// @HEADER

#include "Tpetra_Details_checkPointer.hpp"
#include "Kokkos_Core.hpp"
#ifdef HAVE_TPETRACORE_CUDA
# include "cuda_runtime_api.h"
#endif // HAVE_TPETRACORE_CUDA

namespace Tpetra {
namespace Details {
namespace Impl {

enum class EMemoryType {
  HOST,
  CUDA_HOST_PINNED,
  CUDA,
  CUDA_UVM,
  ERROR
};

EMemoryType getCudaMemoryType (const void* ptr)
{
#ifdef HAVE_TPETRACORE_CUDA
  cudaPointerAttributes attr;
  const cudaError_t err = cudaPointerGetAttributes (&attr, ptr);
  (void) cudaGetLastError (); // reset error state for future CUDA ops

  if (err == cudaErrorInvalidValue) {
    // CUDA 11.0 supports passing in an unregistered host pointer.  In
    // that case, attr.type will be cudaMemoryTypeUnregistered.  CUDA
    // 9.2 doesn't yet have the 'type' field in the
    // cudaPointerAttributes struct, and this function just returns
    // this error code if given an unregistered host pointer.
    return EMemoryType::HOST;
  }
  else if (err != cudaSuccess) {
    return EMemoryType::ERROR;
  }

# if 1
  // cudaPointerAttributes::type doesn't exist yet in CUDA 9.2.  We
  // must use memoryType, which does not distinguish between types of
  // device memory.
  if (attr.memoryType == cudaMemoryTypeDevice) {
    if (attr.isManaged) {
      return EMemoryType::CUDA_UVM;
    }
    else { // if (attr.hostPointer == nullptr) {
      return EMemoryType::CUDA;
    }
    // else {
    //   return EMemoryType::ERROR;
    // }
  }
  else if (attr.memoryType == cudaMemoryTypeHost) {
    if (attr.devicePointer == nullptr) {
      return EMemoryType::HOST; // not device accessible
    }
    else { // device-accessible host memory is pinned
      return EMemoryType::CUDA_HOST_PINNED;
    }
  }
  else {
    return EMemoryType::ERROR;
  }

# else

  const enum cudaMemoryType theType = attr.type;
  if (theType == cudaMemoryTypeManaged) {
    return EMemoryType::CUDA_UVM;
  }
  else if (theType == cudaMemoryTypeDevice) {
    return EMemoryType::CUDA;
  }
  else if (theType == cudaMemoryTypeHost) {
    return EMemoryType::CUDA_HOST_PINNED;
  }
  else {
    return EMemoryType::HOST;
  }

# endif // 1
#else // NOT HAVE_TPETRACORE_CUDA
  return EMemoryType::HOST;
#endif // HAVE_TPETRACORE_CUDA
}

#ifdef HAVE_TPETRACORE_CUDA
bool
CheckPointerAccessibility<Kokkos::Cuda>::
accessible (const void* ptr,
            const Kokkos::Cuda& /* space */)
{
  const EMemoryType type = getCudaMemoryType (ptr);
  return type != EMemoryType::HOST && type != EMemoryType::ERROR;
}
#endif // HAVE_TPETRACORE_CUDA

bool isHostAccessible (const void* ptr)
{
#ifdef HAVE_TPETRACORE_CUDA
  return getCudaMemoryType (ptr) != EMemoryType::CUDA;
#else
  return true;
#endif // HAVE_TPETRACORE_CUDA
}

} // namespace Impl

std::string memorySpaceName (const void* ptr)
{
  using Impl::EMemoryType;

  const EMemoryType type = Impl::getCudaMemoryType (ptr);
  if (type == EMemoryType::CUDA_UVM) {
    return "CudaUVMSpace";
  }
  else if (type == EMemoryType::CUDA) {
    return "CudaSpace";
  }
  else if (type == EMemoryType::CUDA_HOST_PINNED) {
    return "CudaHostPinnedSpace";
  }
  else if (type == EMemoryType::HOST) {
    return "HostSpace";
  }
  else { // EMemoryType::ERROR
    return "ERROR";
  }
}

} // namespace Details
} // namespace Tpetra
