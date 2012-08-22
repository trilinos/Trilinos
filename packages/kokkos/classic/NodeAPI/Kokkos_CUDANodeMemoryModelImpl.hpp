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

#ifndef KOKKOS_CUDA_NODE_MEMORY_MODEL_IMPL_HPP_
#define KOKKOS_CUDA_NODE_MEMORY_MODEL_IMPL_HPP_

#include <cuda.h>
#include <cuda_runtime.h>

#include <iostream>
#include <cstdlib>
#include <stdexcept>

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include "Kokkos_NodeAPIConfigDefs.hpp"
#include "Kokkos_BufferMacros.hpp"
#include "Kokkos_CUDANodeMemoryModel.hpp" // in case someone directly included this implementation file
#include "Kokkos_CUDANodeUtils.hpp"

namespace Kokkos {

  template <class T> inline
  ArrayRCP<T> 
  CUDANodeMemoryModel::allocBuffer(size_t size) {
    // FINISH: if possible, check that there is room; else, boot someone
    T * devptr = NULL;
    const size_t sizeInBytes = sizeof(T)*size;
    if (size > 0) {
      cudaError_t err = cudaMalloc( (void**)&devptr, sizeInBytes );
      TEUCHOS_TEST_FOR_EXCEPTION( err != cudaSuccess, std::runtime_error,
        "Kokkos::CUDANodeMemoryModel::allocBuffer<" 
        << Teuchos::TypeNameTraits<T>::name () << ">: cudaMalloc() returned "
        "error: " << cudaGetErrorString (err) 
        );
#ifdef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_PROFILING
      allocSize_ += sizeInBytes;
#endif
    }
    CUDANodeDeallocator dealloc(sizeInBytes,rcpFromRef(*this));
    const bool OwnsMem = true;
    ArrayRCP<T> buff = arcp<T>(devptr,0,size,dealloc,OwnsMem);
    MARK_COMPUTE_BUFFER(buff);
    return buff;
  }

  template <class T> inline
  void CUDANodeMemoryModel::copyFromBuffer(size_t size, const ArrayRCP<const T> &buffSrc, const ArrayView<T> &hostDest) {
    CHECK_COMPUTE_BUFFER(buffSrc);
    TEUCHOS_TEST_FOR_EXCEPTION( (size_t)buffSrc.size() < size, std::runtime_error,
      "CUDANodeMemoryModel::copyFromBuffer<" 
      << Teuchos::TypeNameTraits<T>::name () 
      << ">: invalid copy.  Device source buffer has size " << buffSrc.size () 
      << ", which is less than the requested copy size " << size << ".");
    TEUCHOS_TEST_FOR_EXCEPTION( (size_t)hostDest.size() < size, std::runtime_error,
      "CUDANodeMemoryModel::copyFromBuffer<" 
      << Teuchos::TypeNameTraits<T>::name () 
      << ">: invalid copy.  Host destination buffer has size " << hostDest.size () 
      << ", which is less than the requested copy size " << size << ".");
#ifdef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_PROFILING
    ++numCopiesD2H_;
    bytesCopiedD2H_ += size*sizeof(T);
#endif
#ifdef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_TRACE
    std::cerr << "copyFromBuffer<" << Teuchos::TypeNameTraits<T>::name() << "> of size " << sizeof(T) * size << std::endl;
#endif
    cudaError_t err = cudaMemcpy( hostDest.getRawPtr(), buffSrc.getRawPtr(), size*sizeof(T), cudaMemcpyDeviceToHost);
    TEUCHOS_TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error,
      "Kokkos::CUDANodeMemoryModel::copyFromBuffer<"
      << Teuchos::TypeNameTraits<T>::name () 
      << ">(): cudaMemcpy() returned error: " << cudaGetErrorString (err) 
      );
  }

  template <class T> inline
  void CUDANodeMemoryModel::copyToBuffer(size_t size, const ArrayView<const T> &hostSrc, const ArrayRCP<T> &buffDest) {
    CHECK_COMPUTE_BUFFER(buffDest);
    TEUCHOS_TEST_FOR_EXCEPTION( (size_t)buffDest.size() < size, std::runtime_error,
      "CUDANodeMemoryModel::copyToBuffer<" 
      << Teuchos::TypeNameTraits<T>::name () 
      << ">: invalid copy.  Device destination buffer has size " << buffDest.size () 
      << ", which is less than the requested copy size " << size << ".");
    TEUCHOS_TEST_FOR_EXCEPTION( (size_t)hostSrc.size() < size, std::runtime_error,
      "CUDANodeMemoryModel::copyToBuffer<" 
      << Teuchos::TypeNameTraits<T>::name () 
      << ">: invalid copy.  Host source buffer has size " << hostSrc.size () 
      << ", which is less than the requested copy size " << size << ".");
#ifdef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_PROFILING
    ++numCopiesH2D_;
    bytesCopiedH2D_ += size*sizeof(T);
#endif
#ifdef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_TRACE
    std::cerr << "copyToBuffer<" << Teuchos::TypeNameTraits<T>::name() << "> of size " << sizeof(T) * size << std::endl;
#endif
    cudaError_t err = cudaMemcpy( buffDest.getRawPtr(), hostSrc.getRawPtr(), size*sizeof(T), cudaMemcpyHostToDevice);
    TEUCHOS_TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error,
      "Kokkos::CUDANodeMemoryModel::copyToBuffer<"
      << Teuchos::TypeNameTraits<T>::name () 
      << ">(): cudaMemcpy() returned error: " << cudaGetErrorString (err)
      );
  }

  template <class T> inline
  void CUDANodeMemoryModel::copyBuffers(size_t size, const ArrayRCP<const T> &buffSrc, const ArrayRCP<T> &buffDest) {
    CHECK_COMPUTE_BUFFER(buffSrc);
    CHECK_COMPUTE_BUFFER(buffDest);

    TEUCHOS_TEST_FOR_EXCEPTION( (size_t)buffDest.size() < size, std::runtime_error,
      "CUDANodeMemoryModel::copyBuffers<" 
      << Teuchos::TypeNameTraits<T>::name () 
      << ">: invalid copy.  Device destination buffer has size " << buffDest.size () 
      << ", which is less than the requested copy size " << size << ".");
    TEUCHOS_TEST_FOR_EXCEPTION( (size_t)buffSrc.size() < size, std::runtime_error,
      "CUDANodeMemoryModel::copyBuffers<" 
      << Teuchos::TypeNameTraits<T>::name () 
      << ">: invalid copy.  Device source buffer has size " << buffSrc.size () 
      << ", which is less than the requested copy size " << size << ".");

#ifdef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_PROFILING
    ++numCopiesD2D_;
    bytesCopiedD2D_ += size*sizeof(T);
#endif
#ifdef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_TRACE
    std::cerr << "copyBuffers<" << Teuchos::TypeNameTraits<T>::name() << "> of size " << sizeof(T) * size << std::endl;
#endif
    cudaError_t err = cudaMemcpy( buffDest.getRawPtr(), buffSrc.getRawPtr(), size*sizeof(T), cudaMemcpyDeviceToDevice);
    TEUCHOS_TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error,
      "Kokkos::CUDANodeMemoryModel::copyBuffers<"
      << Teuchos::TypeNameTraits<T>::name () 
      << ">(): cudaMemcpy() returned error: " << cudaGetErrorString (err)
      );
  }

  template <class T> inline
  ArrayRCP<const T> 
  CUDANodeMemoryModel::viewBuffer(size_t size, ArrayRCP<const T> buff) {
    CHECK_COMPUTE_BUFFER(buff);
    ArrayRCP<T> hostBuff;
    if (size != 0) {
      hostBuff = arcp<T>(size);
#ifdef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_TRACE
      std::cerr << "viewBuffer() -> ";
#endif
      this->template copyFromBuffer<T>(size,buff,hostBuff());
    }
    return hostBuff;
  }

  template <class T> inline
  ArrayRCP<T> 
  CUDANodeMemoryModel::viewBufferNonConst(ReadWriteOption rw, size_t size, const ArrayRCP<T> &buff) {
    CHECK_COMPUTE_BUFFER(buff);
    // Create a copy-back deallocator that copies back to buff.
    CUDANodeCopyBackDeallocator<T> dealloc(buff.persistingView(0,size), rcpFromRef(*this));
    // It allocates a host buffer with the appropriate deallocator embedded.
    ArrayRCP<T> hostBuff = dealloc.alloc();
    if (rw == ReadWrite) {
#ifdef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_TRACE
      std::cerr << "viewBufferNonConst(ReadWrite) -> ";
#endif
      this->template copyFromBuffer<T>(size, buff, hostBuff());
    }  
    else {
#ifdef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_TRACE
      std::cerr << "viewBufferNonConst(WriteOnly)" << std::endl;
#endif
    }
    // else rw == WriteOnly, and we need no copy
    return hostBuff;
  }

  inline void CUDANodeMemoryModel::readyBuffers(ArrayView<ArrayRCP<const char> > buffers, ArrayView<ArrayRCP<char> > ncBuffers) {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
    for (size_t i=0; i < (size_t)buffers.size(); ++i) {
      CHECK_COMPUTE_BUFFER(buffers[i]);
    }
    for (size_t i=0; i < (size_t)ncBuffers.size(); ++i) {
      CHECK_COMPUTE_BUFFER(ncBuffers[i]);
    }
#endif
    (void)buffers;
    (void)ncBuffers;
  }

} // end of namespace Kokkos

#endif
