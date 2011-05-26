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
      TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error,
          "Kokkos::CUDANodeMemoryModel::allocBuffer(): cudaMalloc() returned error: "
          << cudaGetErrorString(err) 
          );
#ifdef HAVE_KOKKOS_CUDA_NODE_MEMORY_PROFILING
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
    TEST_FOR_EXCEPTION( (size_t)buffSrc.size() < size || (size_t)hostDest.size() < size, std::runtime_error,
        "CUDANodeMemoryModel::copyFromBuffer: invalid copy.");
#ifdef HAVE_KOKKOS_CUDA_NODE_MEMORY_PROFILING
    ++numCopiesD2H_;
    bytesCopiedD2H_ += size*sizeof(T);
#endif
#ifdef HAVE_KOKKOS_CUDA_NODE_MEMORY_TRACE
    std::cerr << "copyFromBuffer<" << Teuchos::TypeNameTraits<T>::name() << "> of size " << sizeof(T) * size << std::endl;
#endif
    cudaError_t err = cudaMemcpy( hostDest.getRawPtr(), buffSrc.getRawPtr(), size*sizeof(T), cudaMemcpyDeviceToHost);
    TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error,
        "Kokkos::CUDANodeMemoryModel::copyFromBuffer(): cudaMemcpy() returned error: "
        << cudaGetErrorString(err) 
        );
  }

  template <class T> inline
  void CUDANodeMemoryModel::copyToBuffer(size_t size, const ArrayView<const T> &hostSrc, const ArrayRCP<T> &buffDest) {
    CHECK_COMPUTE_BUFFER(buffDest);
    TEST_FOR_EXCEPTION( hostSrc.size() < size, std::runtime_error, "CUDANodeMemoryModel::copyToBuffer: invalid copy.");
    TEST_FOR_EXCEPTION( buffDest.size() < size, std::runtime_error, "CUDANodeMemoryModel::copyToBuffer: invalid copy.");
#ifdef HAVE_KOKKOS_CUDA_NODE_MEMORY_PROFILING
    ++numCopiesH2D_;
    bytesCopiedH2D_ += size*sizeof(T);
#endif
#ifdef HAVE_KOKKOS_CUDA_NODE_MEMORY_TRACE
    std::cerr << "copyToBuffer<" << Teuchos::TypeNameTraits<T>::name() << "> of size " << sizeof(T) * size << std::endl;
#endif
    cudaError_t err = cudaMemcpy( buffDest.getRawPtr(), hostSrc.getRawPtr(), size*sizeof(T), cudaMemcpyHostToDevice);
    TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error,
        "Kokkos::CUDANodeMemoryModel::copyToBuffer(): cudaMemcpy() returned error: "
        << cudaGetErrorString(err) 
        );
  }

  template <class T> inline
  void CUDANodeMemoryModel::copyBuffers(size_t size, const ArrayRCP<const T> &buffSrc, const ArrayRCP<T> &buffDest) {
    CHECK_COMPUTE_BUFFER(buffSrc);
    CHECK_COMPUTE_BUFFER(buffDest);
    TEST_FOR_EXCEPTION( buffSrc.size() < size || buffDest.size() < size, std::runtime_error,
        "CUDANodeMemoryModel::copyBuffers: invalid copy.");
#ifdef HAVE_KOKKOS_CUDA_NODE_MEMORY_PROFILING
    ++numCopiesD2D_;
    bytesCopiedD2D_ += size*sizeof(T);
#endif
#ifdef HAVE_KOKKOS_CUDA_NODE_MEMORY_TRACE
    std::cerr << "copyBuffers<" << Teuchos::TypeNameTraits<T>::name() << "> of size " << sizeof(T) * size << std::endl;
#endif
    cudaError_t err = cudaMemcpy( buffDest.getRawPtr(), buffSrc.getRawPtr(), size*sizeof(T), cudaMemcpyDeviceToDevice);
    TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error,
        "Kokkos::CUDANodeMemoryModel::copyBuffers(): cudaMemcpy() returned error: "
        << cudaGetErrorString(err) 
        );
  }

  template <class T> inline
  ArrayRCP<const T> 
  CUDANodeMemoryModel::viewBuffer(size_t size, ArrayRCP<const T> buff) {
    CHECK_COMPUTE_BUFFER(buff);
    ArrayRCP<T> hostBuff;
    if (size != 0) {
      hostBuff = arcp<T>(size);
#ifdef HAVE_KOKKOS_CUDA_NODE_MEMORY_TRACE
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
    // create a copy-back deallocator that copies back to "buff"
    CUDANodeCopyBackDeallocator<T> dealloc(buff.persistingView(0,size), rcpFromRef(*this));
    // it allocates a host buffer with the appropriate deallocator embedded
    ArrayRCP<T> hostBuff = dealloc.alloc();
    if (rw == ReadWrite) {
#ifdef HAVE_KOKKOS_CUDA_NODE_MEMORY_TRACE
      std::cerr << "viewBufferNonConst(ReadWrite) -> ";
#endif
      this->template copyFromBuffer<T>(size, buff, hostBuff());
    }  
    else {
#ifdef HAVE_KOKKOS_CUDA_NODE_MEMORY_TRACE
      std::cerr << "viewBufferNonConst(WriteOnly)" << std::endl;
#endif
    }
    // else rw == WriteOnly, and we need no copy
    return hostBuff;
  }

  inline void CUDANodeMemoryModel::readyBuffers(ArrayView<ArrayRCP<const char> > buffers, ArrayView<ArrayRCP<char> > ncBuffers) {
#ifdef HAVE_KOKKOS_DEBUG
    for (size_t i=0; i < buffers.size(); ++i) {
      CHECK_COMPUTE_BUFFER(buffers[i]);
    }
    for (size_t i=0; i < ncBuffers.size(); ++i) {
      CHECK_COMPUTE_BUFFER(ncBuffers[i]);
    }
#endif
    (void)buffers;
    (void)ncBuffers;
  }

} // end of namespace Kokkos

#endif
