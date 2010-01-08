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
#include "Kokkos_CUDA_util_inline_runtime.h"
#include "Kokkos_CUDANodeMemoryModel.hpp" // in case someone directly included this implementation file
#include "Kokkos_CUDANodeUtils.hpp"

namespace Kokkos {

  template <class T> inline
  Teuchos::ArrayRCP<T> 
  CUDANodeMemoryModel::allocBuffer(size_t size) {
    // FINISH: if possible, check that there is room; else, boot someone
    T * devptr = NULL;
    if (size > 0) {
      cutilSafeCallNoSync( cudaMalloc( (void**)&devptr, sizeof(T)*size ) );
    }
    CUDANodeDeallocator dealloc;
    const bool OwnsMem = true;
    const typename Teuchos::ArrayRCP<T>::Ordinal 
        LowerBound = 0,   
        UpperBound = size-1;
    Teuchos::ArrayRCP<T> buff(devptr,LowerBound,UpperBound,dealloc,OwnsMem);
    MARK_COMPUTE_BUFFER(buff);
    return buff;
  }

  template <class T> inline
  void CUDANodeMemoryModel::copyFromBuffer(size_t size, const Teuchos::ArrayRCP<const T> &buffSrc, const Teuchos::ArrayView<T> &hostDest) {
    CHECK_COMPUTE_BUFFER(buffSrc);
    TEST_FOR_EXCEPTION( buffSrc.size() < size || hostDest.size() < size, std::runtime_error,
        "CUDANode::copyFromBuffer: invalid copy.");
    cutilSafeCallNoSync( cudaMemcpy( hostDest.getRawPtr(), buffSrc.getRawPtr(), size*sizeof(T), cudaMemcpyDeviceToHost) );
  }

  template <class T> inline
  void CUDANodeMemoryModel::copyToBuffer(size_t size, const Teuchos::ArrayView<const T> &hostSrc, const Teuchos::ArrayRCP<T> &buffDest) {
    CHECK_COMPUTE_BUFFER(buffDest);
    TEST_FOR_EXCEPTION( hostSrc.size() < size || buffDest.size() < size, std::runtime_error,
        "CUDANode::copyFromBuffer: invalid copy.");
    cutilSafeCallNoSync( cudaMemcpy( buffDest.getRawPtr(), hostSrc.getRawPtr(), size*sizeof(T), cudaMemcpyHostToDevice) );
  }

  template <class T> inline
  void CUDANodeMemoryModel::copyBuffers(size_t size, const Teuchos::ArrayRCP<const T> &buffSrc, const Teuchos::ArrayRCP<T> &buffDest) {
    CHECK_COMPUTE_BUFFER(buffSrc);
    CHECK_COMPUTE_BUFFER(buffDest);
    TEST_FOR_EXCEPTION( buffSrc.size() < size || buffDest.size() < size, std::runtime_error,
        "CUDANode::copyFromBuffer: invalid copy.");
    cutilSafeCallNoSync( cudaMemcpy( buffDest.getRawPtr(), buffSrc.getRawPtr(), size*sizeof(T), cudaMemcpyDeviceToDevice) );
  }

  template <class T> inline
  Teuchos::ArrayRCP<const T> 
  CUDANodeMemoryModel::viewBuffer(size_t size, Teuchos::ArrayRCP<const T> buff) {
    CHECK_COMPUTE_BUFFER(buff);
    Teuchos::ArrayRCP<T> hostBuff;
    if (size != 0) {
      hostBuff = Teuchos::arcp<T>(size);
      this->template copyFromBuffer<T>(size,buff,hostBuff());
    }
    return hostBuff;
  }

  template <class T> inline
  Teuchos::ArrayRCP<T> 
  CUDANodeMemoryModel::viewBufferNonConst(ReadWriteOption rw, size_t size, const Teuchos::ArrayRCP<T> &buff) {
    CHECK_COMPUTE_BUFFER(buff);
    CUDANodeCopyBackDeallocator<T> dealloc(buff.getRawPtr(), size);
    Teuchos::ArrayRCP<T> hostBuff = dealloc.alloc();
    if (rw == ReadWrite) {
      this->template copyFromBuffer<T>(size, buff, hostBuff());
    }  
    // else rw == WriteOnly, and we need no copy
    return hostBuff;
  }

  inline void CUDANodeMemoryModel::readyBuffers(Teuchos::ArrayView<Teuchos::ArrayRCP<const char> > buffers, Teuchos::ArrayView<Teuchos::ArrayRCP<char> > ncBuffers) {
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
