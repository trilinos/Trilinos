#ifndef KOKKOS_CUDANODE_IMPL_HPP_
#define KOKKOS_CUDANODE_IMPL_HPP_

#include "Kokkos_CUDANode.hpp"
#include "Kokkos_CUDANodeUtils.hpp"

#include <cuda.h>
#include <cuda_runtime.h>
#include "Kokkos_CUDA_util_inline_runtime.h"
#include "Kokkos_BufferMacros.hpp"

#include <iostream>
#include <cstdlib>
#include <stdexcept>

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ParameterList.hpp>

namespace Kokkos {

  template <class T> inline
  Teuchos::ArrayRCP<T> CUDANode::allocBuffer(size_t size) {
    // FINISH: if possible, check that there is room; else, boot someone
    T * devptr = NULL;
    if (size > 0) {
      cutilSafeCallNoSync( cudaMalloc( (void**)&devptr, sizeof(T)*size ) );
    }
    CUDANodeDeallocator dealloc;
    Teuchos::ArrayRCP<T> buff(devptr,0,size,dealloc,true);
    MARK_COMPUTE_BUFFER(buff);
    return buff;
  }

  template <class T>
  void CUDANode::copyFromBuffer(size_t size, const Teuchos::ArrayRCP<const T> &buffSrc, const Teuchos::ArrayView<T> &hostDest)
  {
    CHECK_COMPUTE_BUFFER(buffSrc);
    TEST_FOR_EXCEPTION( buffSrc.size() < size || hostDest.size() < size, std::runtime_error,
        "CUDANode::copyFromBuffer: invalid copy.");
    cutilSafeCallNoSync( cudaMemcpy( hostDest.getRawPtr(), buffSrc.getRawPtr(), size*sizeof(T), cudaMemcpyDeviceToHost) );
  }

  template <class T>
  void CUDANode::copyToBuffer(size_t size, const Teuchos::ArrayView<const T> &hostSrc, const Teuchos::ArrayRCP<T> &buffDest)
  {
    CHECK_COMPUTE_BUFFER(buffDest);
    TEST_FOR_EXCEPTION( hostSrc.size() < size || buffDest.size() < size, std::runtime_error,
        "CUDANode::copyFromBuffer: invalid copy.");
    cutilSafeCallNoSync( cudaMemcpy( buffDest.getRawPtr(), hostSrc.getRawPtr(), size*sizeof(T), cudaMemcpyHostToDevice) );
  }

  template <class T>
  void CUDANode::copyBuffers(size_t size, const Teuchos::ArrayRCP<const T> &buffSrc, const Teuchos::ArrayRCP<T> &buffDest) 
  {
    CHECK_COMPUTE_BUFFER(buffSrc);
    CHECK_COMPUTE_BUFFER(buffDest);
    TEST_FOR_EXCEPTION( buffSrc.size() < size || buffDest.size() < size, std::runtime_error,
        "CUDANode::copyFromBuffer: invalid copy.");
    cutilSafeCallNoSync( cudaMemcpy( buffDest.getRawPtr(), buffSrc.getRawPtr(), size*sizeof(T), cudaMemcpyDeviceToDevice) );
  }

  template <class T> inline
  Teuchos::ArrayRCP<const T> CUDANode::viewBuffer(size_t size, Teuchos::ArrayRCP<const T> buff) 
  {
    CHECK_COMPUTE_BUFFER(buff);
    TEST_FOR_EXCEPT(true);
  }

  template <class T> inline
  Teuchos::ArrayRCP<T> CUDANode::viewBufferNonConst(ReadWriteOption rw, size_t size, const Teuchos::ArrayRCP<T> &buff) 
  {
    CHECK_COMPUTE_BUFFER(buff);
    TEST_FOR_EXCEPT(true);
  }

}

#endif
