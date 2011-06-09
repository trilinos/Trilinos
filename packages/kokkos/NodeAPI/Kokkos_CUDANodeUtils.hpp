#ifndef KOKKOS_CUDANODEUTILS_HPP_
#define KOKKOS_CUDANODEUTILS_HPP_

#include <cuda.h>
#include <cuda_runtime.h>

#include "Kokkos_ConfigDefs.hpp"
#define KOKKOS_NO_INCLUDE_INSTANTIATIONS
#include "Kokkos_CUDANodeMemoryModel.hpp"

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>

namespace Kokkos {

  class CUDANodeDeallocator {
    public:
      CUDANodeDeallocator(size_t sizeInBytes, const RCP<CUDANodeMemoryModel> &node);
      void free(void *ptr);
    private:
#ifdef HAVE_KOKKOS_CUDA_NODE_MEMORY_PROFILING
      const RCP<CUDANodeMemoryModel> node_;
      const size_t allocSize_;
#endif
  };

  //! \class CUDANodeCopyBackDeallocator
  /*! \brief Allocator/deallocator with host/device copy-back capability.
  
    Allocates a segment of page-locked memory associated with CUDA
    device memory. Upon deallocation, performs a copy-back of the allocated host
    memory before the host memory is deallocated. This copy-back is only performed
    if the device buffer is still valid (i.e., it hasn't been deallocated).
  */
  template <class T>
  class CUDANodeCopyBackDeallocator {
    public:
      CUDANodeCopyBackDeallocator(const ArrayRCP<T> &buffer, const RCP<CUDANodeMemoryModel> &node);

      //! Allocate the buffer, returning a Teuchos::ArrayRCP of the requested type, with a copy-back to GPU memory occurring at deallocation.
      ArrayRCP<T> alloc()const ;

      void free(void *ptr) const;
    private:
      // we have to keep a copy of this ArrayRCP, to know whether the underlying memory was deleted
      const ArrayRCP<T> devbuf_;
      const RCP<CUDANodeMemoryModel> node_;
#ifdef HAVE_KOKKOS_DEBUG
      mutable T * originalHostPtr_;
#endif
  };

  template <class T>
  CUDANodeCopyBackDeallocator<T>::CUDANodeCopyBackDeallocator(const ArrayRCP<T> &buffer,   
                                                              const RCP<CUDANodeMemoryModel> &node)
  : devbuf_(buffer.create_weak())
  , node_(node)
  { 
#ifdef HAVE_KOKKOS_DEBUG
    TEST_FOR_EXCEPT(node_ == null);
    originalHostPtr_ = NULL;
#endif
  }

  template <class T>
  ArrayRCP<T>
  CUDANodeCopyBackDeallocator<T>::alloc() const {
#ifdef HAVE_KOKKOS_DEBUG
    TEST_FOR_EXCEPTION( originalHostPtr_ != NULL, std::runtime_error,
        Teuchos::typeName(*this) << "::alloc(): alloc() has already been called." );
#endif
    T *hostPtr = NULL;
    // alloc page-locked ("pinned") memory on the host
    cudaError_t err = cudaHostAlloc( (void**)&hostPtr, devbuf_.size()*sizeof(T), cudaHostAllocDefault);
    TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error,
        "Kokkos::CUDANodeCopyBackDeallocator::alloc(): cudaHostAlloc() returned error:\n"
        << cudaGetErrorString(err) 
    );
#ifdef HAVE_KOKKOS_DEBUG
    // save the allocated address for debug checking
    originalHostPtr_ = hostPtr; 
#endif
    // create an ARCP<T> owning this memory, with a copy of *this for the deallocator
    const bool OwnsMem = true;
    return arcp<T>( hostPtr, 0, devbuf_.size(), *this, OwnsMem );
  }

  template <class T>
  void CUDANodeCopyBackDeallocator<T>::free(void *hostPtr) const {
#ifdef HAVE_KOKKOS_DEBUG
    TEST_FOR_EXCEPTION( hostPtr != originalHostPtr_, std::logic_error,
        Teuchos::typeName(*this) << "::free(): pointer to free not consistent with originally allocated pointer." );
    originalHostPtr_ = NULL;
#endif
    // only perform the copy back if the device ptr is still valid
    if (devbuf_.is_valid_ptr()) {
      // create temporary ArrayView for use with copyToBuffer
      // we must disable the lookup, or a debug build of Teuchos will freak out
      ArrayView<const T> tmpav((const T*)hostPtr, devbuf_.size(), Teuchos::RCP_DISABLE_NODE_LOOKUP);
      node_->template copyToBuffer<T>(devbuf_.size(), tmpav, devbuf_);
    }
    cudaError_t err = cudaFreeHost( (void**)hostPtr );
    TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error,
        "Kokkos::CUDANodeCopyBackDeallocator::free(): cudaFreeHost() returned error:\n"
        << cudaGetErrorString(err) 
    );
    hostPtr = NULL;
  }

}

#endif // KOKKOS_CUDANODEUTILS_HPP_
