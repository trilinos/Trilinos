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
#ifdef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_PROFILING
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
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
      mutable T * originalHostPtr_;
#endif
  };

  template <class T>
  CUDANodeCopyBackDeallocator<T>::CUDANodeCopyBackDeallocator(const ArrayRCP<T> &buffer,
                                                              const RCP<CUDANodeMemoryModel> &node)
  : devbuf_(buffer.create_weak())
  , node_(node)
  {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(node_ == null);
    originalHostPtr_ = NULL;
#endif
  }

  template <class T>
  ArrayRCP<T>
  CUDANodeCopyBackDeallocator<T>::alloc() const {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( originalHostPtr_ != NULL, std::runtime_error,
        Teuchos::typeName(*this) << "::alloc(): alloc() has already been called." );
#endif
    T *hostPtr = NULL;
    // alloc page-locked ("pinned") memory on the host
    // TODO: review: instead of cudaHostAllocDefault, this might should be cudaHostAllocWriteCombined
    cudaError_t err = cudaHostAlloc( (void**)&hostPtr, devbuf_.size()*sizeof(T), cudaHostAllocDefault);
    TEUCHOS_TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error,
        "Kokkos::CUDANodeCopyBackDeallocator::alloc(): cudaHostAlloc() returned error:\n"
        << cudaGetErrorString(err)
    );
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
    // save the allocated address for debug checking
    originalHostPtr_ = hostPtr;
#endif
    // create an ARCP<T> owning this memory, with a copy of *this for the deallocator
    const bool OwnsMem = true;
    return arcp<T>( hostPtr, 0, devbuf_.size(), *this, OwnsMem );
  }

  template <class T>
  void CUDANodeCopyBackDeallocator<T>::free(void *hostPtr) const {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( hostPtr != originalHostPtr_, std::logic_error,
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
    cudaError_t err = cudaFreeHost( hostPtr );
    TEUCHOS_TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error,
        "Kokkos::CUDANodeCopyBackDeallocator::free(): cudaFreeHost() returned error:\n"
        << cudaGetErrorString(err)
    );
    hostPtr = NULL;
  }

  //! \class CUDANodeHostPinnedAllocator
  /*! \brief Allocator/deallocator with pinned-host capability.

    Allocates a segment of page-locked memory associated with CUDA
    device memory. Upon deallocation, delete it appropriately.
  */
  template <class T>
  class CUDANodeHostPinnedDeallocator {
  public:
    // Constructor.
    CUDANodeHostPinnedDeallocator();

    /// \brief Allocate a buffer of the requested size.
    ///
    /// \warning This method may be called at most once on a single
    ///   instance.  In a debug build, calling this method more than
    ///   once throws an exception.
    ///
    /// \param sz [in] Number of entries in the buffer.
    ///
    /// \return An ArrayRCP of the requested type and with at least sz
    ///   entries, which copies back to GPU memory at deallocation.
    ArrayRCP<T> alloc (const size_t sz) const ;

    /// \brief Deallocate the memory pointed to by ptr.
    ///
    /// \note In a debug build of Kokkos, ptr must be the same as the
    ///   memory allocated by alloc().
    void free (void *ptr) const;

    private:
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
    /// \brief The memory allocated by the last call to alloc().
    ///
    /// This is NULL if alloc() has not yet been called.
    mutable T* originalHostPtr_;
#endif // HAVE_KOKKOSCLASSIC_DEBUG
  };

  template <class T>
  CUDANodeHostPinnedDeallocator<T>::CUDANodeHostPinnedDeallocator()
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
  : originalHostPtr_(NULL)
#endif // HAVE_KOKKOSCLASSIC_DEBUG
  { }

  template <class T>
  ArrayRCP<T>
  CUDANodeHostPinnedDeallocator<T>::alloc(const size_t sz) const {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(originalHostPtr_ != NULL, std::runtime_error,
      Teuchos::typeName(*this) << "::alloc(): alloc() has already been called." );
#endif
    T *hostPtr = NULL;
    // alloc page-locked ("pinned") memory on the host
    cudaError_t err = cudaHostAlloc( (void**)&hostPtr, sz*sizeof(T), cudaHostAllocDefault);
    TEUCHOS_TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error,
        "Kokkos::CUDANodeHostPinnedDeallocator::alloc(): cudaHostAlloc() returned error:\n"
        << cudaGetErrorString(err)
    );
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
    // save the allocated address for debug checking
    originalHostPtr_ = hostPtr;
#endif
    // create an ARCP<T> owning this memory, with a copy of *this for the deallocator
    const bool OwnsMem = true;
    return arcp<T>( hostPtr, 0, sz, *this, OwnsMem );
  }

  template <class T>
  void CUDANodeHostPinnedDeallocator<T>::free(void *hostPtr) const {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( hostPtr != originalHostPtr_, std::logic_error,
        Teuchos::typeName(*this) << "::free(): pointer to free not consistent with originally allocated pointer." );
    originalHostPtr_ = NULL;
#endif
    cudaError_t err = cudaFreeHost( hostPtr );
    TEUCHOS_TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error,
        "Kokkos::CUDANodeHostPinnedDeallocator::free(): cudaFreeHost() returned error:\n"
        << cudaGetErrorString(err)
    );
  }

}

#endif // KOKKOS_CUDANODEUTILS_HPP_
