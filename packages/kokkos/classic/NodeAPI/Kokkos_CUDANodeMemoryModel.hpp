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

#ifndef KOKKOS_CUDA_NODE_MEMORY_MODEL_HPP_
#define KOKKOS_CUDA_NODE_MEMORY_MODEL_HPP_

#include "Kokkos_NodeAPIConfigDefs.hpp"
#include <Teuchos_RCP.hpp>

// forward declarations of Teuchos classes
namespace std {
  template <typename CharT> class char_traits;
}
namespace Teuchos {
  class ParameterList;
  template <typename T> class ArrayRCP;
  template <typename T> class ArrayView;
  template <typename T> class RCP;
  template <typename CharT, typename Traits> class basic_FancyOStream;
  typedef basic_FancyOStream<char, std::char_traits<char> > FancyOStream;
}

namespace Kokkos {

  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;

  /** \brief A default implementation of the Node memory architecture for Node with a distinct device memory space allocated by the CUDA runtime.
      \ingroup kokkos_node_api
   */
  class CUDANodeMemoryModel {
    public:
      //! Indicates that parallel buffers allocated by this node are not available for use on the host thread.
      static const bool isHostNode = false;
      static const bool isCUDANode = true;

      //@{ Default Constructor 

      //! Default constructor.
      CUDANodeMemoryModel();

      //@}

      //@{ Memory management

      /*! \brief Allocate a parallel buffer, returning it as a pointer ecnapsulated in an ArrayRCP.

          Dereferencing the returned ArrayRCP or its underlying pointer in general results in undefined 
          behavior outside of parallel computations.

          The buffer will be automatically freed by the Node when no more references remain.

          @tparam T The data type of the allocate buffer. This is used to perform alignment and determine the number of bytes to allocate.
          @param[in] size The size requested for the parallel buffer, greater than zero.

          \post The method will return an ArrayRCP encapsulating a pointer. The underlying pointer may be used in parallel computation routines, 
                and is guaranteed to have size large enough to reference \c size number of entries of type \c T.
      */
      template <class T> inline
      ArrayRCP<T> allocBuffer(size_t size);

      /*! \brief Copy data to host memory from a parallel buffer.

          @param[in] size       The number of entries to copy from \c buffSrc to \c hostDest.
          @param[in] buffSrc    The parallel buffer from which to copy.
          @param[out] hostDest  The location in host memory where the data from \c buffSrc is copied to.

          \pre  \c size is non-negative.
          \pre  \c buffSrc has length at least <tt>size</tt>.
          \pre  \c hostDest has length equal to \c size.
          \post On return, entries in the range <tt>[0 , size)</tt> of \c buffSrc have been copied to \c hostDest entries in the range <tt>[0 , size)</tt>.
      */
      template <class T> inline
      void copyFromBuffer(size_t size, const ArrayRCP<const T> &buffSrc, const ArrayView<T> &hostDest);

      /*! \brief Copy data to host memory from a parallel buffer.

          @param[in]  size        The number of entries to copy from \c hostSrc to \c buffDest.
          @param[in]  hostSrc     The location in host memory from where the data is copied.
          @param[out] buffDest    The parallel buffer to which the data is copied.

          \pre  \c size is non-negative.
          \pre  \c hostSrc has length equal to \c size.
          \pre  \c buffSrc has length at least <tt>size</tt>.
          \post On return, entries in the range <tt>[0 , size)</tt> of \c hostSrc are allowed to be written to. The data is guaranteed to be present in \c buffDest before it is used in a parallel computation.
      */
      template <class T> inline
      void copyToBuffer(size_t size, const ArrayView<const T> &hostSrc, const ArrayRCP<T> &buffDest);

      /*! \brief Copy data between buffers.
          
        @param[in]     size     The size of the copy, greater than zero.
        @param[in]     buffSrc  The source buffer, with length at least as large as \c size.
        @param[in,out] buffDest The destination buffer, with length at least as large as \c size.

        \post The data is guaranteed to have been copied before any other usage of buffSrc or buffDest occurs.
      */
      template <class T> inline
      void copyBuffers(size_t size, const ArrayRCP<const T> &buffSrc, const ArrayRCP<T> &buffDest);

      /*! \brief Return a const view of a buffer for use on the host.
          This creates a \c const view of length \c size, constituting the first \c size entries of \c buff, as they exist at the time of view creation.
          host memory allocated for the creation of this view is automatically deleted when no refences to the view remain. 
          \pre <tt>buff.size() >= size</tt>
       */
      template <class T> inline
      ArrayRCP<const T> viewBuffer(size_t size, ArrayRCP<const T> buff);

      /*! \brief Return a non-const view of a buffer for use on the host.

          \param[in] rw Specifies Kokkos::ReadWrite or Kokkos::WriteOnly. If Kokkos::WriteOnly, the contents of the view are undefined when it is created and must 
          be initialized on the host. However, this prevents the potential need for a copy from device to host memory needed to set the view 
          values as when Kokkos::ReadWrite is specified.
          
          This creates a view of length \c size, constituting the first \c size entries of \c buff, as they exist at the time of view creation.

          A non-const view permits changes, which must be copied back to the buffer. This does not occur until all references to the view are deleted.
          If the buffer is deallocated before the view is deleted, then the copy-back does not occur.

          \pre <tt>buff.size() >= size</tt>
       */
      template <class T> inline
      ArrayRCP<T> viewBufferNonConst(ReadWriteOption rw, size_t size, const ArrayRCP<T> &buff);

      inline void readyBuffers(ArrayView<ArrayRCP<const char> > buffers, ArrayView<ArrayRCP<char> > ncBuffers);

      //@}

      //@{ Book-keeping information

      //! \brief Print some statistics regarding node allocation and memory transfer.
      void printStatistics(const RCP< Teuchos::FancyOStream > &os) const;
      
      //! \brief Clear all statistics on memory transfer.
      void clearStatistics();

      //@}

    public:
      size_t allocSize_;
      size_t numCopiesD2H_, numCopiesH2D_, numCopiesD2D_;
      size_t bytesCopiedD2H_, bytesCopiedH2D_, bytesCopiedD2D_;
  };

} // end of namespace Kokkos

#ifndef KOKKOS_NO_INCLUDE_INSTANTIATIONS
#include "Kokkos_CUDANodeMemoryModelImpl.hpp"
#endif

#endif
