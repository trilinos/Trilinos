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

#ifndef KOKKOS_STANDARD_NODE_MEMORY_MODEL_HPP_
#define KOKKOS_STANDARD_NODE_MEMORY_MODEL_HPP_

#include "Kokkos_NodeAPIConfigDefs.hpp"
#include "Kokkos_BufferMacros.hpp"

#include <cstdlib>
#include <algorithm>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>

namespace Kokkos {

  /** \brief A default implementation of the Node memory architecture for a single memory space allocated by standard library calls.
      \ingroup kokkos_node_api
   */
  class StandardNodeMemoryModel {
    public:
      //! Indicates that parallel buffers allocated by this node are available for use on the host thread.
      static const bool isHostNode = true;
      static const bool isCUDANode = false;

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
      ArrayRCP<T> allocBuffer(size_t size) {
        ArrayRCP<T> buff;
        if (size > 0) {
          buff = arcp<T>(size);
        }
        if (isHostNode == false) {
          MARK_COMPUTE_BUFFER(buff);
        }
        return buff;
      }

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
      void copyFromBuffer(size_t size, const ArrayRCP<const T> &buffSrc, const ArrayView<T> &hostDest) {
        if (isHostNode == false) {
          CHECK_COMPUTE_BUFFER(buffSrc);
        }
        ArrayRCP<T> buffDest = arcpFromArrayView(hostDest);
        copyBuffers(size,buffSrc,buffDest);
      }

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
      void copyToBuffer(size_t size, const ArrayView<const T> &hostSrc, const ArrayRCP<T> &buffDest) {
        if (isHostNode == false) {
          CHECK_COMPUTE_BUFFER(buffDest);
        }
        ArrayRCP<const T> buffSrc = arcpFromArrayView(hostSrc);
        copyBuffers<T>(size,buffSrc,buffDest);
      }

      /*! \brief Copy data between buffers.
          
        @param[in]     size     The size of the copy, greater than zero.
        @param[in]     buffSrc  The source buffer, with length at least as large as \c size.
        @param[in,out] buffDest The destination buffer, with length at least as large as \c size.

        \post The data is guaranteed to have been copied before any other usage of buffSrc or buffDest occurs.
      */
      template <class T> inline
      void copyBuffers(size_t size, const ArrayRCP<const T> &buffSrc, const ArrayRCP<T> &buffDest) {
        if (isHostNode == false) {
          CHECK_COMPUTE_BUFFER(buffSrc);
          CHECK_COMPUTE_BUFFER(buffDest);
        }
        ArrayView<const T> av_src = buffSrc(0,size);
        ArrayView<T>       av_dst = buffDest(0,size);
        std::copy(av_src.begin(),av_src.end(),av_dst.begin());
      }

      //! \brief Return a const view of a buffer for use on the host.
      template <class T> inline
      ArrayRCP<const T> viewBuffer(size_t size, ArrayRCP<const T> buff) {
        if (isHostNode == false) {
          CHECK_COMPUTE_BUFFER(buff);
        }
        return buff.persistingView(0,size);
      }

      //! \brief Return a non-const view of a buffer for use on the host.
      template <class T> inline
      ArrayRCP<T> viewBufferNonConst(ReadWriteOption rw, size_t size, const ArrayRCP<T> &buff) {
	(void) rw; // Silence "unused parameter" compiler warning
        if (isHostNode == false) {
          CHECK_COMPUTE_BUFFER(buff);
        }
        return buff.persistingView(0,size);
      }

      inline void readyBuffers(ArrayView<ArrayRCP<const char> > buffers, ArrayView<ArrayRCP<char> > ncBuffers) {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
        if (isHostNode == false) {
          for (size_t i=0; i < (size_t)buffers.size(); ++i) {
            CHECK_COMPUTE_BUFFER(buffers[i]);
          }
          for (size_t i=0; i < (size_t)ncBuffers.size(); ++i) {
            CHECK_COMPUTE_BUFFER(ncBuffers[i]);
          }
        }
#endif
        (void)buffers;
        (void)ncBuffers;
      }


      //@}
  };

} // end of namespace Kokkos

#endif
