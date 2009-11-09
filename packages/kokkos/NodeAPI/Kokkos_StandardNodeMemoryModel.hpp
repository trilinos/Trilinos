#ifndef KOKKOS_STANDARD_NODE_MEMORY_MODEL_HPP_
#define KOKKOS_STANDARD_NODE_MEMORY_MODEL_HPP_

#include <stdlib.h>
#include <algorithm>
#include <Kokkos_NodeAPIConfigDefs.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>

namespace Kokkos {

  /*! A default implementation of the Node memory architecture for Node with a single memory space. */
  class StandardNodeMemoryModel {
    public:

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
      Teuchos::ArrayRCP<T> allocBuffer(size_t size) {
        Teuchos::ArrayRCP<T> ptr;
        if (size > 0) {
          ptr = Teuchos::arcp<T>(size);
        }
        return ptr;
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
      void copyFromBuffer(size_t size, const Teuchos::ArrayRCP<const T> &buffSrc, const Teuchos::ArrayView<T> &hostDest) {
        Teuchos::ArrayRCP<T> buffDest = Teuchos::arcpFromArrayView(hostDest);
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
      void copyToBuffer(size_t size, const Teuchos::ArrayView<const T> &hostSrc, const Teuchos::ArrayRCP<T> &buffDest) {
        Teuchos::ArrayRCP<const T> buffSrc = Teuchos::arcpFromArrayView(hostSrc);
        copyBuffers<T>(size,buffSrc,buffDest);
      }

      /*! \brief Copy data between buffers.
          
        @param[in]     size     The size of the copy, greater than zero.
        @param[in]     buffSrc  The source buffer, with length at least as large as \c size.
        @param[in,out] buffDest The destination buffer, with length at least as large as \c size.

        \post The data is guaranteed to have been copied before any other usage of buffSrc or buffDest occurs.
      */
      template <class T> inline
      void copyBuffers(size_t size, const Teuchos::ArrayRCP<const T> &buffSrc, const Teuchos::ArrayRCP<T> &buffDest) {
        Teuchos::ArrayView<const T> av_src = buffSrc(0,size);
        Teuchos::ArrayView<T>       av_dst = buffDest(0,size);
        std::copy(av_src.begin(),av_src.end(),av_dst.begin());
      }

      template <class T> inline
      Teuchos::ArrayRCP<const T> viewBuffer(size_t size, Teuchos::ArrayRCP<const T> buff) {
        return buff.persistingView(0,size);
      }

      template <class T> inline
      Teuchos::ArrayRCP<T> viewBufferNonConst(ReadWriteOption rw, size_t size, const Teuchos::ArrayRCP<T> &buff) {
        return buff.persistingView(0,size);
      }

      inline void readyBuffers(Teuchos::ArrayView<Teuchos::ArrayRCP<const char> > buffers, Teuchos::ArrayView<Teuchos::ArrayRCP<char> > ncBuffers) {
        (void)buffers;
        (void)ncBuffers;
      }


      //@}
  };

} // end of namespace Kokkos

#endif
