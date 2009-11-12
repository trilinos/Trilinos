#ifndef KOKKOS_CUDANODE_HPP_
#define KOKKOS_CUDANODE_HPP_

#include <cstring>
#include "Kokkos_NodeAPIConfigDefs.hpp"

// forward declarations of Teuchos classes
namespace Teuchos {
  class ParameterList;
  template <typename T> class ArrayRCP;
  template <typename T> class ArrayView;
}

namespace Kokkos {

class CUDANode {
  public:

    CUDANode(Teuchos::ParameterList &pl);

    ~CUDANode();

    //@{ Computational methods

    template <class WDP>
    static void parallel_for(int begin, int end, WDP wdp);

    template <class WDP>
    static typename WDP::ReductionType
    parallel_reduce(int begin, int end, WDP wd);

    //@} 

    //@{ Memory methods

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
    Teuchos::ArrayRCP<T> allocBuffer(size_t size);

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
    void copyFromBuffer(size_t size, const Teuchos::ArrayRCP<const T> &buffSrc, const Teuchos::ArrayView<T> &hostDest);

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
    void copyToBuffer(size_t size, const Teuchos::ArrayView<const T> &hostSrc, const Teuchos::ArrayRCP<T> &buffDest);

    /*! \brief Copy data between buffers.

      @param[in]     size     The size of the copy, greater than zero.
      @param[in]     buffSrc  The source buffer, with length at least as large as \c size.
      @param[in,out] buffDest The destination buffer, with length at least as large as \c size.

      \post The data is guaranteed to have been copied before any other usage of buffSrc or buffDest occurs.
     */
    template <class T> inline
    void copyBuffers(size_t size, const Teuchos::ArrayRCP<const T> &buffSrc, const Teuchos::ArrayRCP<T> &buffDest);

    template <class T> inline
    Teuchos::ArrayRCP<const T> viewBuffer(size_t size, Teuchos::ArrayRCP<const T> buff);

    template <class T> inline
    Teuchos::ArrayRCP<T> viewBufferNonConst(ReadWriteOption rw, size_t size, const Teuchos::ArrayRCP<T> &buff);

    void readyBuffers(Teuchos::ArrayView<Teuchos::ArrayRCP<const char> > buffers, Teuchos::ArrayView<Teuchos::ArrayRCP<char> > ncBuffers);

    //@} 

  private:
    template <class WDP, int FirstLevel>
    void call_reduce(int length, WDP wd, int threads, int blocks, void *d_blkpart);
    // numBlocks_ is 
    // - the number of blocks launched in a call to parallel_for()
    // - not used by reduce1D()
    int numBlocks_;
    // numThreads_ is required to be a power-of-two (our requirement) between 1 and 512 (CUDA's requirement). It is:
    // - the maximum number of threads used by reduce1D()
    // - the number of threads per block in a call to parallel_for()
    int numThreads_;
    // total global device memory, in bytes
    int totalMem_;
};

} // namespace Kokkos

#ifndef KOKKOS_CUDANODE_NO_IMPL
#include "Kokkos_CUDANodeImpl.hpp"
#endif

#endif
