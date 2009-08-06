#ifndef KOKKOS_STANDARD_NODE_MEMORY_MODEL_HPP_
#define KOKKOS_STANDARD_NODE_MEMORY_MODEL_HPP_

#include <stdlib.h>
#include <algorithm>
#include <Kokkos_ConfigDefs.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>

namespace Kokkos {

  class StandardMemoryModel {
    public:

      //@{ Memory management

      template <class T> inline
      Teuchos::ArrayRCP<T> allocBuffer(size_type size) {
        if (size > 0) {
          return Teuchos::arcp<T>(size);
        }
        return Teuchos::null;
      }

      template <class T> inline
      void copyFromBuffer(const Teuchos::ArrayRCP<const T> &buffSrc, size_type src_offset, const Teuchos::ArrayView<T> &hostDest) {
        Teuchos::ArrayRCP<T> buffDest = Teuchos::arcpFromArrayView(hostDest);
        copyBuffers(buffDest.size(),buffSrc,src_offset,buffDest,0);
      }

      template <class T> inline
      void copyToBuffer(Teuchos::ArrayView<const T> hostSrc, Teuchos::ArrayRCP<T> buffDest, size_type dest_offset) {
        Teuchos::ArrayRCP<const T> buffSrc = Teuchos::arcpFromArrayView(hostSrc);
        copyBuffers<T>(buffSrc.size(),buffSrc,0,buffDest,dest_offset);
      }

      template <class T> inline
      void copyBuffers(size_type size, const Teuchos::ArrayRCP<const T> &buffSrc , size_type src_offset, 
                                       const Teuchos::ArrayRCP<      T> &buffDest, size_type dest_offset) {
        Teuchos::ArrayView<const T> av_src = buffSrc(src_offset,size);
        Teuchos::ArrayView<T>       av_dst = buffDest(dest_offset,size);
        std::copy(av_src.begin(),av_src.end(),av_dst.begin());
      }

      template <class T> inline
      Teuchos::ArrayRCP<const T> viewBufferConst(size_type size, Teuchos::ArrayRCP<const T> buff, size_type offset) {
        return buff.persistingView(offset,size);
      }

      template <class T> inline
      Teuchos::ArrayRCP<T> viewBuffer(bool writeOnly, size_type size, const Teuchos::ArrayRCP<T> &buff, size_type offset) {
        return buff.persistingView(offset,size);
      }

      void readyBuffers(const Teuchos::ArrayView<const void *> &buffers, Teuchos::ArrayView<void *> &ncBuffers) {
        (void)buffers;
        (void)ncBuffers;
      }

      //@}
  };

} // end of namespace Kokkos

#endif
