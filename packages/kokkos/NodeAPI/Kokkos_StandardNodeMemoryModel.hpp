#ifndef KOKKOS_STANDARD_NODE_MEMORY_MODEL_HPP_
#define KOKKOS_STANDARD_NODE_MEMORY_MODEL_HPP_

#include <stdlib.h>
#include <algorithm>

namespace Kokkos {

class StandardMemoryModel {
  public:
    template <class T>
    struct buffer {
      typedef T* buffer_t;
    };

    //@{ Memory management

    template <class T> inline
    typename buffer<T>::buffer_t allocBuffer(unsigned int size) 
    {
      if (size > 0) {
        return (T *)malloc(sizeof(T)*size);
      }
      return 0;
    }

    template <class T> inline 
    void freeBuffer(typename buffer<T>::buffer_t buff) {
      cfree(buff);
    }

    template <class T> inline
    void copyFromBuffer(unsigned int size, typename buffer<const T>::buffer_t buff, unsigned int offset, T *ptr) {
      copyBuffers(size,buff,offset,ptr,0);
    }

    template <class T> inline
    void copyToBuffer(unsigned int size, const T *ptr, typename buffer<T>::buffer_t buff, unsigned int offset) {
      copyBuffers(size,ptr,0,buff,offset);
    }

    template <class T> inline
    void copyBuffers(unsigned int size, typename buffer<const T>::buffer_t src, unsigned int src_offset, 
                                       typename buffer<T>::buffer_t dest, unsigned int dest_offset) 
    {
      std::copy(src+src_offset,src+src_offset+size, dest+dest_offset);
    }

    template <class T> inline
    const T * viewBufferConst(unsigned int size, typename buffer<const T>::buffer_t buff, unsigned int offset) {
      return buff+offset;
    }

    template <class T> inline
    T * viewBuffer(bool writeOnly, unsigned int size, typename buffer<T>::buffer_t buff, unsigned int offset) {
      return buff+offset;
    }

    template <class T> inline
    void releaseView(T * /*view_ptr*/) {}

    template <class T> inline
    void releaseView(const T * /*view_ptr*/) {}

    void readyBuffers(const void * const * /*buffers*/, unsigned int /*numBuffers*/,
                            void * const * /*buffers*/, unsigned int /*numBuffers*/) {}

  private:
    void cfree(void *ptr)       {free(ptr);}
    void cfree(const void *ptr) {free(const_cast<void *>(ptr));}

    //@}
};

}

#endif
