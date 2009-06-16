#ifndef KOKKOS_CUDANODE_HPP_
#define KOKKOS_CUDANODE_HPP_

#include <map>

// forward declaration
class CUDANode;

template <class T>
struct CUDABuffer;

#define CONVERT_TO_CUDABUFFER_VOID \
    operator CUDABuffer<void>() const { \
      return CUDABuffer<void>(length_,flags_,devc_ptr_,host_ptr_); \
    }

#define CONVERT_TO_CUDABUFFER_CONSTVOID \
    operator CUDABuffer<const void>() const { \
      return CUDABuffer<const void>(length_,flags_,devc_ptr_,host_ptr_);  \
    }


template <>
struct CUDABuffer<const void> {
  friend class CUDANode;
  template <class T2>
  friend class CUDABuffer;
  public:
    CUDABuffer() { length_ = 0; flags_ = 0; devc_ptr_ = 0; host_ptr_ = 0; }
  private:
    CUDABuffer(int length, char *flags, const void *dptr, const void *hptr)
    {
      length_ = length;
      flags_ = flags;
      devc_ptr_ = dptr;
      host_ptr_ = hptr;
    }
    int length_;
    char *flags_;
    const void *devc_ptr_, *host_ptr_;
};

template <>
struct CUDABuffer<void> {
  friend class CUDANode;
  template <class T2>
  friend class CUDABuffer;
  public:
    CUDABuffer() { length_ = 0; flags_ = 0; devc_ptr_ = 0; host_ptr_ = 0; }
    CONVERT_TO_CUDABUFFER_CONSTVOID
  private:
    CUDABuffer(int length, char *flags, const void *dptr, const void *hptr)
    {
      length_ = length;
      flags_ = flags;
      devc_ptr_ = dptr;
      host_ptr_ = hptr;
    }
    int length_;
    char *flags_;
    const void *devc_ptr_, *host_ptr_;
};

template <class T>
struct CUDABuffer<const T> {
  friend class CUDANode;
  template <class T2>
  friend class CUDABuffer;
  public:
    CUDABuffer() { length_ = 0; flags_  = 0; devc_ptr_ = 0; host_ptr_ = 0; }
    CONVERT_TO_CUDABUFFER_VOID
    CONVERT_TO_CUDABUFFER_CONSTVOID
    const T & operator[](int i) const {return devc_ptr_[i];}
  private:
    CUDABuffer(int length, char *flags, const T *dptr, const T *hptr)
    {
      length_ = length;
      flags_ = flags;
      devc_ptr_ = dptr;
      host_ptr_ = hptr;
    }
    int length_;
    char *flags_;
    const T *devc_ptr_, *host_ptr_;
};

template <class T>
struct CUDABuffer {
  friend class CUDANode;
  public:
    CUDABuffer() { length_ = 0; flags_  = 0; devc_ptr_ = 0; host_ptr_ = 0; }
    CONVERT_TO_CUDABUFFER_VOID
    CONVERT_TO_CUDABUFFER_CONSTVOID
    operator CUDABuffer<const T>() const {
      return CUDABuffer<const T>(length_,flags_,devc_ptr_,host_ptr_); 
    }
    T & operator[](int i) const {return devc_ptr_[i];}
  private:
    CUDABuffer(int length, char *flags, T *dptr, T *hptr)
    {
      length_ = length;
      flags_ = flags;
      devc_ptr_ = dptr;
      host_ptr_ = hptr;
    }
    int length_;
    char *flags_;
    T *devc_ptr_, *host_ptr_;
};

class CUDANode {
  public:
    template <class T>
    struct buffer {
      typedef CUDABuffer<T> buffer_t;
    };

    CUDANode(int device = 0, int numBlocks = -1, int numThreads = 256, int verbose = 1);

    ~CUDANode();

    //@{ Computational methods

    template <class WDP>
    void execute1D(int length, WDP wdp);

    template <class WDP>
    void reduce1D(int length, WDP &wd);

    //@} 

    //@{ Memory methods

    template <class T>
    CUDABuffer<T> allocBuffer(unsigned int size);

    template <class T>
    void freeBuffer(CUDABuffer<T> buff);

    template <class T>
    void copyFromBuffer(unsigned int size, typename buffer<const T>::buffer_t buff, unsigned int offset, T *ptr);

    template <class T>
    void copyToBuffer(unsigned int size, const T *ptr, typename buffer<T>::buffer_t buff, unsigned int offset);

    template <class T>
    void copyBuffers(unsigned int size, typename buffer<const T>::buffer_t src , unsigned int  src_offset, 
                                        typename buffer<      T>::buffer_t dest, unsigned int dest_offset);

    template <class T>
    const T * viewBufferConst(unsigned int size, typename buffer<const T>::buffer_t buff, unsigned int offset);

    template <class T>
    T * viewBuffer(bool writeOnly, unsigned int size, typename buffer<T>::buffer_t buff, unsigned int offset);

    template <class T>
    void releaseView(const T * view_ptr);

    void readyBuffers(CUDABuffer<const void> * const cBuffers,  unsigned int numConstBuffers,
                      CUDABuffer<      void> * const ncBuffers, unsigned int numNonConstBuffers);

    //@} 

  private:
    static std::map<const void *, CUDABuffer<const void> > out_views_;
    static const char CLEAN      = 0,
                      DIRTY_DEVC = 1,
                      DIRTY_HOST = 2;
    template <class WDP, int FirstLevel>
    void call_reduce(int length, WDP wd, int threads, int blocks, void *d_blkpart);
    // numBlocks_ is 
    // - the number of blocks launched in a call to execute1D()
    // - not used by reduce1D()
    int numBlocks_;
    // numThreads_ is required to be a power-of-two (our requirement) between 1 and 512 (CUDA's requirement). It is:
    // - the maximum number of threads used by reduce1D()
    // - the number of threads per block in a call to execute1D()
    int numThreads_;
    // total global device memory, in bytes
    int totalMem_;
};

#endif
