#ifndef STK_ALIGNED_ALLOCATOR_HPP
#define STK_ALIGNED_ALLOCATOR_HPP

#if defined(__APPLE__)
#include <stdlib.h>
#else
#include <stdlib.h>
#include <malloc.h>
#endif

#include <memory>

namespace non_std {

template <class T, size_t Alignment>
struct AlignedAllocator
  : public std::allocator<T> 
{
  typedef typename std::allocator<T>::size_type size_type;
  typedef typename std::allocator<T>::pointer pointer;
  typedef typename std::allocator<T>::const_pointer const_pointer;

  template <class U>
  struct rebind { typedef AlignedAllocator<U,Alignment> other; };

  AlignedAllocator() throw() { }

  AlignedAllocator(const AlignedAllocator& other) throw()
    : std::allocator<T>(other) { }

  template <class U>
  AlignedAllocator(const AlignedAllocator<U,Alignment>&) throw() { }

  ~AlignedAllocator() throw() { }

  inline pointer allocate(size_type n) {
    return allocate(n, const_pointer(0));
  }
    
  inline pointer allocate(size_type n, const_pointer ) {
    void *p;
    if ( posix_memalign(&p, Alignment, n*sizeof(T)) != 0 ) { p = nullptr; throw std::bad_alloc(); }
    return static_cast<pointer>(p);
  }

  inline void deallocate(pointer p, size_type ) {
    free(p);
  }
};

template <class T1, size_t A1, class T2, size_t A2> inline
bool operator == (const AlignedAllocator<T1,A1> &, const AlignedAllocator<T2,A2> &) {
  return true;
}

template <class T1, size_t A1, class T2, size_t A2> inline
bool operator != (const AlignedAllocator<T1,A1> &, const AlignedAllocator<T2,A2> &) {
  return false;
}

} // namespace non_std

#endif

