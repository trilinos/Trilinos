// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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

  AlignedAllocator() noexcept { }

  AlignedAllocator(const AlignedAllocator& other) noexcept
    : std::allocator<T>(other) { }

  template <class U>
  AlignedAllocator(const AlignedAllocator<U,Alignment>&) noexcept { }

  ~AlignedAllocator() noexcept { }

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

