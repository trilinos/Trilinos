/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_Pool_alloc_hpp_
#define _fei_Pool_alloc_hpp_

#include "fei_macros.hpp"
#include "fei_Pool.hpp"
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <new>
#include <stdexcept>
#include <iostream>

fei_Pool* get_fei_mem_pool(size_t n);

/** fei_Pool_alloc (along with its underlying memory-pool class
  fei_Pool) is taken almost verbatim from Stroustrup's book
  "The C++ Programming Language 3rd edition", pages 567-572.

  The difference between this class and the one in the book is that the length
  of allocated memory-chunks must always be the same. That is, the first time
  allocate is called, an internal memory-pool is created with the specified
  alloc-size. That memory-pool can only satisfy allocation requests of that
  same size from then on.
*/
template<typename T>
class fei_Pool_alloc {
 private:
  fei_Pool* mem; //pool of elements
  size_t n_;

 public:
  typedef T value_type;
  typedef std::size_t size_type;
  typedef T* pointer;
  typedef const T* const_pointer;
  typedef T& reference;
  typedef const T& const_reference;
  typedef std::ptrdiff_t difference_type;

  // Boilerplate allocator stuff
  template <typename U>
  struct rebind
  {
    typedef fei_Pool_alloc<U> other;
  };

  pointer address (reference value) const
  {
    return &value;
  }
  const_pointer address (const_reference value) const
  {
    return &value;
  }

  fei_Pool_alloc() throw();
  fei_Pool_alloc(const T&) throw();
  template<typename U> fei_Pool_alloc(const fei_Pool_alloc<U>&) throw()
   : mem(NULL),n_(0) {}

  ~fei_Pool_alloc() throw();

  pointer allocate(size_type n, const void* hint = NULL);
  void deallocate(pointer p, size_type n);

  template<typename U> void construct(U* p, const U& val)
  { new(p) U(val); }

  void construct(pointer p, const T& val)
  { new(p) T(val); }

  template<typename U> void destroy(U* p)
  { p->~U(); }

  void destroy(pointer p)
  { p->~T(); }

  size_type max_size() const throw() { return std::numeric_limits<size_type>::max(); }

};

template<typename T> fei_Pool_alloc<T>::fei_Pool_alloc() throw() : mem(NULL),n_(0) {}
template<typename T> fei_Pool_alloc<T>::fei_Pool_alloc(const T&) throw() : mem(NULL),n_(0) {}

template<typename T> fei_Pool_alloc<T>::~fei_Pool_alloc() throw() {}

template<typename T>
T* fei_Pool_alloc<T>::allocate(size_type n, const void*)
{
  if (n==0) return NULL;

  if (n_ == 0) {
    n_ = n;
    mem = get_fei_mem_pool(n_*sizeof(T));
  }

  if (n != n_) {
    std::cerr << "fei_Pool_alloc ERROR, allocate given bad length ("<<n
       <<"), must be " <<n_<<". throwing exception."<<std::endl;
    throw std::bad_alloc();
  }
  return static_cast<T*>(mem->alloc());
}

template<typename T>
void fei_Pool_alloc<T>::deallocate(pointer p, size_type n)
{
  if (p == NULL || n == 0) return;

  if (n == n_) {
    mem->free(p);
    return;
  }

  std::cerr << "fei_Pool_alloc ERROR, deallocate given bad length ("<<n
    <<"), must be " <<n_<<". aborting."<<std::endl;
  std::abort();
}

template<typename T>
inline bool operator==(const fei_Pool_alloc<T>&,
                       const fei_Pool_alloc<T>&) throw()
{ return true; }
template<typename T>
inline bool operator!=(const fei_Pool_alloc<T>&,
                       const fei_Pool_alloc<T>&) throw()
{ return false; }

#endif

