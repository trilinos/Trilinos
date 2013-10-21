#ifndef STK_UTIL_STK_UTIL_UTIL_CACHE_ALIGNED_ALLOCATOR_HPP
#define STK_UTIL_STK_UTIL_UTIL_CACHE_ALIGNED_ALLOCATOR_HPP

#include <stk_util/util/AllocatorMemoryUsage.hpp>
#include <cstdlib>
#include <limits>
#include <boost/type_traits/is_same.hpp>
#include <boost/static_assert.hpp>

namespace stk {

template <typename T, typename Tag = void, size_t CacheSize = 64 >
class cache_aligned_allocator
{
public:
  BOOST_STATIC_ASSERT(( CacheSize != 0u && !( CacheSize & (CacheSize-1u)) ));

  typedef Tag                         tag;
  typedef stk::allocator_memory_usage<tag> memory_usage;

  // type definitions
  typedef T              value_type;
  typedef T*             pointer;
  typedef const T*       const_pointer;
  typedef T&             reference;
  typedef const T&       const_reference;
  typedef std::size_t    size_type;
  typedef std::ptrdiff_t difference_type;

  // rebind allocator to type U
  template <typename U>
  struct rebind
  {
    typedef cache_aligned_allocator<U,tag> other;
  };

  // constructors
  cache_aligned_allocator() {}

  cache_aligned_allocator(const cache_aligned_allocator&) {}

  template <typename U>
  cache_aligned_allocator (const cache_aligned_allocator<U,tag>&) {}

  // destructor
  ~cache_aligned_allocator() {}

  // return address of values
  static pointer       address(      reference value) { return &value; }
  static const_pointer address(const_reference value) { return &value; }


  // return maximum number of elements that can be allocated
  static size_type max_size()
  {
    return std::numeric_limits<std::size_t>::max() / sizeof(T);
  }

  // allocate but don't initialize num elements of type T
  static pointer allocate(size_type num, const void* = 0)
  {
    size_t size = num * sizeof(T);

    memory_usage::allocate(size);

    pointer ptr = NULL;
#if defined( __INTEL_COMPILER )
    ptr = (pointer)_mm_malloc(size, CacheSize);
#else
    posix_memalign( (void**)&ptr, CacheSize, size );
#endif

    return ptr;
  }

  // deallocate storage p of deleted elements
  static void deallocate(pointer p, size_type num)
  {
    memory_usage::deallocate(num * sizeof(T));
#if defined( __INTEL_COMPILER )
    _mm_free(p);
#else
    free(p);
#endif
  }

  // initialize elements of allocated storage p with value value
  static void construct(pointer p, const T& value)
  {
    new((void*)p)T(value);
  }

  // destroy elements of initialized storage p
  static void destroy(pointer p)
  {
    p->~T();
  }
};

// return that all specializations of the cache_aligned_allocator with the same allocator and same tag are interchangeable
template <typename T1, typename T2, typename Tag1, typename Tag2>
inline bool operator==(const cache_aligned_allocator<T1,Tag1>&, const cache_aligned_allocator<T2,Tag2>&)
{ return boost::is_same<Tag1,Tag2>::value; }

template <typename T1, typename T2, typename Tag1, typename Tag2>
inline bool operator!=(const cache_aligned_allocator<T1,Tag1>&, const cache_aligned_allocator<T2,Tag2>&)
{ return !boost::is_same<Tag1,Tag2>::value; }

} // namespace stk

#endif // STK_UTIL_STK_UTIL_UTIL_CACHE_ALIGNED_ALLOCATOR_HPP

