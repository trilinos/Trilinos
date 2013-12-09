#ifndef STK_UTIL_STK_UTIL_UTIL_PAGE_ALIGNED_ALLOCATOR_HPP
#define STK_UTIL_STK_UTIL_UTIL_PAGE_ALIGNED_ALLOCATOR_HPP

#include <stk_util/util/AllocatorMemoryUsage.hpp>

#include <cstdlib>
#include <limits>

#include <boost/type_traits/is_same.hpp>

namespace stk {

namespace detail {

struct page_aligned_allocator_impl
{
  static const size_t  half_page_size;

  static void * allocate(size_t num_bytes);
  static void deallocate(void * ptr, size_t num_bytes);

  static const int     m_mmap_flags;
  static const int     m_mmap_protection;
};

} // namespace detail

template <typename T, typename Tag = void>
class page_aligned_allocator
{
public:

  typedef Tag                         tag;
  typedef allocator_memory_usage<tag> memory_usage;

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
    typedef page_aligned_allocator<U,tag> other;
  };


  // constructors
  page_aligned_allocator() {}

  page_aligned_allocator(const page_aligned_allocator&) {}

  template <typename U>
  page_aligned_allocator (const page_aligned_allocator<U,tag>&) {}

  // destructor
  ~page_aligned_allocator() {}

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

    pointer ret;

    if (use_page_aligned_memory(size)) {
      ret = static_cast<pointer>(detail::page_aligned_allocator_impl::allocate(size));
    }
    else {
      ret = static_cast<pointer>(malloc(size));
    }
    return ret;
  }

  // deallocate storage p of deleted elements
  static void deallocate(pointer p, size_type num)
  {
    size_t size = num * sizeof(T);

    memory_usage::deallocate(size);

    if (use_page_aligned_memory(size)) {
      detail::page_aligned_allocator_impl::deallocate(p, size);
    }
    else {
      free(p);
    }
  }

  // initialize elements of allocated storage p with value value
  static void construct(pointer p, const T& value)
  {
    new(p)T(value);
  }

  // destroy elements of initialized storage p
  static void destroy(pointer p)
  {
    p->~T();
  }

private:

  static bool use_page_aligned_memory(size_type num_bytes)
  {
    return num_bytes >= detail::page_aligned_allocator_impl::half_page_size;
  }

};

// return that all specializations of the page_aligned_allocator with the same allocator and same tag are interchangeable
template <typename T1, typename T2, typename Tag1, typename Tag2>
inline bool operator==(const page_aligned_allocator<T1,Tag1>&, const page_aligned_allocator<T2,Tag2>&)
{ return boost::is_same<Tag1,Tag2>::value; }

template <typename T1, typename T2, typename Tag1, typename Tag2>
inline bool operator!=(const page_aligned_allocator<T1,Tag1>&, const page_aligned_allocator<T2,Tag2>&)
{ return !boost::is_same<Tag1,Tag2>::value; }

} // namespace stk

#endif //STK_UTIL_STK_UTIL_UTIL_PAGE_ALIGNED_ALLOCATOR_HPP
