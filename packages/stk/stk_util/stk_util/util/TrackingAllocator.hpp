#ifndef STK_UTIL_STK_UTIL_UTIL_TRACKING_ALLOCATOR_HPP
#define STK_UTIL_STK_UTIL_UTIL_TRACKING_ALLOCATOR_HPP

#include <stk_util/util/AllocatorMemoryUsage.hpp>

#include <cstdlib>
#include <limits>

#include <boost/type_traits/is_same.hpp>



namespace stk {


template <typename T, typename Tag = void>
class tracking_allocator
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
    typedef tracking_allocator<U,tag> other;
  };

  // constructors
  tracking_allocator() {}

  tracking_allocator(const tracking_allocator&) {}

  template <typename U>
  tracking_allocator (const tracking_allocator<U,tag>&) {}

  // destructor
  ~tracking_allocator() {}

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

    return (pointer) malloc(size);
  }

  // deallocate storage p of deleted elements
  static void deallocate(pointer p, size_type num)
  {
    memory_usage::deallocate(num * sizeof(T));
    free(p);
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

// return that all specializations of the tracking_allocator with the same allocator and same tag are interchangeable
template <typename T1, typename T2, typename Tag1, typename Tag2>
inline bool operator==(const tracking_allocator<T1,Tag1>&, const tracking_allocator<T2,Tag2>&)
{ return boost::is_same<Tag1,Tag2>::value; }

template <typename T1, typename T2, typename Tag1, typename Tag2>
inline bool operator!=(const tracking_allocator<T1,Tag1>&, const tracking_allocator<T2,Tag2>&)
{ return !boost::is_same<Tag1,Tag2>::value; }


} // namespace stk

#endif //STK_UTIL_STK_UTIL_UTIL_TRACKING_ALLOCATOR_HPP
