#ifndef STK_UTIL_STK_UTIL_UTIL_TRACKING_ALLOCATOR_HPP
#define STK_UTIL_STK_UTIL_UTIL_TRACKING_ALLOCATOR_HPP

#include <limits>
#include <boost/type_traits/is_same.hpp>

namespace stk {

// Based off of code example in "The C++ Standard Library - A Tutorial and Reference"

namespace detail {

template <typename Tag = void>
struct memory_usage
{
  static size_t  peak_memory;
  static size_t  current_memory;
  static size_t  num_allocations;
};

template <typename Tag>
size_t memory_usage<Tag>::peak_memory = 0;

template <typename Tag>
size_t memory_usage<Tag>::current_memory = 0;

template <typename Tag>
size_t memory_usage<Tag>::num_allocations = 0;


typedef memory_usage<void> default_memory_usage;

} // namespace details


template <typename T, typename Tag = void>
class tracking_allocator
{
public:

  typedef Tag tag_type;
  typedef detail::memory_usage<Tag> memory_usage;

  static size_t current_memory()  { return memory_usage::current_memory; }
  static size_t peak_memory()     { return memory_usage::peak_memory; }
  static size_t num_allocations() { return memory_usage::num_allocations; }

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
    typedef tracking_allocator<U,Tag> other;
  };

  // return address of values
  pointer       address(      reference value) const { return &value; }
  const_pointer address(const_reference value) const { return &value; }

  // constructors
  tracking_allocator() {}

  tracking_allocator(const tracking_allocator&) {}

  template <typename U>
  tracking_allocator (const tracking_allocator<U,Tag>&) {}

  // destructor
  ~tracking_allocator() {}

  // return maximum number of elements that can be allocated
  size_type max_size() const
  {
    return std::numeric_limits<std::size_t>::max() / sizeof(T);
  }

  // allocate but don't initialize num elements of type T
  pointer allocate(size_type num, const void* = 0)
  {
    size_t size = num * sizeof(T);

    ++memory_usage::num_allocations;
    memory_usage::current_memory += size;
    memory_usage::peak_memory =   (memory_usage::peak_memory < memory_usage::current_memory)
                                ? memory_usage::current_memory
                                : memory_usage::peak_memory;

    pointer ret = (pointer)(::operator new(size));

    return ret;
  }

  // initialize elements of allocated storage p with value value
  void construct(pointer p, const T& value)
  {
    new((void*)p)T(value);
  }

  // destroy elements of initialized storage p
  void destroy(pointer p)
  {
    p->~T();
  }

  // deallocate storage p of deleted elements
  void deallocate(pointer p, size_type num)
  {
    size_t size = num * sizeof(T);
    memory_usage::current_memory -= size;

    ::operator delete((void*)p);
  }

};



// return that all specializations of this allocator with same tag are interchangeable
template <typename T1, typename T2, typename Tag1, typename Tag2>
bool operator==(const tracking_allocator<T1,Tag1>&, const tracking_allocator<T2,Tag2>&)
{ return boost::is_same<Tag1,Tag2>::value; }

template <typename T1, typename T2, typename Tag1, typename Tag2>
bool operator!=(const tracking_allocator<T1,Tag1>&, const tracking_allocator<T2,Tag2>&)
{ return !boost::is_same<Tag1,Tag2>::value; }

} // namespace stk

#endif //STK_UTIL_STK_UTIL_UTIL_TRACKING_ALLOCATOR_HPP
