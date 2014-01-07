#ifndef KOKKOS_TEST_UNORDERED_MAP_HPP
#define KOKKOS_TEST_UNORDERED_MAP_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <impl/Kokkos_Timer.hpp>


namespace Test {

namespace Impl {

  template <typename MapType>
  struct test_insert_close
  {
    typedef test_insert_close<MapType> self_type;

    typedef MapType map_type;
    typedef typename MapType::device_type device_type;

    map_type m_map;
    uint32_t m_num_insert;
    uint32_t m_num_duplicates;

    test_insert_close(map_type map, uint32_t num_inserts, uint32_t num_duplicates)
      : m_map(map)
      , m_num_insert(num_inserts)
      , m_num_duplicates(num_duplicates)
    {
      Kokkos::parallel_for(num_inserts, *this);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(typename device_type::size_type i) const
    {
      m_map.insert(i/m_num_duplicates, uint32_t(i));
    }
  };

  template <typename MapType>
  struct test_insert_far
  {
    typedef MapType map_type;
    typedef typename MapType::device_type device_type;

    map_type m_map;
    uint32_t m_num_insert;
    uint32_t m_num_duplicates;

    test_insert_far(map_type map, uint32_t num_inserts, uint32_t num_duplicates)
      : m_map(map)
      , m_num_insert(num_inserts)
      , m_num_duplicates(num_duplicates)
    {
      Kokkos::parallel_for(num_inserts, *this);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(typename device_type::size_type i) const
    {
      m_map.insert(i%(m_num_insert/m_num_duplicates), uint32_t(i));
    }
  };

  template <typename MapType>
  struct test_find
  {
    typedef MapType map_type;
    typedef typename MapType::device_type device_type;
    typedef uint32_t value_type;

    map_type m_map;
    uint32_t m_num_insert;
    uint32_t m_num_duplicates;

    test_find(map_type map, uint32_t num_inserts, uint32_t num_duplicates, value_type & errors)
      : m_map(map)
      , m_num_insert(num_inserts)
      , m_num_duplicates(num_duplicates)
    {
      Kokkos::parallel_reduce(num_inserts, *this, errors);
      device_type::fence();
    }

    KOKKOS_INLINE_FUNCTION
    static void init( value_type & dst)
    {
      dst = 0;
    }

    KOKKOS_INLINE_FUNCTION
    static void join( volatile value_type & dst, const volatile value_type & src)
    { dst += src; }

    KOKKOS_INLINE_FUNCTION
    void operator()(typename device_type::size_type i, value_type & errors) const
    {
      const uint32_t max_i = (m_num_insert + m_num_duplicates -1u) / m_num_duplicates;
      const bool expect_to_find_i = (i < max_i);

      typename map_type::const_pointer ptr = m_map.find(i);

      if (expect_to_find_i && !ptr)  ++errors;
      if (!expect_to_find_i && ptr)  ++errors;
    }
  };

  template <typename MapType>
  struct test_insert_mark_pending_delete
  {
    typedef MapType map_type;
    typedef typename MapType::device_type device_type;

    map_type m_map;
    uint32_t m_num_insert;
    uint32_t m_num_duplicates;

    test_insert_mark_pending_delete(map_type map, uint32_t num_inserts, uint32_t num_duplicates)
      : m_map(map)
      , m_num_insert(num_inserts)
      , m_num_duplicates(num_duplicates)
    {
      Kokkos::parallel_for(num_inserts, *this);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(typename device_type::size_type i) const
    {
      m_map.insert(i%(m_num_insert/m_num_duplicates), uint32_t(i));
      m_map.mark_pending_delete(i/m_num_duplicates);
    }
  };

  template <typename MapType>
  struct mark_pending_delete_functor
  {
    typedef MapType map_type;
    typedef typename MapType::device_type device_type;

    map_type m_map;

    mark_pending_delete_functor(map_type map)
      : m_map(map)
    {
      Kokkos::parallel_for(map.capacity(), *this);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(typename device_type::size_type i) const
    {
      typename map_type::const_pointer ptr = m_map.get_value(i);
      if (ptr != NULL)
        m_map.mark_pending_delete(ptr->first);
    }
  };


} // namespace Impl




template <typename Device>
void test_insert_close(  uint32_t num_nodes
                       , uint32_t num_inserts
                       , uint32_t num_duplicates
                      )
{
  typedef Kokkos::UnorderedMap<uint32_t,uint32_t, Device> map_type;
  typedef Kokkos::UnorderedMap<const uint32_t,const uint32_t, Device> const_map_type;

  const uint32_t expected_inserts = (num_inserts + num_duplicates -1u) / num_duplicates;

  map_type map(num_nodes);
  Device::fence();
  ASSERT_NO_THROW(map.check_sanity());

  Impl::test_insert_close<map_type> test_insert_close(map, num_inserts, num_duplicates);
  Device::fence();
  map.check_sanity();

  const uint32_t map_size = map.size();
  const bool failed_inserts = map.failed_inserts();

  ASSERT_FALSE( failed_inserts );

  if (!failed_inserts) {
    ASSERT_EQ(map_size, expected_inserts);

    uint32_t find_errors = 0;
    Impl::test_find<const_map_type> test_find(map, num_inserts, num_duplicates, find_errors);
    ASSERT_EQ( find_errors, 0u);
  }
}

template <typename Device>
void test_insert_far(  uint32_t num_nodes
                       , uint32_t num_inserts
                       , uint32_t num_duplicates
                      )
{
  typedef Kokkos::UnorderedMap<uint32_t,uint32_t, Device> map_type;
  typedef Kokkos::UnorderedMap<const uint32_t,const uint32_t, Device> const_map_type;

  const uint32_t expected_inserts = (num_inserts + num_duplicates -1u) / num_duplicates;

  map_type map(num_nodes);
  Device::fence();
  ASSERT_NO_THROW(map.check_sanity());

  Impl::test_insert_far<map_type> test_insert_far(map, num_inserts, num_duplicates);
  Device::fence();
  map.check_sanity();

  const uint32_t map_size = map.size();
  const bool failed_inserts = map.failed_inserts();

  ASSERT_FALSE( failed_inserts );

  if (!failed_inserts) {
    ASSERT_EQ(map_size, expected_inserts);

    uint32_t find_errors = 0;
    Impl::test_find<const_map_type> test_find(map, num_inserts, num_duplicates, find_errors);
    ASSERT_EQ( find_errors, 0u);
  }
}

template <typename Device>
void test_failed_insert( uint32_t num_nodes)
{
  typedef Kokkos::UnorderedMap<uint32_t,uint32_t, Device> map_type;
  typedef Kokkos::UnorderedMap<const uint32_t,const uint32_t, Device> const_map_type;

  map_type map(num_nodes);
  Device::fence();
  ASSERT_NO_THROW(map.check_sanity());

  Impl::test_insert_far<map_type> test_insert_far(map, 2u*num_nodes, 1u);
  Device::fence();
  ASSERT_TRUE( map.failed_inserts() );
  ASSERT_THROW(map.check_sanity(), std::runtime_error);
}

template <typename Device>
void test_insert_mark_pending_delete(  uint32_t num_nodes
                                     , uint32_t num_inserts
                                     , uint32_t num_duplicates
                                    )
{
  typedef Kokkos::UnorderedMap<uint32_t,uint32_t, Device> map_type;
  typedef Kokkos::UnorderedMap<const uint32_t,const uint32_t, Device> const_map_type;

  const uint32_t expected_inserts = (num_inserts + num_duplicates -1u) / num_duplicates;

  map_type map(num_nodes);
  Device::fence();
  ASSERT_NO_THROW(map.check_sanity());

  Impl::test_insert_mark_pending_delete<map_type> test_insert_far(map, num_inserts, num_duplicates);
  Device::fence();
  map.check_sanity();

  const uint32_t map_size = map.size();
  const uint32_t pending_delete = map.pending_delete();
  const uint32_t failed_inserts = map.failed_inserts();

  ASSERT_FALSE( failed_inserts );

  if (failed_inserts == 0u) {
    ASSERT_EQ(map_size, 0u);
    ASSERT_EQ(pending_delete, expected_inserts);

    // keys not deleted so should be able to find them
    {
      uint32_t find_errors = 0;
      Impl::test_find<const_map_type> test_find(map, num_inserts, num_duplicates, find_errors);
      ASSERT_EQ( find_errors, 0u);
    }

    // remove the keys marked pending_delete
    map.remove_pending_delete();
    ASSERT_EQ(map.pending_delete(), 0u);
    ASSERT_EQ(map.unused(), map.capacity());

    // keys deleted so should NOT be able to find them
    {
      uint32_t find_errors = 0;
      Impl::test_find<const_map_type> test_find(map, num_inserts, num_duplicates, find_errors);
      ASSERT_EQ( find_errors, expected_inserts);
    }
  }
}

template <typename Device>
void test_assignement_operators( uint32_t num_nodes )
{
  typedef Kokkos::UnorderedMap<uint32_t,uint32_t, Device> map_type;
  typedef Kokkos::UnorderedMap<const uint32_t, uint32_t, Device> non_insertable_map_type;
  typedef Kokkos::UnorderedMap<const uint32_t, const uint32_t, Device> const_map_type;

  map_type map(num_nodes);
  Device::fence();
  ASSERT_NO_THROW(map.check_sanity());

  {
    Impl::test_insert_far<map_type> test_insert_far(map, num_nodes, 1);
    Device::fence();
    map.check_sanity();
    ASSERT_EQ( map.size(), num_nodes);
  }

  non_insertable_map_type nmap = map;

  {
    Impl::mark_pending_delete_functor<non_insertable_map_type> mark_delete(nmap);
    Device::fence();
    map.check_sanity();
    ASSERT_EQ( map.size(), 0u);
    ASSERT_EQ( map.pending_delete(), num_nodes);
  }

  const_map_type cmap = nmap;

  cmap.remove_pending_delete();
  ASSERT_EQ( map.size(), 0u);
  ASSERT_EQ( map.pending_delete(), 0u);
}

template <typename Device>
void test_deep_copy( uint32_t num_nodes )
{
  typedef typename Device::host_mirror_device_type host_type ;

  typedef Kokkos::UnorderedMap<uint32_t,uint32_t, Device> map_type;
  typedef Kokkos::UnorderedMap<const uint32_t, uint32_t, Device> non_insertable_map_type;
  typedef Kokkos::UnorderedMap<const uint32_t, const uint32_t, Device> const_map_type;

  typedef Kokkos::UnorderedMap<uint32_t, uint32_t, host_type> host_map_type;
  typedef Kokkos::UnorderedMap<const uint32_t, const uint32_t, host_type> const_host_map_type;

  map_type map(num_nodes);
  Device::fence();
  ASSERT_NO_THROW(map.check_sanity());

  {
    Impl::test_insert_far<map_type> test_insert_far(map, num_nodes, 2);
    Device::fence();
    map.check_sanity();
    ASSERT_EQ( map.size(), num_nodes/2u);
  }

  host_map_type hmap;
  Kokkos::deep_copy(hmap, map);

  map_type mmap;
  Kokkos::deep_copy(mmap, hmap);

  const_map_type cmap = mmap;

  cmap.shrink_to_fit();

  ASSERT_EQ( cmap.size(), num_nodes/2u);

  uint32_t find_errors = 0;
  Impl::test_find<const_map_type> test_find(cmap, num_nodes/2u, 1, find_errors);
  ASSERT_EQ( find_errors, 0u);

}

} // namespace Test

#endif //KOKKOS_TEST_UNORDERED_MAP_HPP
