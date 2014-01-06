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
  struct test_erase_close
  {
    typedef test_erase_close<MapType> self_type;

    typedef MapType map_type;
    typedef typename MapType::device_type device_type;

    map_type m_map;
    uint32_t m_num_erase;
    uint32_t m_num_duplicates;

    test_erase_close(map_type map, uint32_t num_erases, uint32_t num_duplicates)
      : m_map(map)
      , m_num_erase(num_erases)
      , m_num_duplicates(num_duplicates)
    {
      Kokkos::parallel_for(num_erases, *this);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(typename device_type::size_type i) const
    {
      m_map.erase(i/m_num_duplicates);
    }
  };

  template <typename MapType>
  struct test_erase_far
  {
    typedef MapType map_type;
    typedef typename MapType::device_type device_type;

    map_type m_map;
    uint32_t m_num_erase;
    uint32_t m_num_duplicates;

    test_erase_far(map_type map, uint32_t num_erases, uint32_t num_duplicates)
      : m_map(map)
      , m_num_erase(num_erases)
      , m_num_duplicates(num_duplicates)
    {
      Kokkos::parallel_for(num_erases, *this);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(typename device_type::size_type i) const
    {
      m_map.erase(i%(m_num_erase/m_num_duplicates));
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

      const bool exists = m_map.exists(i);

      if (expect_to_find_i && !exists)  ++errors;
      if (!expect_to_find_i && exists)  ++errors;
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

  Impl::test_insert_close<map_type> test_insert_close(map, num_inserts, num_duplicates);
  Device::fence();

  const uint32_t map_size = map.size();
  const bool failed_inserts = map.failed_inserts();

  ASSERT_FALSE( failed_inserts );

  if (!failed_inserts) {
    ASSERT_EQ(map_size, expected_inserts);

    {
      uint32_t find_errors = 0;
      Impl::test_find<const_map_type> test_find(map, num_inserts, num_duplicates, find_errors);
      ASSERT_EQ( find_errors, 0u);
    }

    map.begin_erase();
    Impl::test_erase_close<map_type> erase_close(map, num_inserts, num_duplicates);
    map.end_erase();
    ASSERT_EQ(map.size(), 0u);

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

  Impl::test_insert_far<map_type> test_insert_far(map, num_inserts, num_duplicates);
  Device::fence();

  const uint32_t map_size = map.size();
  const bool failed_inserts = map.failed_inserts();

  ASSERT_FALSE( failed_inserts );

  if (!failed_inserts) {
    ASSERT_EQ(map_size, expected_inserts);

    {
      uint32_t find_errors = 0;
      Impl::test_find<const_map_type> test_find(map, num_inserts, num_duplicates, find_errors);
      ASSERT_EQ( find_errors, 0u);
    }

    map.begin_erase();
    Impl::test_erase_far<map_type> erase_far(map, num_inserts, num_duplicates);
    map.end_erase();
    ASSERT_EQ(map.size(), 0u);
  }
}

template <typename Device>
void test_failed_insert( uint32_t num_nodes)
{
  typedef Kokkos::UnorderedMap<uint32_t,uint32_t, Device> map_type;
  typedef Kokkos::UnorderedMap<const uint32_t,const uint32_t, Device> const_map_type;

  map_type map(num_nodes);
  Device::fence();

  Impl::test_insert_far<map_type> test_insert_far(map, 2u*num_nodes, 1u);
  Device::fence();
  //ASSERT_TRUE( map.failed_inserts() );
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

  {
    Impl::test_insert_far<map_type> test_insert_far(map, num_nodes, 1);
    Device::fence();
    ASSERT_EQ( map.size(), num_nodes);
  }

  host_map_type hmap;
  Kokkos::deep_copy(hmap, map);

  map_type mmap;
  Kokkos::deep_copy(mmap, hmap);

  const_map_type cmap = mmap;

  ASSERT_EQ( cmap.size(), num_nodes);

  uint32_t find_errors = 0;
  Impl::test_find<const_map_type> test_find(cmap, num_nodes/2u, 1, find_errors);
  ASSERT_EQ( find_errors, 0u);

}

} // namespace Test

#endif //KOKKOS_TEST_UNORDERED_MAP_HPP
