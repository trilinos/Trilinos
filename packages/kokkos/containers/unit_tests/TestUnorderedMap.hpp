#ifndef KOKKOS_TEST_UNORDERED_MAP_HPP
#define KOKKOS_TEST_UNORDERED_MAP_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <impl/Kokkos_Timer.hpp>


namespace Test {

namespace Impl {

template <typename MapType, bool Near = false>
struct TestInsert
{
  typedef MapType map_type;
  typedef typename map_type::device_type device_type;
  typedef uint32_t value_type;

  map_type map;
  uint32_t inserts;
  uint32_t collisions;

  TestInsert( map_type arg_map, uint32_t arg_inserts, uint32_t arg_collisions)
    : map(arg_map)
    , inserts(arg_inserts)
    , collisions(arg_collisions)
  {}

  void apply( bool rehash_on_fail = true )
  {
    Kokkos::Impl::Timer wall_clock ;
    wall_clock.reset();

    uint32_t failed_count = 0;
    int loop_count = 0;
    do {
      ++loop_count;

      failed_count = 0;
      Kokkos::parallel_reduce(inserts, *this, failed_count);

      if (rehash_on_fail && failed_count > 0u) {
        const uint32_t new_capacity = map.capacity() + ((map.capacity()*3ull)/20u) + failed_count/collisions ;
        map.rehash( new_capacity );
        map.set_full_search();
      }
    } while (rehash_on_fail && failed_count > 0u);
  }


  KOKKOS_INLINE_FUNCTION
  void init( value_type & failed_count ) const { failed_count = 0; }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type & failed_count, const volatile value_type & count ) const
  { failed_count += count; }

  KOKKOS_INLINE_FUNCTION
  void operator()(uint32_t i, value_type & failed_count) const
  {
    const uint32_t key = Near ? i/collisions : i%(inserts/collisions);
    if (map.insert(key,i).failed()) ++failed_count;
  }

};

  template <typename MapType, bool Near>
  struct TestErase
  {
    typedef TestErase<MapType, Near> self_type;

    typedef MapType map_type;
    typedef typename MapType::device_type device_type;

    map_type m_map;
    uint32_t m_num_erase;
    uint32_t m_num_duplicates;

    TestErase(map_type map, uint32_t num_erases, uint32_t num_duplicates)
      : m_map(map)
      , m_num_erase(num_erases)
      , m_num_duplicates(num_duplicates)
    {
      Kokkos::parallel_for(num_erases, *this);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(typename device_type::size_type i) const
    {
      if (Near) {
        m_map.erase(i/m_num_duplicates);
      }
      else {
        m_map.erase(i%(m_num_erase/m_num_duplicates));
      }

    }
  };

  template <typename MapType>
  struct TestFind
  {
    typedef MapType map_type;
    typedef typename MapType::device_type device_type;
    typedef uint32_t value_type;

    map_type m_map;
    uint32_t m_num_insert;
    uint32_t m_num_duplicates;

    TestFind(map_type map, uint32_t num_inserts, uint32_t num_duplicates, value_type & errors)
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
void test_insert( uint32_t num_nodes , uint32_t num_inserts , uint32_t num_duplicates , bool near )
{
  typedef Kokkos::UnorderedMap<uint32_t,uint32_t, Device> map_type;
  typedef Kokkos::UnorderedMap<const uint32_t,const uint32_t, Device> const_map_type;

  const uint32_t expected_inserts = (num_inserts + num_duplicates -1u) / num_duplicates;

  map_type map(num_nodes);
  map.set_full_search();
  Device::fence();

  if (near) {
    Impl::TestInsert<map_type,true> test_insert(map, num_inserts, num_duplicates);
    test_insert.apply();
  } else
  {
    Impl::TestInsert<map_type,false> test_insert(map, num_inserts, num_duplicates);
    test_insert.apply();
  }
  Device::fence();

  const uint32_t map_size = map.size();

  ASSERT_FALSE( map.failed_insert());
  {
    EXPECT_EQ(map_size, expected_inserts);

    {
      uint32_t find_errors = 0;
      Impl::TestFind<const_map_type> test_find(map, num_inserts, num_duplicates, find_errors);
      EXPECT_EQ( find_errors, 0u);
    }

    map.begin_erase();
    Impl::TestErase<map_type,false> erase_far(map, num_inserts, num_duplicates);
    map.end_erase();
    EXPECT_EQ(map.size(), 0u);
  }
}

template <typename Device>
void test_failed_insert( uint32_t num_nodes)
{
  typedef Kokkos::UnorderedMap<uint32_t,uint32_t, Device> map_type;

  map_type map(num_nodes);
  Impl::TestInsert<map_type> test_insert(map, 2u*num_nodes, 1u);
  test_insert.apply(false /*don't rehash on fail*/);
  Device::fence();

  EXPECT_TRUE( map.failed_insert() );
}



template <typename Device>
void test_deep_copy( uint32_t num_nodes )
{
  typedef typename Device::host_mirror_device_type host_type ;

  typedef Kokkos::UnorderedMap<uint32_t,uint32_t, Device> map_type;
  typedef Kokkos::UnorderedMap<const uint32_t, const uint32_t, Device> const_map_type;

  typedef Kokkos::UnorderedMap<uint32_t, uint32_t, host_type> host_map_type;

  map_type map((num_nodes*5ull)/4u);
  map.set_full_search();
  Device::fence();

  {
    Impl::TestInsert<map_type> test_insert(map, num_nodes, 1);
    test_insert.apply();
    Device::fence();
    EXPECT_EQ( map.size(), num_nodes);
  }

  host_map_type hmap;
  Kokkos::deep_copy(hmap, map);

  map_type mmap;
  Kokkos::deep_copy(mmap, hmap);

  const_map_type cmap = mmap;

  EXPECT_EQ( cmap.size(), num_nodes);

  uint32_t find_errors = 0;
  Impl::TestFind<const_map_type> test_find(cmap, num_nodes/2u, 1, find_errors);
  EXPECT_EQ( find_errors, 0u);

}

} // namespace Test

#endif //KOKKOS_TEST_UNORDERED_MAP_HPP
