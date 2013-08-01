#ifndef KOKKOS_TEST_UNORDERED_MAP_INSERT_HPP
#define KOKKOS_TEST_UNORDERED_MAP_INSERT_HPP

#include <gtest/gtest.h>
#include <iostream>

namespace Test {

namespace Impl {

  template <typename MapType>
  struct test_insert
  {
    typedef MapType map_type;
    typedef typename MapType::device_type device_type;

    map_type m_map;
    uint32_t m_num_insert;
    uint32_t m_num_duplicates;

    test_insert(map_type map, uint32_t num_inserts, uint32_t num_duplicates)
      : m_map(map)
      , m_num_insert(num_inserts)
      , m_num_duplicates(num_duplicates)
    {
      Kokkos::parallel_for(num_inserts, *this);
      device_type::fence();
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(typename device_type::size_type i) const
    {
      m_map.insert(i/m_num_duplicates);
    }
  };

  template <typename MapType>
  struct test_insert_2
  {
    typedef MapType map_type;
    typedef typename MapType::device_type device_type;

    map_type m_map;
    uint32_t m_num_insert;
    uint32_t m_num_duplicates;

    test_insert_2(map_type map, uint32_t num_inserts, uint32_t num_duplicates)
      : m_map(map)
      , m_num_insert(num_inserts)
      , m_num_duplicates(num_duplicates)
    {
      Kokkos::parallel_for(num_inserts, *this);
      device_type::fence();
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(typename device_type::size_type i) const
    {
      m_map.insert(i%(m_num_insert/m_num_duplicates));
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

} // namespace Imp

template <typename Device>
struct test_unordered_map_insert
{
  typedef Device device_type;

  typedef Kokkos::unordered_map<uint32_t, void, device_type> map_type;
  typedef Kokkos::unordered_map<const uint32_t, void, device_type> const_map_type;

  test_unordered_map_insert(  uint32_t num_nodes
                            , uint32_t num_inserts
                            , uint32_t num_duplicates
                           )
  {
    const uint32_t expected_inserts = (num_inserts + num_duplicates -1u) / num_duplicates;

    {
      map_type map(num_nodes);

      EXPECT_EQ(map.capacity(), num_nodes);

      Impl::test_insert<map_type> test_insert(map, num_inserts, num_duplicates);

      const uint32_t map_size = map.size();
      const uint32_t num_failed_inserts = map.num_failed_inserts();


      if (num_failed_inserts == 0u) {
        EXPECT_EQ(map_size, expected_inserts);

        if (map_size != expected_inserts) {
          uint32_t num_duplicates_keys = map.count_duplicate_keys();
          EXPECT_EQ(num_duplicates_keys, 0u);
        }
        // no keys have been marked for deletions so capacity == size + unused
        //EXPECT_EQ(map.capacity(), map_size + map.unused());

        uint32_t find_errors = 0;
        const_map_type const_map(map);
        Impl::test_find<const_map_type> test_find(const_map, num_inserts, num_duplicates, find_errors);
        //Impl::test_find<map_type> test_find(map, num_inserts, num_duplicates, find_errors);
        EXPECT_EQ( find_errors, 0u);
      }
      else { // num_failed_inserts > 0u
        EXPECT_LT( map_size, expected_inserts);
      }

      // insert should not fail
      if ( (9u*num_nodes)/10u >= num_inserts) {
        EXPECT_EQ(map.num_failed_inserts(), 0u);
      }

      if (num_nodes < num_inserts) {
        EXPECT_GT(map.num_failed_inserts(), 0u);
      }
    }

    {
      map_type map(num_nodes);

      EXPECT_EQ(map.capacity(), num_nodes);

      Impl::test_insert_2<map_type> test_insert(map, num_inserts, num_duplicates);

      const uint32_t map_size = map.size();
      const uint32_t num_failed_inserts = map.num_failed_inserts();


      if (num_failed_inserts == 0u) {
        EXPECT_EQ(map_size, expected_inserts);

        if (map_size != expected_inserts) {
          uint32_t num_duplicates_keys = map.count_duplicate_keys();
          EXPECT_EQ(num_duplicates_keys, 0u);
        }
        // no keys have been marked for deletions so capacity == size + unused
        EXPECT_EQ(map.capacity(), map_size + map.unused());

        uint32_t find_errors = 0;
        const_map_type const_map(map);
        Impl::test_find<const_map_type> test_find(const_map, num_inserts, num_duplicates, find_errors);
        //Impl::test_find<map_type> test_find(map, num_inserts, num_duplicates, find_errors);
        EXPECT_EQ( find_errors, 0u);
      }
      else { // num_failed_inserts > 0u
        EXPECT_LT( map_size, expected_inserts);
      }

      // insert should not fail
      if ( (9u*num_nodes)/10u >= num_inserts) {
        EXPECT_EQ(map.num_failed_inserts(), 0u);
      }

      if (num_nodes < num_inserts) {
        EXPECT_GT(map.num_failed_inserts(), 0u);
      }

    }
  }


};


} // namespace Test

#endif //KOKKOS_TEST_UNORDERED_MAP_INSERT_HPP
