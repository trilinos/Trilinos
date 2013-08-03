#ifndef KOKKOS_TEST_UNORDERED_MAP_INSERT_HPP
#define KOKKOS_TEST_UNORDERED_MAP_INSERT_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <impl/Kokkos_Timer.hpp>

#include <TestTimes.hpp>


#if defined(__CUDACC__)
//#define CUDA_ONE_THREAD_PER_BLOCK
#endif

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

    static const int max_grid_dim = 0x0000FFFF;
    static const int grid_shift = 16;

    test_insert_close(map_type map, uint32_t num_inserts, uint32_t num_duplicates)
      : m_map(map)
      , m_num_insert(num_inserts)
      , m_num_duplicates(num_duplicates)
    {
#if defined (CUDA_ONE_THREAD_PER_BLOCK)
      const dim3 block(1,1,1);

      int x = static_cast<int>(num_inserts & static_cast<uint32_t>(max_grid_dim));
      int y = (num_inserts >> grid_shift) ? (num_inserts >> grid_shift) : 1;

      dim3 grid(x, y, 1);
      Kokkos::Impl::CudaParallelLaunch< self_type >( *this , grid  , block , 0 );
#else
      Kokkos::parallel_for(num_inserts, *this);
#endif
    }


#if defined (CUDA_ONE_THREAD_PER_BLOCK)
    __device__ inline
    void operator()(void) const
    {
      typename device_type::size_type i = blockIdx.x + (blockIdx.y << grid_shift);
      if (i < m_num_insert) m_map.insert(i/m_num_duplicates, uint32_t(i));
    }
#endif

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

} // namespace Impl




template <typename Device>
void test_insert_close(  uint32_t num_nodes
                       , uint32_t num_inserts
                       , uint32_t num_duplicates
                       , map_test_times & test_times
                      )
{
  typedef Kokkos::unordered_map<uint32_t,uint32_t, Device> map_type;
  typedef Kokkos::unordered_map<const uint32_t,const uint32_t, Device> const_map_type;

  const uint32_t expected_inserts = (num_inserts + num_duplicates -1u) / num_duplicates;

  Kokkos::Impl::Timer timer;

  map_type map(num_nodes);
  Device::fence();

  test_times.construct += timer.seconds();

  timer.reset();
  ASSERT_NO_THROW(map.check_sanity());
  test_times.santity_check += timer.seconds();

  timer.reset();
  Impl::test_insert_close<map_type> test_insert_close(map, num_inserts, num_duplicates);
  Device::fence();
  test_times.insert += timer.seconds();

#if 1
  timer.reset();
  map.check_sanity();
  test_times.santity_check += timer.seconds();
#else
  try {
    map.check_sanity();
  }
  catch( std::runtime_error & err)
  {
    std::cout << err.what();
    ASSERT_TRUE(false);
  }
#endif

  const uint32_t map_size = map.size();
  const uint32_t num_failed_inserts = map.num_failed_inserts();

  ASSERT_EQ( num_failed_inserts, 0u);

  if (num_failed_inserts == 0u) {
    ASSERT_EQ(map_size, expected_inserts);

    timer.reset();
    uint32_t find_errors = 0;
    Impl::test_find<const_map_type> test_find(map, num_inserts, num_duplicates, find_errors);
    test_times.find += timer.seconds();
    ASSERT_EQ( find_errors, 0u);
  }
}

template <typename Device>
void test_insert_far(  uint32_t num_nodes
                       , uint32_t num_inserts
                       , uint32_t num_duplicates
                       , map_test_times & test_times
                      )
{
  typedef Kokkos::unordered_map<uint32_t,uint32_t, Device> map_type;
  typedef Kokkos::unordered_map<const uint32_t,const uint32_t, Device> const_map_type;

  const uint32_t expected_inserts = (num_inserts + num_duplicates -1u) / num_duplicates;

  Kokkos::Impl::Timer timer;

  map_type map(num_nodes);
  Device::fence();

  test_times.construct += timer.seconds();

  timer.reset();
  ASSERT_NO_THROW(map.check_sanity());
  test_times.santity_check += timer.seconds();

  timer.reset();
  Impl::test_insert_far<map_type> test_insert_far(map, num_inserts, num_duplicates);
  Device::fence();
  test_times.insert += timer.seconds();

#if 1
  timer.reset();
  map.check_sanity();
  test_times.santity_check += timer.seconds();
#else
  try {
    map.check_sanity();
  }
  catch( std::runtime_error & err)
  {
    std::cout << err.what();
    ASSERT_TRUE(false);
  }
#endif

  const uint32_t map_size = map.size();
  const uint32_t num_failed_inserts = map.num_failed_inserts();

  ASSERT_EQ( num_failed_inserts, 0u);

  if (num_failed_inserts == 0u) {
    ASSERT_EQ(map_size, expected_inserts);

    timer.reset();
    uint32_t find_errors = 0;
    Impl::test_find<const_map_type> test_find(map, num_inserts, num_duplicates, find_errors);
    test_times.find += timer.seconds();
    ASSERT_EQ( find_errors, 0u);
  }
}


} // namespace Test

#endif //KOKKOS_TEST_UNORDERED_MAP_INSERT_HPP
