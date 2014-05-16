
#include <stddef.h>                     // for size_t
#include <algorithm>                    // for swap, sort, equal
#include <boost/timer.hpp>              // for timer
#include <cstdlib>                      // for rand, srand, RAND_MAX
#include <deque>                        // for _Deque_iterator, operator-, etc
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <gtest/gtest.h>
#include <stk_util/util/BlockVector.hpp>  // for block_vector<>::iterator, etc
#include <stk_util/util/RadixSort.hpp>  // for radix_sort_unsigned
#include <stk_util/util/RadixSort2.hpp>  // for radix_sort
#include <vector>                       // for vector, swap
#include "gtest/gtest.h"                // for AssertHelper, ASSERT_EQ, etc



TEST( radix_sort, sort_unsigned )
{
  size_t size = 1000;
  std::vector<unsigned> a(size);

  srand(0);
  for (size_t i=0; i<size; ++i) a[i] = std::rand();

  std::vector<unsigned> b(a);

  std::sort(a.begin(), a.end());

  stk::radix_sort(b);

  ASSERT_EQ( a.size(), b.size());
  ASSERT_TRUE( std::equal(a.begin(), a.end(), b.begin()));
}

TEST( radix_sort, sort_signed )
{
  size_t size = 1000;
  std::vector<int> a(size);

  srand(0);
  for (size_t i=0; i<size; ++i) a[i] = std::rand() - (RAND_MAX >> 1);

  std::vector<int> b(a);

  std::sort(a.begin(), a.end());

  stk::radix_sort(b);

  ASSERT_EQ( a.size(), b.size());
  ASSERT_TRUE( std::equal(a.begin(), a.end(), b.begin()));
}

TEST( radix_sort, sort_float )
{
  size_t size = 1000;
  std::vector<float> a(size);

  srand(0);
  for (size_t i=0; i<size; ++i)
    a[i] = std::rand() / float(RAND_MAX);

  std::vector<float> b(a);

  std::sort(a.begin(), a.end());

  stk::radix_sort(b);

  ASSERT_EQ( a.size(), b.size());
  ASSERT_TRUE( std::equal(a.begin(), a.end(), b.begin()));
}

TEST( radix_sort, sort_signed_double )
{
  size_t size = 1000;
  std::vector<double> a(size);

  srand(0);
  for (size_t i=0; i<size; ++i)
    a[i] = std::rand() / double(RAND_MAX) -0.5;

  std::vector<double> b(a);

  std::sort(a.begin(), a.end());

  stk::radix_sort(b);

  ASSERT_EQ( a.size(), b.size());
  ASSERT_TRUE( std::equal(a.begin(), a.end(), b.begin()));
}


namespace {

struct RadixTag {};
struct SortTag {};

} //namespace



TEST( radix_sort, DISABLED_performance_vs_radix_sort_unsigned_vector )
{
  typedef std::vector<size_t> Vector;

  size_t size = 1ull << 24;
  Vector a(size);
  for (size_t i=0; i<size; ++i) a[i] = std::rand();

  Vector b(a);

  boost::timer timer;

  timer.restart();
  stk::util::radix_sort_unsigned( &*a.begin(), a.size());
  double sort_time = timer.elapsed();

  timer.restart();
  stk::radix_sort(b);
  double radix_time = timer.elapsed();

  std::cout.precision(3);

  std::cout << "Compare stk::util::radix_sort_unsigned to stk::radix_sort => vector<size_t> " << std::endl;
  std::cout << "Size = 2^24" << std::endl;

  std::cout << "  stk::util::radix_sort_unsigned : " << sort_time << " secs " << std::endl;
  std::cout << "                 stk::radix_sort : " << radix_time << " secs, " << sort_time/radix_time << "x speedup" << std::endl;
}

TEST( radix_sort, DISABLED_performance_vs_std_sort_vector )
{
  typedef std::vector<size_t> Vector;

  size_t size = 1ull << 24;
  Vector a(size);
  for (size_t i=0; i<size; ++i) a[i] = std::rand();

  Vector b(a);

  boost::timer timer;

  timer.restart();
  std::sort(a.begin(), a.end());
  double sort_time = timer.elapsed();

  timer.restart();
  stk::radix_sort(b);
  double radix_time = timer.elapsed();

  std::cout << "Compare std::sort to stk::radix_sort => vector<size_t>" << std::endl;
  std::cout << "Size = 2^24" << std::endl;

  std::cout << "        std::sort : " << sort_time << " secs " << std::endl;
  std::cout << "  stk::radix_sort : " << radix_time << " secs, " << sort_time/radix_time  << "x speedup" << std::endl;

}

TEST( radix_sort, DISABLED_performance_vs_std_sort_block_vector )
{
  typedef stk::block_vector<size_t> Vector;

  size_t size = 1ull << 24;
  Vector a(size);
  for (size_t i=0; i<size; ++i) a[i] = std::rand();

  Vector b(a);

  boost::timer timer;

  timer.restart();
  std::sort(a.begin(), a.end());
  double sort_time = timer.elapsed();

  timer.restart();
  stk::radix_sort(b);
  double radix_time = timer.elapsed();

  std::cout << "Compare std::sort to stk::radix_sort => block_vector<size_t>" << std::endl;
  std::cout << "Size = 2^24" << std::endl;

  std::cout << "        std::sort : " << sort_time << " secs " << std::endl;
  std::cout << "  stk::radix_sort : " << radix_time << " secs, " << sort_time/radix_time << "x speedup"  << std::endl;

}

TEST( radix_sort, DISABLED_performance_vs_std_sort_deque )
{
  typedef std::deque<size_t> Vector;

  size_t size = 1ull << 24;
  Vector a(size);
  for (size_t i=0; i<size; ++i) a[i] = std::rand();

  Vector b(a);

  boost::timer timer;

  timer.restart();
  std::sort(a.begin(), a.end());
  double sort_time = timer.elapsed();

  timer.restart();
  stk::radix_sort(b);
  double radix_time = timer.elapsed();

  std::cout << "Compare std::sort to stk::radix_sort => deque<size_t>" << std::endl;
  std::cout << "Size = 2^24" << std::endl;

  std::cout << "        std::sort : " << sort_time << " secs " << std::endl;
  std::cout << "  stk::radix_sort : " << radix_time << " secs, " << sort_time/radix_time << "x speedup"  << std::endl;
}

