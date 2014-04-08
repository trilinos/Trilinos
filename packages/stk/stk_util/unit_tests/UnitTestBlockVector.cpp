
#include <stddef.h>                     // for size_t
#include <algorithm>                    // for equal
#include <boost/timer.hpp>              // for timer
#include <cstdlib>                      // for rand, srand, RAND_MAX
#include <deque>                        // for deque, operator!=
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/util/BlockVector.hpp>
#include <stk_util/util/CacheAlignedAllocator.hpp>
#include <stk_util/util/TrackingAllocator.hpp>  // for tracking_allocator
#include <stk_util/util/human_bytes.hpp>  // for human_bytes
#include <string>                       // for operator<<, string
#include <utility>                      // for pair, operator==, make_pair
#include <vector>                       // for vector
#include "gtest/gtest.h"                // for AssertHelper, ASSERT_EQ, etc
#include "stk_util/util/AllocatorMemoryUsage.hpp"


STKUNIT_UNIT_TEST( block_vector, basic )
{
  typedef stk::block_vector<unsigned> bvector;
  typedef std::vector<unsigned> vector;

  const size_t block_size = bvector::block_size;

  // test default constructor
  {
    bvector v;
    ASSERT_EQ(0u, v.size());
    ASSERT_EQ(0u, v.capacity());
    ASSERT_TRUE(v.empty());
  }

  {
    // test constructor
    bvector v(block_size*2+1, 1u);
    ASSERT_EQ(block_size*2u+1u, v.size());
    ASSERT_EQ(block_size*3, v.capacity());

    for (size_t i=0, ie=v.size(); i<ie; ++i) EXPECT_EQ(1u, v[i]);

    // test copy constructor
    bvector w(v);
    ASSERT_EQ(v.size(), w.size());
    ASSERT_EQ(v.capacity(), w.capacity());
    ASSERT_TRUE(v == w);

    // test assignment operator
    bvector x;
    x = w;
    ASSERT_EQ(w.size(), x.size());
    ASSERT_EQ(w.capacity(), x.capacity());
    ASSERT_TRUE(w == x);

    // test comparisons
    x.push_back(1);
    ASSERT_TRUE(w < x);
    ASSERT_TRUE(x > w);
    ASSERT_TRUE(x != w);

    w.back() = 2;
    ASSERT_TRUE(v <= w);
    ASSERT_TRUE(w >= v);
    ASSERT_TRUE(v < w);

    // test clear
    w.clear();
    ASSERT_EQ(0u, w.size());
    ASSERT_EQ(block_size*3, w.capacity());

    // test resize
    w.resize(block_size*3, 5);
    ASSERT_EQ(block_size*3, w.size());
    ASSERT_EQ(block_size*3, w.capacity());
    // test forward iterator
    for (bvector::iterator itr = w.begin(), end_itr =w.end(); itr != end_itr; ++itr) {
      EXPECT_EQ(5u, *itr);
    }

    w.resize(block_size-1);
    ASSERT_EQ(block_size-1u, w.size());
    // test reverse iterator
    for (bvector::const_reverse_iterator itr = w.rbegin(), end_itr =w.rend(); itr != end_itr; ++itr) {
      ASSERT_EQ(5u, *itr);
    }

    // test shrink_to_fit
    w.shrink_to_fit();
    ASSERT_EQ(block_size, w.capacity());

    ASSERT_EQ(5u, w.front());

    // test assign
    {
      vector vec(2000);
      for (size_t i=0, ie=vec.size(); i<ie; ++i) {
        vec[i] = i;
      }
      x.assign(vec.begin(), vec.end());

      ASSERT_EQ(vec.size(), x.size());
      for (size_t i=0, ie=x.size(); i<ie; ++i) {
        ASSERT_EQ(i, x.at(i));
      }
    }

    // test swap
    v = x;
    w.clear();
    w.swap(v);

    ASSERT_TRUE(x == w);
    ASSERT_TRUE(v.empty());

    // test pop_back
    for (size_t i=0, ie = w.size(); i<ie; ++i) w.pop_back();
    ASSERT_TRUE(w.empty());

    {
      size_t curr_size = x.size();
      x.pop_back(curr_size/2);

      ASSERT_EQ(curr_size/2, x.size());
    }
  }

  // compare behavior to vector
  {
    vector v;
    bvector b;

    // cover multiple blocks
    for (int i=0; i<2000; ++i) {
      v.push_back(i);
      b.push_back(i);
    }

    ASSERT_EQ( b.size(), v.size());
    ASSERT_TRUE( std::equal(b.cbegin(), b.cend(), v.begin()));

    b.clear();
    b.assign(v.begin(), v.end());
    ASSERT_EQ( b.size(), v.size());
    ASSERT_TRUE( std::equal(b.cbegin(), b.cend(), v.begin()));
    ASSERT_TRUE( std::equal(b.rbegin(), b.rend(), v.rbegin()));
    ASSERT_TRUE( std::equal(v.rbegin(), v.rend(), b.rbegin()));

    v.insert(v.begin()+1000, 5000);
    b.insert(b.begin()+1000, 5000);
    ASSERT_EQ( b.size(), v.size());
    ASSERT_TRUE( std::equal(b.cbegin(), b.cend(), v.begin()));

    v.erase(v.begin()+1000);
    b.erase(b.begin()+1000);
    ASSERT_EQ( b.size(), v.size());
    ASSERT_TRUE( std::equal(b.cbegin(), b.cend(), v.begin()));

    v.insert(v.begin()+500, 1000, 0);
    b.insert(b.begin()+500, 1000, 0);
    ASSERT_EQ( b.size(), v.size());
    ASSERT_TRUE( std::equal(b.cbegin(), b.cend(), v.begin()));

    v.erase(v.begin()+500, v.begin()+1500);
    b.erase(b.begin()+500, b.begin()+1500);
    ASSERT_EQ( b.size(), v.size());
    ASSERT_TRUE( std::equal(b.cbegin(), b.cend(), v.begin()));

    unsigned array[5] = { 100u, 200u, 300u, 400u, 500u };

    v.insert(v.begin(), array, array+5);
    b.insert(b.begin(), array, array+5);
    ASSERT_EQ( b.size(), v.size());
    ASSERT_TRUE( std::equal(b.cbegin(), b.cend(), v.begin()));

    for (int i=0; i<5; ++i) {
      v.erase(v.begin());
      b.erase(b.begin());
    }
    ASSERT_EQ( b.size(), v.size());
    ASSERT_TRUE( std::equal(b.cbegin(), b.cend(), v.begin()));
  }
}

namespace {

template <typename T>
struct Random
{
  T operator()() const
  { return static_cast<T>(std::rand()); }
};

template <>
struct Random<float>
{
  float operator()() const
  { return static_cast<float>(std::rand()) / RAND_MAX; }
};

template <>
struct Random<double>
{
  double operator()() const
  { return static_cast<double>(std::rand()) / RAND_MAX; }
};


template <typename T, typename U>
inline std::pair<T,U> & operator+=(std::pair<T,U> & lhs, const std::pair<T,U> & rhs)
{
  lhs.first += rhs.first;
  lhs.second += rhs.second;
  return lhs;
}

template <typename T, typename U>
struct Random< std::pair<T,U> >
{
  std::pair<T,U> operator()() const
  {
    Random<T> trand;
    Random<U> urand;
    return std::make_pair(trand(),urand());
  }
};

struct RunTest
{
  double fill_time;
  double iterator_time;
  double for_time;

  RunTest()
    : fill_time()
    , iterator_time()
    , for_time()
  {}

  template <typename Array>
  void apply( Array & a ,size_t num)
  {
    typedef typename Array::value_type value_type;
    boost::timer timer;

    timer.restart();
    Random<typename Array::value_type> random;
    for (size_t i=0; i<num; ++i) {
      a.push_back(random());
    }
    fill_time += timer.elapsed();

    timer.restart();

    value_type sum = value_type();
    for (typename Array::const_iterator itr = a.begin(), eitr = a.end(); itr != eitr; ++itr) {
      sum += *itr;
    }
    iterator_time += timer.elapsed();

    timer.restart();
    value_type sum2 = value_type();
    for (size_t i=0, ie=a.size(); i < ie; ++i) {
      sum2 += a[i];
    }
    for_time += timer.elapsed();

    EXPECT_TRUE( sum == sum2);
  }

  double total_time() const
  {
    return fill_time + iterator_time + for_time;
  }
};

struct VectorTag {};
struct DequeTag {};
struct BlockVectorTag {};

template < typename Vector, typename VectorMemory
          ,typename BlockVector, typename BlockVectorMemory
          ,typename Deque, typename DequeMemory
         >
void run_performance_test(std::string const& type, size_t size, size_t stride, size_t size_base, size_t size_exp)
{
  VectorMemory::reset();
  BlockVectorMemory::reset();
  DequeMemory::reset();

  RunTest vector_test;
  RunTest deque_test;
  RunTest block_vector_test;

  BlockVector b;
  Vector v;
  Deque d;

  std::cout << std::endl;
  std::cout << "   ";
  for (size_t i=0; i<=size; i += stride)
  {
    srand(i);
    block_vector_test.apply(b, stride);

    srand(i);
    vector_test.apply(v, stride);

    srand(i);
    deque_test.apply(d, stride);
    int percent = static_cast<int>((100.0 * i) / size);

    std::cout << "\b\b\b" << percent << '%';
  }

  std::cout << std::endl;

  std::cout.precision(2);

  std::cout << "Comparison to vector<" << type << ">" << std::endl;
  std::cout << "Size = " << size_base << "^" << size_exp << std::endl;
  std::cout << "Minimal Memory : " << stk::human_bytes(size*sizeof(typename Vector::value_type)) << std::endl;
  std::cout << std::endl;

  std::cout << "  vector<" << type << "> " << std::endl;
  std::cout << "        push_back : " << vector_test.fill_time / vector_test.fill_time << "x" << std::endl;
  std::cout << "         iterator : " << vector_test.iterator_time / vector_test.iterator_time << "x" << std::endl;
  std::cout << "         for loop : " << vector_test.for_time / vector_test.for_time << "x" << std::endl;
  std::cout << "       TOTAL time : " << vector_test.total_time() << " seconds" << std::endl;
  std::cout << "           MEMORY : " << stk::human_bytes(VectorMemory::peak_memory()) << std::endl;
  std::cout << std::endl;

  std::cout << "  block_vector<" << type << "> " << std::endl;
  std::cout << "        push_back : " << block_vector_test.fill_time / vector_test.fill_time << "x" << std::endl;
  std::cout << "         iterator : " << block_vector_test.iterator_time / vector_test.iterator_time << "x" << std::endl;
  std::cout << "         for loop : " << block_vector_test.for_time / vector_test.for_time << "x" << std::endl;
  std::cout << "       TOTAL time : " << block_vector_test.total_time() << " seconds" << std::endl;
  std::cout << "           MEMORY : " << stk::human_bytes(BlockVectorMemory::peak_memory()) << std::endl;
  std::cout << std::endl;

  std::cout << "  deque<" << type << "> " << std::endl;
  std::cout << "        push_back : " << deque_test.fill_time / vector_test.fill_time << "x" << std::endl;
  std::cout << "         iterator : " << deque_test.iterator_time / vector_test.iterator_time << "x" << std::endl;
  std::cout << "         for loop : " << deque_test.for_time / vector_test.for_time << "x" << std::endl;
  std::cout << "       TOTAL time : " << deque_test.total_time() << " seconds" << std::endl;
  std::cout << "           MEMORY : " << stk::human_bytes(DequeMemory::peak_memory()) << std::endl;
  std::cout << std::endl;

  std::cout << std::endl;

}

} // namespace

STKUNIT_UNIT_TEST( block_vector, DISABLED_performance_int )
{
  typedef std::vector<int, stk::tracking_allocator<int, VectorTag> > Vector;
  typedef std::deque<int, stk::tracking_allocator<int, DequeTag> > Deque;
  typedef stk::block_vector<int, 1024, stk::cache_aligned_allocator<int, BlockVectorTag> > BlockVector;

  typedef stk::allocator_memory_usage<VectorTag>      VectorMemory;
  typedef stk::allocator_memory_usage<DequeTag>       DequeMemory;
  typedef stk::allocator_memory_usage<BlockVectorTag> BlockVectorMemory;


  BlockVector b;
  Vector v;
  Deque d;

  const size_t size   = 100000000;
  const size_t stride = 1000000;

  run_performance_test<Vector,VectorMemory,BlockVector,BlockVectorMemory,Deque,DequeMemory>("int", size, stride, 10, 8);
}

STKUNIT_UNIT_TEST( block_vector, DISABLED_performance_size_t )
{
  typedef std::vector<size_t, stk::tracking_allocator<size_t, VectorTag> > Vector;
  typedef std::deque<size_t, stk::tracking_allocator<size_t, DequeTag> > Deque;
  typedef stk::block_vector<size_t, 512, stk::cache_aligned_allocator<size_t, BlockVectorTag> > BlockVector;

  typedef stk::allocator_memory_usage<VectorTag>      VectorMemory;
  typedef stk::allocator_memory_usage<DequeTag>       DequeMemory;
  typedef stk::allocator_memory_usage<BlockVectorTag> BlockVectorMemory;


  BlockVector b;
  Vector v;
  Deque d;

  const size_t size   = 100000000;
  const size_t stride = 1000000;

  run_performance_test<Vector,VectorMemory,BlockVector,BlockVectorMemory,Deque,DequeMemory>("size_t", size, stride, 10, 8);
}

STKUNIT_UNIT_TEST( block_vector, DISABLED_performance_pair_double_size_t )
{
  typedef std::pair<double,size_t> MyPair;
  typedef std::vector<MyPair, stk::tracking_allocator<MyPair, VectorTag> > Vector;
  typedef std::deque<MyPair, stk::tracking_allocator<MyPair, DequeTag> > Deque;
  typedef stk::block_vector<MyPair, 256, stk::cache_aligned_allocator<MyPair, BlockVectorTag> > BlockVector;

  typedef stk::allocator_memory_usage<VectorTag>      VectorMemory;
  typedef stk::allocator_memory_usage<DequeTag>       DequeMemory;
  typedef stk::allocator_memory_usage<BlockVectorTag> BlockVectorMemory;


  BlockVector b;
  Vector v;
  Deque d;

  const size_t shift = 24;
  const size_t size   = 1ull << shift;
  const size_t stride = size >> 10 ;

  run_performance_test<Vector,VectorMemory,BlockVector,BlockVectorMemory,Deque,DequeMemory>("pair<double,size_t>", size, stride, 2, shift);
}
