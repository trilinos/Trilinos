/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

/// \file Test_Common_Sorting.hpp
/// \brief Tests for radixSort and bitonicSort in KokkoKernels_Sorting.hpp

#ifndef KOKKOSKERNELS_SORTINGTEST_HPP
#define KOKKOSKERNELS_SORTINGTEST_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>
#include <KokkosKernels_SimpleUtils.hpp>
#include <KokkosKernels_Utils.hpp>
#include <KokkosKernels_Sorting.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_Complex.hpp>
#include <cstdlib>

//Generate n randomized counts with mean <avg>.
//Then prefix-sum into randomOffsets.
//This simulates a CRS rowmap or other batched sorting scenario
template<typename OrdView, typename ExecSpace>
size_t generateRandomOffsets(OrdView randomCounts, OrdView randomOffsets, size_t n, size_t avg)
{
  srand(54321);
  auto countsHost = Kokkos::create_mirror_view(randomCounts);
  size_t total = 0;
  for(size_t i = 0; i < n; i++)
  {
    if(avg == 0)
      countsHost(i) = 0;
    else
      countsHost(i) = 0.5 + rand() % (avg * 2);
    total += countsHost(i);
  }
  Kokkos::deep_copy(randomCounts, countsHost);
  Kokkos::deep_copy(randomOffsets, randomCounts);
  KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<OrdView, ExecSpace>(n, randomOffsets);
  return total;
}

struct Coordinates
{
  KOKKOS_INLINE_FUNCTION Coordinates()
  {
     x = 0;
     y = 0;
     z = 0;
  }
  KOKKOS_INLINE_FUNCTION Coordinates(double x_, double y_, double z_)
  {
     x = x_;
     y = y_;
     z = z_;
  }
  double x;
  double y;
  double z;
};

/*
getRandom generates a random object for sorting.
the general version just produces integers - below
are some specializations
*/

template<typename T>
T getRandom()
{
  return rand() % Kokkos::ArithTraits<T>::max();
}

//Generate a uniform double between (-5, 5)
template<>
double getRandom<double>()
{
  return -5 + (10.0 * rand()) / RAND_MAX;
}

template<>
Coordinates getRandom<Coordinates>()
{
  return Coordinates(getRandom<double>(), getRandom<double>(), getRandom<double>());
}

//Specialize for Kokkos::complex, with the real and imaginary parts different
template<typename Key, typename Value>
struct kvHash
{
  Value operator()(const Key& k)
  {
    return (Value) (3 * k + 4);
  }
};

template<typename Key>
struct kvHash<Key, Kokkos::complex<double>>
{
  Kokkos::complex<double> operator()(const Key& k)
  {
    return Kokkos::complex<double>(3 * k + 4, k - 10.4);
  }
};

template<typename View>
void fillRandom(View v)
{
  srand(12345);
  typedef typename View::value_type Value;
  auto vhost = Kokkos::create_mirror_view(v);
  for(size_t i = 0; i < v.extent(0); i++)
    vhost(i) = getRandom<Value>();
  Kokkos::deep_copy(v, vhost);
}

template<typename KeyView, typename ValView>
void fillRandom(KeyView k, ValView v)
{
  srand(23456);
  typedef typename KeyView::value_type Key;
  typedef typename ValView::value_type Value;
  auto khost = Kokkos::create_mirror_view(k);
  auto vhost = Kokkos::create_mirror_view(v);
  for(size_t i = 0; i < v.extent(0); i++)
  {
    khost(i) = getRandom<Key>();
    vhost(i) = kvHash<Key, Value>()(khost(i));
  }
  Kokkos::deep_copy(k, khost);
  Kokkos::deep_copy(v, vhost);
}

template<typename KeyView, typename OrdView>
struct TestSerialRadixFunctor
{
  using Key = typename KeyView::value_type;
  using UnsignedKey = typename std::make_unsigned<Key>::type;

  TestSerialRadixFunctor(KeyView& keys_, KeyView& keysAux_, OrdView& counts_, OrdView& offsets_)
    : keys(keys_), keysAux(keysAux_), counts(counts_), offsets(offsets_)
  {}
  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    int off = offsets(i);
    KokkosKernels::Impl::SerialRadixSort<int, UnsignedKey>(
        (UnsignedKey*) keys.data() + off, (UnsignedKey*) keysAux.data() + off, counts(i));
  }
  KeyView keys;
  KeyView keysAux;
  OrdView counts;
  OrdView offsets;
};

template<typename KeyView, typename ValView, typename OrdView>
struct TestSerialRadix2Functor
{
  //Sort by keys, while permuting values
  using Key = typename KeyView::value_type;
  using UnsignedKey = typename std::make_unsigned<Key>::type;
  using Value = typename ValView::value_type;

  TestSerialRadix2Functor(KeyView& keys_, KeyView& keysAux_, ValView& values_, ValView& valuesAux_, OrdView& counts_, OrdView& offsets_)
    : keys(keys_), keysAux(keysAux_), values(values_), valuesAux(valuesAux_), counts(counts_), offsets(offsets_)
  {}
  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    int off = offsets(i);
    KokkosKernels::Impl::SerialRadixSort2<int, UnsignedKey, Value>(
        (UnsignedKey*) keys.data() + off, (UnsignedKey*) keysAux.data() + off,
        values.data() + off, valuesAux.data() + off, counts(i));
  }
  KeyView keys;
  KeyView keysAux;
  ValView values;
  ValView valuesAux;
  OrdView counts;
  OrdView offsets;
};

template<typename ExecSpace, typename Key>
void testSerialRadixSort(size_t k, size_t subArraySize)
{
  //Create a view of randomized data
  typedef typename ExecSpace::memory_space mem_space;
  typedef Kokkos::View<int*, mem_space> OrdView;
  typedef Kokkos::View<Key*, mem_space> KeyView;
  OrdView counts("Subarray Sizes", k);
  OrdView offsets("Subarray Offsets", k);
  //Generate k sub-array sizes, each with size about 20
  size_t n = generateRandomOffsets<OrdView, ExecSpace>(counts, offsets, k, subArraySize);
  auto countsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), counts);
  auto offsetsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), offsets);
  KeyView keys("Radix sort testing data", n);
  fillRandom(keys);
  //Sort using std::sort on host to do correctness test
  Kokkos::View<Key*, Kokkos::HostSpace> gold("Host sorted", n);
  Kokkos::deep_copy(gold, keys);
  for(size_t i = 0; i < k; i++)
  {
    Key* begin = gold.data() + offsetsHost(i);
    Key* end = begin + countsHost(i);
    std::sort(begin, end);
  }
  KeyView keysAux("Radix sort aux data", n);
  //Run the sorting on device in all sub-arrays in parallel
  typedef Kokkos::RangePolicy<ExecSpace> range_policy;
  Kokkos::parallel_for(range_policy(0, k),
        TestSerialRadixFunctor<KeyView, OrdView>(keys, keysAux, counts, offsets));
  //Copy result to host
  auto keysHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), keys);
  for(size_t i = 0; i < n; i++)
  {
    ASSERT_EQ(keysHost(i), gold(i));
  }
}

template<typename ExecSpace, typename Key, typename Value>
void testSerialRadixSort2(size_t k, size_t subArraySize)
{
  //Create a view of randomized data
  typedef typename ExecSpace::memory_space mem_space;
  typedef Kokkos::View<int*, mem_space> OrdView;
  typedef Kokkos::View<Key*, mem_space> KeyView;
  typedef Kokkos::View<Value*, mem_space> ValView;
  OrdView counts("Subarray Sizes", k);
  OrdView offsets("Subarray Offsets", k);
  //Generate k sub-array sizes, each with size about 20
  size_t n = generateRandomOffsets<OrdView, ExecSpace>(counts, offsets, k, subArraySize);
  auto countsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), counts);
  auto offsetsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), offsets);
  KeyView keys("Radix test keys", n);
  ValView data("Radix test data", n);
  //The keys are randomized
  fillRandom(keys, data);
  KeyView keysAux("Radix sort aux keys", n);
  ValView dataAux("Radix sort aux data", n);
  //Run the sorting on device in all sub-arrays in parallel
  typedef Kokkos::RangePolicy<ExecSpace> range_policy;
  //Deliberately using a weird number for vector length
  Kokkos::parallel_for(range_policy(0, k),
        TestSerialRadix2Functor<KeyView, ValView, OrdView>(keys, keysAux, data, dataAux, counts, offsets));
  //Sort using std::sort on host to do correctness test
  Kokkos::View<Key*, Kokkos::HostSpace> gold("Host sorted", n);
  Kokkos::deep_copy(gold, keys);
  for(size_t i = 0; i < k; i++)
  {
    Key* begin = gold.data() + offsetsHost(i);
    Key* end = begin + countsHost(i);
    std::sort(begin, end);
  }
  //Copy results to host
  auto keysHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), keys);
  auto dataHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), data);
  //Make sure keys are sorted exactly (stability of sort doesn't matter)
  for(size_t i = 0; i < n; i++)
  {
    ASSERT_EQ(keysHost(i), gold(i));
  }
  //Make sure the hashes of each key still matches the corresponding value
  for(size_t i = 0; i < n; i++)
  {
    auto correctHash = kvHash<Key, Value>()(keysHost(i));
    ASSERT_EQ(dataHost(i), correctHash);
  }
}

template<typename ValView, typename OrdView>
struct TestTeamBitonicFunctor
{
  typedef typename ValView::value_type Value;

  TestTeamBitonicFunctor(ValView& values_, OrdView& counts_, OrdView& offsets_)
    : values(values_), counts(counts_), offsets(offsets_)
  {}

  template<typename TeamMem>
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMem t) const
  {
    int i = t.league_rank();
    KokkosKernels::Impl::TeamBitonicSort<int, Value, TeamMem>(values.data() + offsets(i), counts(i), t);
  }

  ValView values;
  OrdView counts;
  OrdView offsets;
};

template<typename KeyView, typename ValView, typename OrdView>
struct TestTeamBitonic2Functor
{
  typedef typename KeyView::value_type Key;
  typedef typename ValView::value_type Value;

  TestTeamBitonic2Functor(KeyView& keys_, ValView& values_, OrdView& counts_, OrdView& offsets_)
    : keys(keys_), values(values_), counts(counts_), offsets(offsets_)
  {}

  template<typename TeamMem>
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMem t) const
  {
    int i = t.league_rank();
    KokkosKernels::Impl::TeamBitonicSort2<int, Key, Value, TeamMem>(keys.data() + offsets(i), values.data() + offsets(i), counts(i), t);
  }

  KeyView keys;
  ValView values;
  OrdView counts;
  OrdView offsets;
};

template<typename ExecSpace, typename Scalar>
void testTeamBitonicSort(size_t k, size_t subArraySize)
{
  //Create a view of randomized data
  typedef typename ExecSpace::memory_space mem_space;
  typedef Kokkos::View<int*, mem_space> OrdView;
  typedef Kokkos::View<Scalar*, mem_space> ValView;
  OrdView counts("Subarray Sizes", k);
  OrdView offsets("Subarray Offsets", k);
  //Generate k sub-array sizes, each with size about 20
  size_t n = generateRandomOffsets<OrdView, ExecSpace>(counts, offsets, k, subArraySize);
  auto countsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), counts);
  auto offsetsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), offsets);
  ValView data("Bitonic sort testing data", n);
  fillRandom(data);
  //Run the sorting on device in all sub-arrays in parallel
  Kokkos::parallel_for(Kokkos::TeamPolicy<ExecSpace>(k, Kokkos::AUTO()),
      TestTeamBitonicFunctor<ValView, OrdView>(data, counts, offsets));
  //Copy result to host
  auto dataHost = Kokkos::create_mirror_view(data);
  Kokkos::deep_copy(dataHost, data);
  //Sort using std::sort on host to do correctness test
  Kokkos::View<Scalar*, Kokkos::HostSpace> gold("Host sorted", n);
  Kokkos::deep_copy(gold, data);
  for(size_t i = 0; i < k; i++)
  {
    Scalar* begin = gold.data() + offsetsHost(i);
    Scalar* end = begin + countsHost(i);
    std::sort(begin, end);
  }
  for(size_t i = 0; i < n; i++)
  {
    ASSERT_EQ(dataHost(i), gold(i));
  }
}

template<typename ExecSpace, typename Key, typename Value>
void testTeamBitonicSort2(size_t k, size_t subArraySize)
{
  //Create a view of randomized data
  typedef typename ExecSpace::memory_space mem_space;
  typedef Kokkos::View<int*, mem_space> OrdView;
  typedef Kokkos::View<Key*, mem_space> KeyView;
  typedef Kokkos::View<Value*, mem_space> ValView;
  OrdView counts("Subarray Sizes", k);
  OrdView offsets("Subarray Offsets", k);
  //Generate k sub-array sizes, each with size about 20
  size_t n = generateRandomOffsets<OrdView, ExecSpace>(counts, offsets, k, subArraySize);
  auto countsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), counts);
  auto offsetsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), offsets);
  KeyView keys("Bitonic test keys", n);
  ValView data("Bitonic test data", n);
  //The keys are randomized
  fillRandom(keys, data);
  //Run the sorting on device in all sub-arrays in parallel, just using vector loops
  //Deliberately using a weird number for vector length
  Kokkos::parallel_for(Kokkos::TeamPolicy<ExecSpace>(k, Kokkos::AUTO()),
      TestTeamBitonic2Functor<KeyView, ValView, OrdView>(keys, data, counts, offsets));
  //Sort using std::sort on host to do correctness test
  Kokkos::View<Key*, Kokkos::HostSpace> gold("Host sorted", n);
  Kokkos::deep_copy(gold, keys);
  for(size_t i = 0; i < k; i++)
  {
    Key* begin = gold.data() + offsetsHost(i);
    Key* end = begin + countsHost(i);
    std::sort(begin, end);
  }
  //Copy results to host
  auto keysHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), keys);
  auto dataHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), data);
  //Make sure keys are sorted exactly (stability of sort doesn't matter)
  for(size_t i = 0; i < n; i++)
  {
    ASSERT_EQ(keysHost(i), gold(i));
  }
  //Make sure the hashes of each key still matches the corresponding value
  for(size_t i = 0; i < n; i++)
  {
    auto correctHash = kvHash<Key, Value>()(keysHost(i));
    ASSERT_EQ(dataHost(i), correctHash);
  }
}

template<typename View>
struct CheckSortedFunctor
{
  CheckSortedFunctor(View& v_)
    : v(v_) {}
  KOKKOS_INLINE_FUNCTION void operator()(int i, int& lval) const
  {
    if(v(i) > v(i + 1))
      lval = 0;
  }
  View v;
};

template<typename ExecSpace, typename Scalar>
void testBitonicSort(size_t n)
{
  //Create a view of randomized data
  typedef typename ExecSpace::memory_space mem_space;
  typedef Kokkos::View<Scalar*, mem_space> ValView;
  ValView data("Bitonic sort testing data", n);
  fillRandom(data);
  KokkosKernels::Impl::bitonicSort<ValView, ExecSpace, int>(data);
  int ordered = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, n - 1),
      CheckSortedFunctor<ValView>(data), Kokkos::Min<int>(ordered));
  ASSERT_TRUE(ordered);
}

//Check that the view is weakly ordered according to Comparator:
//Comparator never says that element i+1 belongs before element i.
template<typename View, typename Comparator>
struct CheckOrderedFunctor
{
  CheckOrderedFunctor(View& v_)
    : v(v_) {}
  KOKKOS_INLINE_FUNCTION void operator()(int i, int& lval) const
  {
    Comparator comp;
    if(comp(v(i + 1), v(i)))
      lval = 0;
  }
  View v;
};

template<typename Scalar>
struct CompareDescending
{
  KOKKOS_INLINE_FUNCTION bool operator()(const Scalar lhs, const Scalar rhs) const
  {
    return lhs > rhs;
  }
};

template<typename ExecSpace>
void testBitonicSortDescending()
{
  typedef char Scalar;
  typedef CompareDescending<Scalar> Comp;
  //Create a view of randomized data
  typedef typename ExecSpace::memory_space mem_space;
  typedef Kokkos::View<Scalar*, mem_space> ValView;
  size_t n = 12521;
  ValView data("Bitonic sort testing data", n);
  fillRandom(data);
  KokkosKernels::Impl::bitonicSort<ValView, ExecSpace, int, Comp>(data);
  int ordered = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, n - 1),
      CheckOrderedFunctor<ValView, Comp>(data), Kokkos::Min<int>(ordered));
  ASSERT_TRUE(ordered);
}

struct LexCompare
{
  KOKKOS_INLINE_FUNCTION bool operator()(const Coordinates lhs, const Coordinates rhs) const
  {
    if(lhs.x < rhs.x)
      return true;
    else if(lhs.x > rhs.x)
      return false;
    else if(lhs.y < rhs.y)
      return true;
    else if(lhs.y > rhs.y)
      return false;
    else if(lhs.z < rhs.z)
      return true;
    return false;
  }
};

template<typename ExecSpace>
void testBitonicSortLexicographic()
{
  typedef Coordinates Scalar;
  //Create a view of randomized data
  typedef typename ExecSpace::memory_space mem_space;
  typedef Kokkos::View<Scalar*, mem_space> ValView;
  size_t n = 9521;
  ValView data("Bitonic sort testing data", n);
  fillRandom(data);
  KokkosKernels::Impl::bitonicSort<ValView, ExecSpace, int, LexCompare>(data);
  int ordered = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, n - 1),
      CheckOrderedFunctor<ValView, LexCompare>(data), Kokkos::Min<int>(ordered));
  ASSERT_TRUE(ordered);
}

TEST_F(TestCategory, common_serial_radix) {
  //Test serial radix over some contiguous small arrays
  //1st arg is #arrays, 2nd arg is max subarray size
  size_t numArrays = 100;
  for(size_t arrayMax = 0; arrayMax < 1000; arrayMax = 1 + 4 * arrayMax)
  {
    testSerialRadixSort<TestExecSpace, char>(numArrays, arrayMax);
    testSerialRadixSort<TestExecSpace, int>(numArrays, arrayMax);
  }
}

TEST_F(TestCategory, common_serial_radix2) {
  //Test serial radix over some contiguous small arrays
  //1st arg is #arrays, 2nd arg is max subarray size
  size_t numArrays = 100;
  for(size_t arrayMax = 0; arrayMax < 1000; arrayMax = 1 + 4 * arrayMax)
  {
    testSerialRadixSort2<TestExecSpace, char, int>(numArrays, arrayMax);
    testSerialRadixSort2<TestExecSpace, int, double>(numArrays, arrayMax);
    testSerialRadixSort2<TestExecSpace, int, Kokkos::complex<double>>(numArrays, arrayMax);
  }
}

TEST_F(TestCategory, common_team_bitonic) {
  //Test team-level bitonic over some contiguous medium arrays
  //1st arg is #arrays, 2nd arg is max subarray size
  size_t numArrays = 20;
  for(size_t arrayMax = 0; arrayMax < 10000; arrayMax = 1 + 4 * arrayMax)
  {
    testTeamBitonicSort<TestExecSpace, char>(numArrays, arrayMax);
    testTeamBitonicSort<TestExecSpace, int>(numArrays, arrayMax);
  }
}

TEST_F(TestCategory, common_team_bitonic2) {
  //Test team-level bitonic over some contiguous medium arrays
  //1st arg is #arrays, 2nd arg is max subarray size
  size_t numArrays = 20;
  for(size_t arrayMax = 0; arrayMax < 10000; arrayMax = 1 + 4 * arrayMax)
  {
    testTeamBitonicSort2<TestExecSpace, char, int>(numArrays, arrayMax);
    testTeamBitonicSort2<TestExecSpace, int, double>(numArrays, arrayMax);
    testTeamBitonicSort2<TestExecSpace, int, Kokkos::complex<double>>(numArrays, arrayMax);
  }
}

TEST_F( TestCategory, common_device_bitonic) {
  //Test device-level bitonic with some larger arrays
  testBitonicSort<TestExecSpace, char>(243743);
  testBitonicSort<TestExecSpace, char>(2157);
  testBitonicSort<TestExecSpace, char>(424);
  testBitonicSort<TestExecSpace, char>(5);
  testBitonicSort<TestExecSpace, int>(92314);
  testBitonicSort<TestExecSpace, int>(123);
  testBitonicSort<TestExecSpace, double>(60234);
  testBitonicSort<TestExecSpace, double>(53);
  //Test custom comparator: ">" instead of "<" to sort descending
  testBitonicSortDescending<TestExecSpace>();
  //Test custom comparator: lexicographic comparison of 3-element struct
  testBitonicSortLexicographic<TestExecSpace>();
}

#endif

