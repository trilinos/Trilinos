#include <Akri_DistributedVector.hpp>
#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

namespace krino {

// Dot must give bitwise identical results regardless of how a vector is ordered or distributed
// across processors, so that serial and parallel runs of the optimizers stay identical.  These
// vectors have one large entry and many small entries whose individual contributions are below
// half a double ulp of the large entry, but whose combined contribution is exactly one double ulp.
// Summing the small entries onto the large one one at a time loses them; summing them first (or on
// another processor) and then adding the large one keeps them.  Dot rounds its result to float, so
// both orders and distributions give the same result.

static double large_entry() { return std::ldexp(1., 20); }
static double small_entry() { return std::ldexp(1., -41); }
static constexpr size_t numSmallEntries = 512;

static std::vector<double> large_part() { return {large_entry()}; }
static std::vector<double> small_part() { return std::vector<double>(numSmallEntries, small_entry()); }

static std::vector<double> concatenate(const std::vector<double> & first, const std::vector<double> & second)
{
  std::vector<double> result = first;
  result.insert(result.end(), second.begin(), second.end());
  return result;
}

static std::vector<double> whole_vector_large_entry_first() { return concatenate(large_part(), small_part()); }
static std::vector<double> whole_vector_large_entry_last() { return concatenate(small_part(), large_part()); }
static double expected_dot_with_ones() { return float(large_entry() + numSmallEntries * small_entry()); }

static double serial_dot_with_ones(const std::vector<double> & values)
{
  const DistributedVector x(stk::parallel_machine_null(), values, {});
  const DistributedVector ones(stk::parallel_machine_null(), std::vector<double>(values.size(), 1.), {});
  return Dot(x, ones);
}

TEST(DistributedVectorDot, largeAndManySmallEntriesInDifferentOrders_sameResult)
{
  EXPECT_EQ(expected_dot_with_ones(), serial_dot_with_ones(whole_vector_large_entry_first()));
  EXPECT_EQ(expected_dot_with_ones(), serial_dot_with_ones(whole_vector_large_entry_last()));
}

TEST(DistributedVectorDot, entriesWithMixedSignsAndNegativeSum_correctResult)
{
  const DistributedVector x(stk::parallel_machine_null(), {-10., 4., 0.25, 1.e30}, {});
  const DistributedVector y(stk::parallel_machine_null(), {1., 1., -2., 0.}, {});
  EXPECT_EQ(-6.5, Dot(x, y));
}

static std::vector<double> random_values_with_mixed_signs_and_magnitudes(const size_t numValues)
{
  std::mt19937 generator(12345);
  std::uniform_real_distribution<double> log10Magnitude(-20., 5.);
  std::bernoulli_distribution isNegative(0.5);
  std::vector<double> values(numValues);
  for (auto & value : values)
    value = (isNegative(generator) ? -1. : 1.) * std::pow(10., log10Magnitude(generator));
  return values;
}

TEST(DistributedVectorDot, randomValuesWithMixedSignsAndMagnitudesShuffled_bitwiseIdentical)
{
  const std::vector<double> values = random_values_with_mixed_signs_and_magnitudes(10000);
  std::vector<double> shuffledValues = values;
  std::shuffle(shuffledValues.begin(), shuffledValues.end(), std::mt19937(54321));
  EXPECT_EQ(serial_dot_with_ones(values), serial_dot_with_ones(shuffledValues));
}

TEST(DistributedVectorDot, largeAndManySmallEntriesDistributedOnTwoProcs_sameAsSerialDot)
{
  const stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(comm) != 2) GTEST_SKIP();

  const std::vector<double> localValues = (stk::parallel_machine_rank(comm) == 0) ? large_part() : small_part();
  const DistributedVector x(comm, localValues, {});
  const DistributedVector ones(comm, std::vector<double>(localValues.size(), 1.), {});

  EXPECT_EQ(expected_dot_with_ones(), Dot(x, ones));
  EXPECT_EQ(serial_dot_with_ones(whole_vector_large_entry_first()), Dot(x, ones));
}

}
