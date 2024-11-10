// Copyright(C) 1999-2020, 2022, 2023, 2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ioss_Sort.h"
#include <cassert>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <fmt/core.h>
#include <iostream>
#include <random>
#include <stddef.h>
#include <stdint.h>
#include <vector>

namespace {
  const int sawtooth = 1;
  const int do_rand  = 2;
  const int stagger  = 3;
  const int plateau  = 4;
  const int shuffle  = 5;

  std::random_device rd;
  std::mt19937_64    rng(rd());

  template <typename INT> bool verify_sorted(const std::vector<INT> &v)
  {
    auto it = std::adjacent_find(v.begin(), v.end(), std::greater<INT>());
    if (it != v.end()) {
      std::cerr << "Unsorted at position " << it - v.begin() + 1 << "\n";
      return false;
    }
    return true;
  }

  std::vector<int64_t> generate_vector(int dist, size_t n, size_t m)
  {
    std::vector<int64_t> x;
    x.reserve(n);
    size_t i = 0, j = 0, k = 1;
    switch (dist) {
    case sawtooth:
      for (; i < n; i++) {
        x.push_back(i % m);
      }
      break;
    case do_rand:
      for (; i < n; i++) {
        x.push_back(rng() % m);
      }
      break;
    case stagger:
      for (; i < n; i++) {
        x.push_back((i * m + i) % n);
      }
      break;
    case plateau:
      for (; i < n; i++) {
        x.push_back(std::min(i, m));
      }
      break;
    case shuffle:
      for (; i < n; i++) {
        x.push_back((rng() % m) != 0u ? (j += 2) : (k += 2));
      }
      break;
    }
    return x;
  }
} // namespace

TEST_CASE("empty")
{
  std::vector<int64_t> x{};
  Ioss::sort(x);
  REQUIRE(verify_sorted(x));
}

TEST_CASE("single-element")
{
  int64_t              n = GENERATE(100, 1023, 1024, 1025, (2 << 16) - 1, 2 << 16, (2 << 16) + 1);
  std::vector<int64_t> x{n};
  Ioss::sort(x);
  REQUIRE(verify_sorted(x));
}

TEST_CASE("sort")
{
  auto   dist = GENERATE_COPY(sawtooth, do_rand, stagger, plateau, shuffle);
  size_t n    = GENERATE(100, 1023, 1024, 1025, (2 << 16) - 1, 2 << 16, (2 << 16) + 1);

  // 'm' shapes the values in the vector; all sorts are on vectors of size 'n'
  size_t m = GENERATE(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768,
                      65536, 131072, 262144);

  if (m < 2 * n) {
    std::vector<int64_t> x = generate_vector(dist, n, m);
    REQUIRE(x.size() == n);

    SECTION("output")
    {
      // Just printing output; no test...
      std::string type[] = {"sawtooth", "do_rand", "stagger", "plateau", "shuffle"};
      fmt::print("Size: {:8}, Shape = {:8}, Type = {:12}\n", n, m, type[dist - 1]);
    }

    DYNAMIC_SECTION("as generated" << n << m << dist)
    {
      Ioss::sort(x); // Copy of x
      REQUIRE(verify_sorted(x));
    }

    DYNAMIC_SECTION("reversed" << n << m << dist)
    {
      std::reverse(x.begin(), x.end()); // Reversed
      Ioss::sort(x);
      REQUIRE(verify_sorted(x));
    }

    DYNAMIC_SECTION("front-half reversed" << n << m << dist)
    {
      std::reverse(&x[0], &x[n / 2]); // Front half reversed
      Ioss::sort(x);
      REQUIRE(verify_sorted(x));
    }

    DYNAMIC_SECTION("back-half reversed" << n << m << dist)
    {
      std::reverse(&x[n / 2], &x[n]); // Back half reversed
      Ioss::sort(x);
      REQUIRE(verify_sorted(x));
      DYNAMIC_SECTION("already sorted" << n << m << dist)
      {
        //        REQUIRE(verify_sorted(x));
        Ioss::sort(x); // Already sorted
        REQUIRE(verify_sorted(x));
      }
    }

    DYNAMIC_SECTION("dithered" << n << m << dist)
    {
      for (size_t p = 0; p < n; p++) {
        x[p] += p % 5;
      }
      Ioss::sort(x); // Dithered
      REQUIRE(verify_sorted(x));
    }
  }
}
