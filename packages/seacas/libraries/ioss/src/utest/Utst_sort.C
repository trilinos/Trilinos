// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <Ioss_Sort.h>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>

namespace {
  template <typename INT> bool verify_sorted(std::vector<INT> v)
  {
    auto it = std::adjacent_find(v.begin(), v.end(), std::greater<INT>());
    if (it != v.end()) {
      std::cerr << "Unsorted at position " << it - v.begin() + 1 << "\n";
      return false;
    }
    return true;
  }
} // namespace

int main()
{
  std::random_device rd;
  std::mt19937_64    rng(rd());

  const int sawtooth = 1;
  const int do_rand  = 2;
  const int stagger  = 3;
  const int plateau  = 4;
  const int shuffle  = 5;

  for (size_t n : {100, 1023, 1024, 1025, (2 << 16) - 1, 2 << 16, (2 << 16) + 1}) {
    std::cerr << "\nSize: " << n << ": ";
    for (size_t m = 1; m < 2 * n; m *= 2) {
      std::cerr << m;
      for (auto dist : {sawtooth, do_rand, stagger, plateau, shuffle}) {
        std::cerr << ".";
        std::vector<int64_t> x(n);

        size_t i = 0, j = 0, k = 1;
        switch (dist) {
        case sawtooth:
          for (; i < n; i++) {
            x[i] = i % m;
          }
          break;
        case do_rand:
          for (; i < n; i++) {
            x[i] = rng() % m;
          }
          break;
        case stagger:
          for (; i < n; i++) {
            x[i] = (i * m + i) % n;
          }
          break;
        case plateau:
          for (; i < n; i++) {
            x[i] = std::min(i, m);
          }
          break;
        case shuffle:
          for (; i < n; i++) {
            x[i] = (rng() % m) != 0u ? (j += 2) : (k += 2);
          }
          break;
        }

        auto xsave = x;
        Ioss::qsort(x); // Copy of x
        assert(verify_sorted(x));

        x = xsave;
        std::reverse(x.begin(), x.end()); // Reversed
        Ioss::qsort(x);
        assert(verify_sorted(x));

        x = xsave;
        std::reverse(&x[0], &x[n / 2]); // Front half reversed
        Ioss::qsort(x);
        assert(verify_sorted(x));

        x = xsave;
        std::reverse(&x[n / 2], &x[n]); // Back half reversed
        Ioss::qsort(x);
        assert(verify_sorted(x));

        Ioss::qsort(x); // Already sorted
        assert(verify_sorted(x));

        x = xsave;
        for (size_t p = 0; p < n; p++) {
          x[p] += p % 5;
        }
        Ioss::qsort(x); // Dithered
        assert(verify_sorted(x));
      }
    }
  }
  std::cerr << "\nDone\n";
}
