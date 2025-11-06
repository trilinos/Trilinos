// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSKERNELS_PERFTEST_BENCHMARK_UTILS_HPP
#define KOKKOSKERNELS_PERFTEST_BENCHMARK_UTILS_HPP

namespace KokkosKernelsBenchmark {

class WrappedBool {
 public:
  WrappedBool(const bool &val) : val_(val) {}

  operator bool() const { return val_; }

 protected:
  bool val_;
};

class DieOnError : public WrappedBool {
 public:
  DieOnError(const bool &val) : WrappedBool(val) {}
};
class SkipOnError : public WrappedBool {
 public:
  SkipOnError(const bool &val) : WrappedBool(val) {}
};

}  // namespace KokkosKernelsBenchmark

#endif  // KOKKOSKERNELS_PERFTEST_BENCHMARK_UTILS_HPP