/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
*/

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