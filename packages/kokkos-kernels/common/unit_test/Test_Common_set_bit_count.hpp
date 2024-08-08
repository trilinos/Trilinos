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

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include "KokkosKernels_BitUtils.hpp"
#include "KokkosKernels_SimpleUtils.hpp"
#include "KokkosKernels_PrintUtils.hpp"
#include <string>
#include <stdexcept>

// const char *input_filename = "sherman1.mtx";
// const char *input_filename = "Si2.mtx";
// const char *input_filename = "wathen_30_30.mtx";
// const size_t expected_num_cols = 9906;

using namespace KokkosKernels;
using namespace KokkosKernels::Impl;

namespace Test {

template <typename view_type>
struct ppctest {
  view_type view;
  typename view_type::non_const_type out_view;
  ppctest(view_type view_, typename view_type::non_const_type out_view_) : view(view_), out_view(out_view_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t row) const { out_view(row) = pop_count(view(row)); }
};

template <typename view_type>
struct ppccheck {
  view_type view;
  typename view_type::non_const_type out_view;
  ppccheck(view_type view_, typename view_type::non_const_type out_view_) : view(view_), out_view(out_view_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t row) const {
    typename view_type::non_const_value_type myval = view(row);
    int num_el2                                    = 0;
    for (; myval; num_el2++) {
      myval = myval & (myval - 1);  // clear the least significant bit set
    }
    out_view(row) = num_el2;
  }
};

template <typename view_type, typename execution_space>
view_type get_array_bit_count(view_type view) {
  typename view_type::non_const_type out_view("out", view.extent(0));

  typedef Kokkos::RangePolicy<execution_space> my_exec_space;
  Kokkos::parallel_for("KokkosKernels::Common::Test::GetArrayBitCount", my_exec_space(0, view.extent(0)),
                       ppctest<view_type>(view, out_view));
  Kokkos::fence();
  return out_view;
}

template <typename view_type, typename execution_space>
view_type check_array_bit_count(view_type view) {
  typename view_type::non_const_type out_view("out", view.extent(0));

  typedef Kokkos::RangePolicy<execution_space> my_exec_space;
  Kokkos::parallel_for("KokkosKernels::Common::Test::CheckArrayBitCount", my_exec_space(0, view.extent(0)),
                       ppccheck<view_type>(view, out_view));
  Kokkos::fence();
  return out_view;
}

template <typename view_type>
struct ffstest {
  view_type view;
  typename view_type::non_const_type out_view;
  ffstest(view_type view_, typename view_type::non_const_type out_view_) : view(view_), out_view(out_view_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t row) const {
    if (view(row) != 0) {
      out_view(row) = least_set_bit(view(row)) - 1;
    } else
      out_view(row) = 0;
  }
};

template <typename view_type>
struct ffscheck {
  view_type view;
  typename view_type::non_const_type out_view;
  ffscheck(view_type view_, typename view_type::non_const_type out_view_) : view(view_), out_view(out_view_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t row) const {
    typename view_type::non_const_value_type myval = view(row);
    typename view_type::non_const_value_type unit  = 1;
    out_view(row)                                  = 0;
    for (int i = 0; i < 64; ++i) {
      if (myval & unit << i) {
        out_view(row) = i;
        break;
      }
    }
  }
};

template <typename view_type, typename execution_space>
view_type get_ffs(view_type view) {
  typename view_type::non_const_type out_view("out", view.extent(0));

  typedef Kokkos::RangePolicy<execution_space> my_exec_space;
  Kokkos::parallel_for("KokkosKernels::Common::Test::GetFFS", my_exec_space(0, view.extent(0)),
                       ffstest<view_type>(view, out_view));
  Kokkos::fence();
  return out_view;
}

template <typename view_type, typename execution_space>
view_type check_ffs(view_type view) {
  typename view_type::non_const_type out_view("out", view.extent(0));

  typedef Kokkos::RangePolicy<execution_space> my_exec_space;
  Kokkos::parallel_for("KokkosKernels::Common::Test::CheckFFS", my_exec_space(0, view.extent(0)),
                       ffscheck<view_type>(view, out_view));
  Kokkos::fence();
  return out_view;
}

}  // namespace Test

template <typename lno_t, typename device>
void test_set_bit_count() {
  const int array_size = 1000000;
  typedef Kokkos::View<lno_t *, device> myview;
  typedef typename myview::non_const_type nonconstview;

  nonconstview count_bit_view("count_bit_view", array_size);

  typename nonconstview::HostMirror hview = Kokkos::create_mirror_view(count_bit_view);

  for (int i = 0; i < array_size; ++i) {
    hview(i) = lno_t(rand()) * lno_t(rand());
  }

  Kokkos::deep_copy(count_bit_view, hview);

  // KokkosKernels::Impl::kk_print_1Dview(count_bit_view);

  myview out1 = Test::get_array_bit_count<myview, typename device::execution_space>(count_bit_view);
  myview out2 = Test::check_array_bit_count<myview, typename device::execution_space>(count_bit_view);
  // KokkosKernels::Impl::kk_print_1Dview(out1);
  // KokkosKernels::Impl::kk_print_1Dview(out2);

  bool is_identical = KokkosKernels::Impl::kk_is_identical_view<myview, myview, typename myview::value_type,
                                                                typename device::execution_space>(out1, out2, 0);
  EXPECT_TRUE(is_identical);
}

template <typename lno_t, typename device>
void test_ffs() {
  const int array_size = 1000000;
  typedef Kokkos::View<lno_t *, device> myview;
  typedef typename myview::non_const_type nonconstview;

  nonconstview count_bit_view("count_bit_view", array_size);

  typename nonconstview::HostMirror hview = Kokkos::create_mirror_view(count_bit_view);

  for (int i = 0; i < array_size; ++i) {
    hview(i) = lno_t(rand()) * lno_t(rand());
  }

  Kokkos::deep_copy(count_bit_view, hview);

  // KokkosKernels::Impl::kk_print_1Dview(count_bit_view);

  myview out1 = Test::get_ffs<myview, typename device::execution_space>(count_bit_view);
  myview out2 = Test::check_ffs<myview, typename device::execution_space>(count_bit_view);
  // KokkosKernels::Impl::kk_print_1Dview(out1);
  // KokkosKernels::Impl::kk_print_1Dview(out2);

  bool is_identical = KokkosKernels::Impl::kk_is_identical_view<myview, myview, typename myview::value_type,
                                                                typename device::execution_space>(out1, out2, 0);
  EXPECT_TRUE(is_identical);
}

TEST_F(TestCategory, common_set_bit_count) {
  test_set_bit_count<int, TestDevice>();
  test_set_bit_count<unsigned int, TestDevice>();
  test_set_bit_count<const int, TestDevice>();
  test_set_bit_count<const unsigned int, TestDevice>();

  test_set_bit_count<long, TestDevice>();
  test_set_bit_count<unsigned long, TestDevice>();
  test_set_bit_count<const long, TestDevice>();
  test_set_bit_count<const unsigned long, TestDevice>();

  test_set_bit_count<long long, TestDevice>();
  test_set_bit_count<unsigned long long, TestDevice>();
  test_set_bit_count<const long long, TestDevice>();
  test_set_bit_count<const unsigned long long, TestDevice>();
}

TEST_F(TestCategory, common_ffs) {
  test_ffs<int, TestDevice>();
  test_ffs<unsigned int, TestDevice>();
  test_ffs<const int, TestDevice>();
  test_ffs<const unsigned int, TestDevice>();

  test_ffs<long, TestDevice>();
  test_ffs<unsigned long, TestDevice>();
  test_ffs<const long, TestDevice>();
  test_ffs<const unsigned long, TestDevice>();

  test_ffs<long long, TestDevice>();
  test_ffs<unsigned long long, TestDevice>();
  test_ffs<const long long, TestDevice>();
  test_ffs<const unsigned long long, TestDevice>();
}
