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

#ifndef TEST_COMMON_IOTA_HPP
#define TEST_COMMON_IOTA_HPP

#include <type_traits>

#include "KokkosKernels_Iota.hpp"

template <typename T, typename Size>
void test_iota_constructor() {
  // empty iota
  {
    Iota<T, Size> i;
    EXPECT_EQ(i.size(), 0);
  }

  // basic iota
  {
    Iota<T, Size> ten(10);
    EXPECT_EQ(ten.size(), 10);
    for (size_t i = 0; i < ten.size(); ++i) {
      EXPECT_EQ(ten(i), i);
    }
  }

  // iota with negative offset
  if constexpr (std::is_signed_v<T>) {
    Iota<T, Size> three(3, -7);
    EXPECT_EQ(three.size(), 3);
    for (size_t i = 0; i < three.size(); ++i) {
      EXPECT_EQ(three(i), T(i) - T(7));
    }
  }

  // iota with positive offset
  {
    Iota<T, Size> three(3, 2);
    EXPECT_EQ(three.size(), 3);
    for (size_t i = 0; i < three.size(); ++i) {
      EXPECT_EQ(three(i), i + 2);
    }
  }

  // negative sizes are capped at 0
  if constexpr (std::is_signed_v<Size>) {
    {
      Iota<T, Size> i(-7);
      EXPECT_EQ(i.size(), 0);
    }
    {
      Iota<T, Size> i(-1, 2);
      EXPECT_EQ(i.size(), 0);
    }
  }
}

template <typename T, typename Size>
void test_iota_rank() {
  EXPECT_EQ((Iota<T, Size>::rank), 1);
}

template <typename T, typename Size>
void test_iota_non_const_value_type() {
  static_assert(std::is_same_v<typename Iota<T, Size>::non_const_value_type, T>,
                "Iota's non-const value type should be same as non-const type provided");
  static_assert(std::is_same_v<typename Iota<const T, Size>::non_const_value_type, T>,
                "Iota's non-const value type should be same as non-const version of "
                "const type provided");
}

template <typename T, typename Size>
void test_iota_subview() {
  // get the 7th and 8th elements of an Iota
  Iota<T, Size> ten(10, 1);                    // 1..<11
  Iota<T, Size> sub(ten, Kokkos::pair{7, 9});  // 8, 9

  EXPECT_EQ(sub.size(), 2);
  EXPECT_EQ(sub(0), 8);
  EXPECT_EQ(sub(1), 9);
}

template <typename T, typename Size>
void test_is_iota() {
  static_assert(KokkosKernels::Impl::is_iota_v<Iota<T, Size>>, "Iota should be an Iota");
  static_assert(!KokkosKernels::Impl::is_iota_v<int>, "int should not be an Iota");
}

template <typename T, typename Size>
void test_iota() {
  test_is_iota<T, Size>();
  test_iota_constructor<T, Size>();
  test_iota_rank<T, Size>();
  test_iota_non_const_value_type<T, Size>();
  test_iota_subview<T, Size>();
}

TEST_F(TestCategory, common_iota) {
  test_iota<int32_t, int32_t>();
  test_iota<int64_t, int32_t>();
  test_iota<uint32_t, int32_t>();
  test_iota<uint64_t, int32_t>();

  test_iota<int32_t, int64_t>();
  test_iota<int64_t, int64_t>();
  test_iota<uint32_t, int64_t>();
  test_iota<uint64_t, int64_t>();

  test_iota<int32_t, uint32_t>();
  test_iota<int64_t, uint32_t>();
  test_iota<uint32_t, uint32_t>();
  test_iota<uint64_t, uint32_t>();

  test_iota<int32_t, uint64_t>();
  test_iota<int64_t, uint64_t>();
  test_iota<uint32_t, uint64_t>();
  test_iota<uint64_t, uint64_t>();
}

#endif  // TEST_COMMON_IOTA_HPP
