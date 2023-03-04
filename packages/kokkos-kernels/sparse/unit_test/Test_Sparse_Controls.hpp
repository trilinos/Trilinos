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

#ifndef TEST_SPARSE_CONTROLS_HPP
#define TEST_SPARSE_CONTROLS_HPP

#include "KokkosKernels_Controls.hpp"

void test_controls_empty() {
  KokkosKernels::Experimental::Controls c;
  EXPECT_EQ(c.isParameter(""), false);
  EXPECT_EQ(c.getParameter(""), "");
  EXPECT_EQ(c.getParameter("", "default"), "default");
}

void test_controls_set() {
  KokkosKernels::Experimental::Controls c;
  c.setParameter("key", "value");
  EXPECT_EQ(c.isParameter("key"), true);
  EXPECT_EQ(c.getParameter("key"), "value");
  EXPECT_EQ(c.getParameter("key", "default"), "value");

  EXPECT_EQ(c.isParameter(""), false);
  EXPECT_EQ(c.getParameter(""), "");
  EXPECT_EQ(c.getParameter("", "default"), "default");
}

TEST_F(TestCategory, controls_empty) { test_controls_empty(); }
TEST_F(TestCategory, controls_set) { test_controls_set(); }

#endif  // TEST_SPARSE_CONTROLS_HPP
