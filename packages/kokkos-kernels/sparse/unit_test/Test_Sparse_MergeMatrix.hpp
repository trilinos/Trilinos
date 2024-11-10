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

#ifndef TEST_COMMON_MERGE_MATRIX_HPP
#define TEST_COMMON_MERGE_MATRIX_HPP

#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include "KokkosKernels_Iota.hpp"
#include "KokkosSparse_merge_matrix.hpp"

namespace Test_Sparse_MergeMatrix {

template <typename View>
View from_std_vec(const std::string &label, const std::vector<typename View::non_const_value_type> &vec) {
  Kokkos::View<const typename View::value_type *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> uvec(vec.data(),
                                                                                                   vec.size());
  View result(label, uvec.size());
  Kokkos::deep_copy(result, uvec);
  return result;
}

template <typename MMD, typename View>
struct CopyMmdToView {
  CopyMmdToView(const View &dst, const MMD &src) : dst_(dst), src_(src) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(size_t i) const { dst_(i) = src_(i); }

 private:
  View dst_;
  MMD src_;
};

template <typename MMD>
void expect_mmd_entries(const MMD &mmd, const std::vector<typename MMD::non_const_value_type> &expected) {
  using execution_space = typename MMD::execution_space;
  using Policy          = Kokkos::RangePolicy<execution_space>;
  using View            = Kokkos::View<typename MMD::non_const_value_type *, execution_space>;

  // size is as expected
  EXPECT_EQ(mmd.size(), expected.size());

  // values are as expected
  View view("mmd-values", mmd.size());
  execution_space space;
  Kokkos::parallel_for(Policy(space, 0, mmd.size()), CopyMmdToView(view, mmd));
  auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), view);
  space.fence();
  for (size_t i = 0; i < host.size(); ++i) {
    EXPECT_EQ(host(i), expected[i]);
  }
}

/*! \brief merge-matrix of two empty views

    Matrix is 0x0.
    Only diagonal 0 exists, and it should be size 0.
*/
template <typename AEntry, typename BEntry, typename ExecSpace>
void view_view_empty_empty() {
  using AView = Kokkos::View<AEntry *, ExecSpace>;
  using BView = Kokkos::View<BEntry *, ExecSpace>;
  using MMD   = KokkosSparse::Impl::MergeMatrixDiagonal<AView, BView>;

  AView a("view-view-empty-empty-a", 0);
  BView b("view-view-empty-empty-b", 0);
  expect_mmd_entries(MMD(a, b, 0), {});
}

/*! \brief merge-matrix of one empty view

    Matrix is Nx0.
    N diagonals exist, all length 0
*/
template <typename AEntry, typename BEntry, typename ExecSpace>
void view_view_full_empty() {
  using AView = Kokkos::View<AEntry *, ExecSpace>;
  using BView = Kokkos::View<BEntry *, ExecSpace>;
  using MMD   = KokkosSparse::Impl::MergeMatrixDiagonal<AView, BView>;

  size_t aNonzero = 5;
  AView a("view-view-full-empty-a", aNonzero);
  BView b("view-view-full-empty-b", 0);

  for (size_t diagonal = 0; diagonal < a.size() + b.size() - 1; ++diagonal) {
    expect_mmd_entries(MMD(a, b, diagonal), {});
  }
}

/*! \brief merge-matrix of one empty view

    Matrix is 0xN.
    N diagonals exist, all length 0
*/
template <typename AEntry, typename BEntry, typename ExecSpace>
void view_view_empty_full() {
  using AView = Kokkos::View<AEntry *, ExecSpace>;
  using BView = Kokkos::View<BEntry *, ExecSpace>;
  using MMD   = KokkosSparse::Impl::MergeMatrixDiagonal<AView, BView>;

  AView a("view-view-empty-full-a", 0);
  BView b = from_std_vec<BView>("view-view-empty-full-b", {0, 1, 2, 3});

  for (size_t diagonal = 0; diagonal < a.size() + b.size() - 1; ++diagonal) {
    expect_mmd_entries(MMD(a, b, diagonal), {});
  }
}

template <typename AView, typename BView>
std::tuple<AView, BView> view_view_case_all_zero() {
  // M[i,j] = 1 iff A[i] > B[j]
  //   B 0 1 2 3
  // A   -------
  // 0 | 0 0 0 0
  // 0 | 0 0 0 0
  // 0 | 0 0 0 0
  // 0 | 0 0 0 0
  AView a = from_std_vec<AView>("view-view-case-all-zero-a", {0, 0, 0, 0});
  BView b = from_std_vec<BView>("view-view-case-all-zero-b", {0, 1, 2, 3});

  return std::make_tuple(a, b);
}

template <typename AView, typename BView>
std::tuple<AView, BView> view_view_case_all_one() {
  // M[i,j] = 1 iff A[i] > B[j]
  //   B 0 0 0 0
  // A   -------
  // 1 | 1 1 1 1
  // 2 | 1 1 1 1
  // 3 | 1 1 1 1
  // 4 | 1 1 1 1
  AView a = from_std_vec<AView>("view-view-case-all-one-a", {1, 2, 3, 4});
  BView b = from_std_vec<BView>("view-view-case-all-one-b", {0, 0, 0, 0});

  return std::make_tuple(a, b);
}

template <typename AView, typename BView>
std::tuple<AView, BView> view_view_case_1() {
  // M[i,j] = 1 iff A[i] > B[j]
  //   B 0 1 2 3
  // A   -------
  // 1 | 1 0 0 0
  // 2 | 1 1 0 0
  // 3 | 1 1 1 0
  // 4 | 1 1 1 1
  AView a = from_std_vec<AView>("view-view-case-1-a", {1, 2, 3, 4});
  BView b = from_std_vec<BView>("view-view-case-1-b", {0, 1, 2, 3});

  // diagonal 0: {}
  // 1: {1}
  // 2: {1,0}
  // 3: {1,1,0}
  // 4: {1,1,0,0}
  // 5: {1,1,0}
  // 6: {1,0}
  // 7: {1}

  return std::make_tuple(a, b);
}

template <typename AView, typename BView>
std::tuple<AView, BView> view_view_case_2() {
  // M[i,j] = 1 iff A[i] > B[j]
  //   B 0 2 2 8 8 8
  // A   -----------
  // 1 | 1 0 0 0 0 0
  // 2 | 1 0 0 0 0 0
  // 9 | 1 1 1 1 1 1
  AView a = from_std_vec<AView>("view-view-case-2-a", {1, 2, 9});
  BView b = from_std_vec<BView>("view-view-case-2-b", {0, 2, 2, 8, 8, 8});
  // 0: {}
  // 1: {1}
  // 2: {1,0}
  // 3: {1,0,0}
  // 4: {1,0,0}
  // 5: {1,0,0}
  // 6: {1,0,0}
  // 7: {1,0}
  // 8: {1}
  return std::make_tuple(a, b);
}

template <typename AView, typename BView>
std::tuple<AView, BView> view_view_case_3() {
  using AEntry = typename AView::non_const_value_type;
  // M[i,j] = 1 iff A[i] > B[j]
  //    B 0 2 7
  // A    -----
  // -1 | 0 0 0
  //  9 | 1 1 1
  //  9 | 1 1 1
  AView a = from_std_vec<AView>("view-view-case-3-a", {AEntry(-1), AEntry(9), AEntry(9)});
  BView b = from_std_vec<BView>("view-view-case-3-b", {0, 2, 7});
  // 0: {}
  // 1: {0}
  // 2: {1,0}
  // 3: {1,1,0}
  // 4: {1,1}
  // 5: {1}
  return std::make_tuple(a, b);
}

template <typename AView, typename BView>
std::tuple<AView, BView> view_view_case_4() {
  using BEntry = typename BView::non_const_value_type;

  // M[i,j] = 1 iff A[i] > B[j]
  //   B -3 -1 7
  // A   -------
  // 1 | 1   1 0
  // 6 | 1   1 0
  // 6 | 1   1 0
  AView a = from_std_vec<AView>("view-view-case-4-a", {1, 6, 6});
  BView b = from_std_vec<BView>("view-view-case-4-b", {BEntry(-3), BEntry(-1), 7});
  // 0: {}
  // 1: {1}
  // 2: {1,1}
  // 3: {1,1,0}
  // 4: {1,0}
  // 5: {0}
  return std::make_tuple(a, b);
}

template <typename AView, typename BView>
std::tuple<AView, BView> view_view_case_5() {
  using AEntry = typename AView::non_const_value_type;
  using BEntry = typename BView::non_const_value_type;

  // M[i,j] = 1 iff A[i] > B[j]
  //    B -2 0 1
  // A    -------
  //  -3 | 0 0 0
  //  -2 | 0 0 0
  //  2  | 1 1 1
  AView a = from_std_vec<AView>("view-view-case-5-a", {AEntry{-3}, AEntry{-2}, AEntry{2}});
  BView b = from_std_vec<BView>("view-view-case-5-b", {BEntry{-2}, BEntry{0}, BEntry{1}});
  // 0: {}
  // 1: {0}
  // 2: {0,0}
  // 3: {0,0,0}
  // 4: {1,0}
  // 5: {1}
  return std::make_tuple(a, b);
}

/*! \brief merge-matrix of two views

    Matrix is MxN.
    M+N-1 diagonals exist.
*/
template <typename AEntry, typename BEntry, typename ExecSpace>
void view_view_full_full() {
  using AView          = Kokkos::View<AEntry *, ExecSpace>;
  using BView          = Kokkos::View<BEntry *, ExecSpace>;
  using MMD            = KokkosSparse::Impl::MergeMatrixDiagonal<AView, BView>;
  using mmd_value_type = typename MMD::non_const_value_type;

  {
    auto [a, b] = view_view_case_all_zero<AView, BView>();
    for (size_t diagonal = 0; diagonal < a.size() + b.size() - 1; ++diagonal) {
      MMD mmd(a, b, diagonal);
      // every matrix entry on this diagonal is 0
      expect_mmd_entries(mmd, std::vector<mmd_value_type>(mmd.size(), mmd_value_type(0)));
    }
  }
  {
    auto [a, b] = view_view_case_all_one<AView, BView>();
    for (size_t diagonal = 0; diagonal < a.size() + b.size() - 1; ++diagonal) {
      MMD mmd(a, b, diagonal);
      // every matrix entry on this diagonal is 0
      expect_mmd_entries(mmd, std::vector<mmd_value_type>(mmd.size(), mmd_value_type(1)));
    }
  }
  {
    auto [a, b] = view_view_case_1<AView, BView>();
    expect_mmd_entries(MMD(a, b, 0), {});
    expect_mmd_entries(MMD(a, b, 1), {1});
    expect_mmd_entries(MMD(a, b, 2), {1, 0});
    expect_mmd_entries(MMD(a, b, 3), {1, 1, 0});
    expect_mmd_entries(MMD(a, b, 4), {1, 1, 0, 0});
    expect_mmd_entries(MMD(a, b, 5), {1, 1, 0});
    expect_mmd_entries(MMD(a, b, 6), {1, 0});
    expect_mmd_entries(MMD(a, b, 7), {1});
  }
  {
    auto [a, b] = view_view_case_2<AView, BView>();
    expect_mmd_entries(MMD(a, b, 0), {});
    expect_mmd_entries(MMD(a, b, 1), {1});
    expect_mmd_entries(MMD(a, b, 2), {1, 0});
    expect_mmd_entries(MMD(a, b, 3), {1, 0, 0});
    expect_mmd_entries(MMD(a, b, 4), {1, 0, 0});
    expect_mmd_entries(MMD(a, b, 5), {1, 0, 0});
    expect_mmd_entries(MMD(a, b, 6), {1, 0, 0});
    expect_mmd_entries(MMD(a, b, 7), {1, 0});
    expect_mmd_entries(MMD(a, b, 8), {1});
  }
  if constexpr (std::is_signed_v<AEntry>) {
    auto [a, b] = view_view_case_3<AView, BView>();
    expect_mmd_entries(MMD(a, b, 0), {});
    expect_mmd_entries(MMD(a, b, 1), {0});
    expect_mmd_entries(MMD(a, b, 2), {1, 0});
    expect_mmd_entries(MMD(a, b, 3), {1, 1, 0});
    expect_mmd_entries(MMD(a, b, 4), {1, 1});
    expect_mmd_entries(MMD(a, b, 5), {1});
  }
  if constexpr (std::is_signed_v<BEntry>) {
    auto [a, b] = view_view_case_4<AView, BView>();
    expect_mmd_entries(MMD(a, b, 0), {});
    expect_mmd_entries(MMD(a, b, 1), {1});
    expect_mmd_entries(MMD(a, b, 2), {1, 1});
    expect_mmd_entries(MMD(a, b, 3), {1, 1, 0});
    expect_mmd_entries(MMD(a, b, 4), {1, 0});
    expect_mmd_entries(MMD(a, b, 5), {0});
  }
  if constexpr (std::is_signed_v<AEntry> && std::is_signed_v<BEntry>) {
    auto [a, b] = view_view_case_5<AView, BView>();
    expect_mmd_entries(MMD(a, b, 0), {});
    expect_mmd_entries(MMD(a, b, 1), {0});
    expect_mmd_entries(MMD(a, b, 2), {0, 0});
    expect_mmd_entries(MMD(a, b, 3), {1, 0, 0});
    expect_mmd_entries(MMD(a, b, 4), {1, 0});
    expect_mmd_entries(MMD(a, b, 5), {1});
  }
}

template <typename AEntry, typename BEntry, typename ExecSpace>
void test_view_view() {
  view_view_empty_empty<AEntry, BEntry, ExecSpace>();
  view_view_full_empty<AEntry, BEntry, ExecSpace>();
  view_view_empty_full<AEntry, BEntry, ExecSpace>();
  view_view_full_full<AEntry, BEntry, ExecSpace>();
}

/*! \brief merge-matrix of an empty view and empty iota

    Matrix is 0x0.
    Only diagonal 0 exists, and it should be size 0.
*/
template <typename AEntry, typename BEntry, typename ExecSpace>
void view_iota_empty_empty() {
  using AView = Kokkos::View<AEntry *, ExecSpace>;
  using BView = KokkosKernels::Impl::Iota<BEntry>;
  using MMD   = KokkosSparse::Impl::MergeMatrixDiagonal<AView, BView>;

  AView a("view-iota-empty-empty-a", 0);
  BView b(0);
  EXPECT_EQ(MMD(a, b, 0).size(), 0);
}

/*! \brief merge-matrix of a full view and empty iota

    Matrix is Nx0.
    N diagonals exist, all length 0
*/
template <typename AEntry, typename BEntry, typename ExecSpace>
void view_iota_full_empty() {
  using AView = Kokkos::View<AEntry *, ExecSpace>;
  using BView = KokkosKernels::Impl::Iota<BEntry>;
  using MMD   = KokkosSparse::Impl::MergeMatrixDiagonal<AView, BView>;

  size_t aNonzero = 5;
  AView a("view-iota-full-empty-a", aNonzero);
  BView b(0);

  for (size_t diagonal = 0; diagonal < a.size() + b.size() - 1; ++diagonal) {
    EXPECT_EQ(MMD(a, b, diagonal).size(), 0);
  }
}

/*! \brief merge-matrix of and empty view and a full iota

    Matrix is 0xN.
    N diagonals exist, all length 0
*/
template <typename AEntry, typename BEntry, typename ExecSpace>
void view_iota_empty_full() {
  using AView = Kokkos::View<AEntry *, ExecSpace>;
  using BView = KokkosKernels::Impl::Iota<BEntry>;
  using MMD   = KokkosSparse::Impl::MergeMatrixDiagonal<AView, BView>;

  AView a("view-iota-empty-full-a", 0);
  BView b(4);

  for (size_t diagonal = 0; diagonal < a.size() + b.size() - 1; ++diagonal) {
    EXPECT_EQ(MMD(a, b, diagonal).size(), 0);
  }
}

template <typename AView, typename BView>
std::tuple<AView, BView> view_iota_case_all_zero() {
  // M[i,j] = 1 iff A[i] > B[j]
  //   B 0 1 2 3
  // A   -------
  // 0 | 0 0 0 0
  // 0 | 0 0 0 0
  // 0 | 0 0 0 0
  // 0 | 0 0 0 0
  AView a = from_std_vec<AView>("view-iota-case-all-zero-a", {0, 0, 0, 0});
  BView b(4);

  return std::make_tuple(a, b);
}

template <typename AView, typename BView>
std::tuple<AView, BView> view_iota_case_all_one() {
  // M[i,j] = 1 iff A[i] > B[j]
  //   B 0 1 2 3
  // A   -------
  // 5 | 1 1 1 1
  // 6 | 1 1 1 1
  // 7 | 1 1 1 1
  // 8 | 1 1 1 1
  AView a = from_std_vec<AView>("view-iota-case-all-one-a", {5, 6, 7, 8});
  BView b(4);

  return std::make_tuple(a, b);
}

template <typename AView, typename BView>
std::tuple<AView, BView> view_iota_case_1() {
  // M[i,j] = 1 iff A[i] > B[j]
  //   B 0 1 2 3
  // A   -------
  // 1 | 1 0 0 0
  // 2 | 1 1 0 0
  // 3 | 1 1 1 0
  // 4 | 1 1 1 1
  AView a = from_std_vec<AView>("view-iota-case-1-a", {1, 2, 3, 4});
  BView b(4);

  // diagonal 0: {}
  // 1: {1}
  // 2: {1,0}
  // 3: {1,1,0}
  // 4: {1,1,0,0}
  // 5: {1,1,0}
  // 6: {1,0}
  // 7: {1}

  return std::make_tuple(a, b);
}

/*! \brief merge-matrix of a full view with a full iota

    Matrix is MxN.
    M+N-1 diagonals exist.
*/
template <typename AEntry, typename BEntry, typename ExecSpace>
void view_iota_full_full() {
  using AView          = Kokkos::View<AEntry *, ExecSpace>;
  using BView          = KokkosKernels::Impl::Iota<BEntry>;
  using MMD            = KokkosSparse::Impl::MergeMatrixDiagonal<AView, BView>;
  using mmd_value_type = typename MMD::non_const_value_type;

  {
    auto [a, b] = view_iota_case_all_zero<AView, BView>();
    for (size_t diagonal = 0; diagonal < a.size() + b.size() - 1; ++diagonal) {
      MMD mmd(a, b, diagonal);
      // every matrix entry on this diagonal is 0
      expect_mmd_entries(mmd, std::vector<mmd_value_type>(mmd.size(), mmd_value_type(0)));
    }
  }
  {
    auto [a, b] = view_iota_case_all_one<AView, BView>();
    for (size_t diagonal = 0; diagonal < a.size() + b.size() - 1; ++diagonal) {
      MMD mmd(a, b, diagonal);
      // every matrix entry on this diagonal is 1
      expect_mmd_entries(mmd, std::vector<mmd_value_type>(mmd.size(), mmd_value_type(1)));
    }
  }
  {
    auto [a, b] = view_iota_case_1<AView, BView>();
    expect_mmd_entries(MMD(a, b, 0), {});
    expect_mmd_entries(MMD(a, b, 1), {1});
    expect_mmd_entries(MMD(a, b, 2), {1, 0});
    expect_mmd_entries(MMD(a, b, 3), {1, 1, 0});
    expect_mmd_entries(MMD(a, b, 4), {1, 1, 0, 0});
    expect_mmd_entries(MMD(a, b, 5), {1, 1, 0});
    expect_mmd_entries(MMD(a, b, 6), {1, 0});
    expect_mmd_entries(MMD(a, b, 7), {1});
  }
}

template <typename AEntry, typename BEntry, typename ExecSpace>
void test_view_iota() {
  view_iota_empty_empty<AEntry, BEntry, ExecSpace>();
  view_iota_full_empty<AEntry, BEntry, ExecSpace>();
  view_iota_empty_full<AEntry, BEntry, ExecSpace>();
  view_iota_full_full<AEntry, BEntry, ExecSpace>();
}

template <typename AEntry, typename BEntry, typename ExecSpace>
void test_rank() {
  {
    using AView = Kokkos::View<AEntry *, ExecSpace>;
    using BView = Kokkos::View<BEntry *, ExecSpace>;
    using MMD   = KokkosSparse::Impl::MergeMatrixDiagonal<AView, BView>;
    static_assert(MMD::rank == 1, "MergeMatrixDiagonal should look like a rank-1 view");
  }

  {
    using AView = Kokkos::View<AEntry *, ExecSpace>;
    using BView = KokkosKernels::Impl::Iota<BEntry>;
    using MMD   = KokkosSparse::Impl::MergeMatrixDiagonal<AView, BView>;
    static_assert(MMD::rank == 1, "MergeMatrixDiagonal should look like a rank-1 view");
  }
}

template <typename AEntry, typename BEntry, typename ExecSpace>
void test_merge_matrix() {
  test_rank<AEntry, BEntry, ExecSpace>();
  test_view_view<AEntry, BEntry, ExecSpace>();
  test_view_iota<AEntry, BEntry, ExecSpace>();
}

}  // namespace Test_Sparse_MergeMatrix

TEST_F(TestCategory, common_merge_matrix) {
  // clang-format off
  Test_Sparse_MergeMatrix::test_merge_matrix<int32_t, int32_t, TestDevice>();
  Test_Sparse_MergeMatrix::test_merge_matrix<int64_t, int32_t, TestDevice>();
  Test_Sparse_MergeMatrix::test_merge_matrix<uint32_t, int32_t, TestDevice>();

  Test_Sparse_MergeMatrix::test_merge_matrix<int32_t, int64_t, TestDevice>();
  Test_Sparse_MergeMatrix::test_merge_matrix<int64_t, int64_t, TestDevice>();
  Test_Sparse_MergeMatrix::test_merge_matrix<uint32_t, int64_t, TestDevice>();

  Test_Sparse_MergeMatrix::test_merge_matrix<int32_t, uint32_t, TestDevice>();
  Test_Sparse_MergeMatrix::test_merge_matrix<int64_t, uint32_t, TestDevice>();
  Test_Sparse_MergeMatrix::test_merge_matrix<uint32_t, uint32_t, TestDevice>();
  Test_Sparse_MergeMatrix::test_merge_matrix<uint64_t, uint32_t, TestDevice>();

  Test_Sparse_MergeMatrix::test_merge_matrix<uint32_t, uint64_t, TestDevice>();
  Test_Sparse_MergeMatrix::test_merge_matrix<uint64_t, uint64_t, TestDevice>();

  Test_Sparse_MergeMatrix::test_merge_matrix<float, float, TestDevice>();
  Test_Sparse_MergeMatrix::test_merge_matrix<float, double, TestDevice>();
  Test_Sparse_MergeMatrix::test_merge_matrix<double, float, TestDevice>();
  Test_Sparse_MergeMatrix::test_merge_matrix<double, double, TestDevice>();

  // test some select integer / float combos
  Test_Sparse_MergeMatrix::test_merge_matrix<float, uint64_t, TestDevice>();
  Test_Sparse_MergeMatrix::test_merge_matrix<int32_t, double, TestDevice>();

  // no generally safe way to compare all possible values of these types
  // Test_Sparse_MergeMatrix::test_merge_matrix<int64_t, uint64_t, TestDevice>();
  // Test_Sparse_MergeMatrix::test_merge_matrix<int32_t, uint64_t, TestDevice>();
  // Test_Sparse_MergeMatrix::test_merge_matrix<uint64_t, int32_t, TestDevice>();

  // clang-format on
}

#endif  // TEST_COMMON_MERGE_MATRIX_HPP
