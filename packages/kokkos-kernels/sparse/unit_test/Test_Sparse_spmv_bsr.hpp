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

/*! \file Test_Sparse_spmv_bsr.hpp

  Test the following 256 combos for at least a few matcies.

  Algorithms              Alpha     Beta     Block Sizes    Modes
  (none)                  0         0        1              N
  native              x   1      x  1     x  2           x  T
  experimental_bsr_tc     -1        -1       5              C
                          3.7       -1.5     9              H

  There are also a subset of tests on larger matrices

  Multivector products are also tested for these cases with 1 and 7 vectors
*/

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include <KokkosKernels_TestUtils.hpp>
#include <KokkosKernels_Test_Structured_Matrix.hpp>
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosKernels_Utils.hpp>
#include "KokkosKernels_default_types.hpp"
#include <KokkosKernels_NaN.hpp>

#include "KokkosSparse_spmv.hpp"
#include "KokkosSparse_BsrMatrix.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_crs_to_bsr_impl.hpp"
#include "KokkosSparse_bsr_to_crs_impl.hpp"
#include "KokkosSparse_Utils.hpp"

using kokkos_complex_double = Kokkos::complex<double>;
using kokkos_complex_float  = Kokkos::complex<float>;

namespace Test_Spmv_Bsr {

/*! \brief Maximum value used to fill A */
template <typename T>
constexpr T max_a() {
  T discard, maxVal;
  KokkosKernels::Impl::getRandomBounds(10.0, discard, maxVal);
  return maxVal;
}

/*! \brief Maximum value used to fill X */
template <typename T>
constexpr T max_x() {
  T discard, maxVal;
  KokkosKernels::Impl::getRandomBounds(10.0, discard, maxVal);
  return maxVal;
}

/*! \brief Maximum value used to fill Y */
template <typename T>
constexpr T max_y() {
  T discard, maxVal;
  KokkosKernels::Impl::getRandomBounds(10.0, discard, maxVal);
  return maxVal;
}

/*! \brief whether the mode transposes the matrix*/
inline bool mode_is_transpose(const char *mode) { return mode[0] == 'T' || mode[0] == 'H'; }

/*! \brief Get the max nonzeros (not max nonzero _blocks_) per row of Op(A) */
template <typename Bsr>
inline size_t opMaxNnzPerRow(const Bsr &A, bool trans) {
  if (trans) {
    auto At = KokkosSparse::Impl::transpose_bsr_matrix(A);
    return At.blockDim() *
           (size_t)KokkosSparse::Impl::graph_max_degree<typename Bsr::execution_space, typename Bsr::ordinal_type>(
               At.graph.row_map);
  } else {
    return A.blockDim() *
           (size_t)KokkosSparse::Impl::graph_max_degree<typename Bsr::execution_space, typename Bsr::ordinal_type>(
               A.graph.row_map);
  }
}

/*! \brief 0x0 matrix */
template <typename Bsr>
Bsr bsr_corner_case_0_by_0(const int blockSize) {
  return Bsr("empty", 0, 0, 0, nullptr, nullptr, nullptr, blockSize);
}

/*! \brief 0x1 matrix */
template <typename Bsr>
Bsr bsr_corner_case_0_by_1(const int blockSize) {
  return Bsr("empty", 0, blockSize, 0, nullptr, nullptr, nullptr, blockSize);
}

/*! \brief 1x0 matrix */
template <typename Bsr>
Bsr bsr_corner_case_1_by_0(const int blockSize) {
  return Bsr("empty", blockSize, 0, 0, nullptr, nullptr, nullptr, blockSize);
}

template <typename Bsr>
Bsr bsr_random(const int blockSize, const int blockRows, const int blockCols) {
  using scalar_type  = typename Bsr::non_const_value_type;
  using ordinal_type = typename Bsr::non_const_ordinal_type;
  using size_type    = typename Bsr::non_const_size_type;
  using Crs          = KokkosSparse::CrsMatrix<scalar_type, ordinal_type, typename Bsr::device_type, void, size_type>;
  using Graph        = typename Crs::staticcrsgraph_type;

  // construct a random Crs Matrix
  Test::RandCsMatrix<scalar_type, Kokkos::LayoutLeft, typename Bsr::device_type, ordinal_type, size_type> rcs(
      blockRows, blockCols, scalar_type(0), max_a<scalar_type>());

  const auto colids = Kokkos::subview(rcs.get_ids(), Kokkos::make_pair(size_type(0), rcs.get_nnz()));
  const auto vals   = Kokkos::subview(rcs.get_vals(), Kokkos::make_pair(size_type(0), rcs.get_nnz()));
  Graph graph(colids, rcs.get_map());
  Crs crs("crs", blockCols, vals, graph);

  // expand to Bsr matrix
  return KokkosSparse::Impl::expand_crs_to_bsr<Bsr>(crs, blockSize);
}

/*! \brief test a specific spmv

*/
template <typename Handle, typename Bsr, typename Crs, typename XVector, typename YVector,
          typename Alpha = typename Bsr::non_const_value_type, typename Beta = typename Bsr::non_const_value_type>
void test_spmv(Handle *handle, const char *mode, const Alpha &alpha, const Beta &beta, const Bsr &a, const Crs &acrs,
               size_t maxNnzPerRow, const XVector &x, const YVector &y) {
  using scalar_type  = typename Bsr::non_const_value_type;
  using ordinal_type = typename Bsr::non_const_ordinal_type;
  using KATS         = Kokkos::ArithTraits<scalar_type>;
  using mag_type     = typename KATS::mag_type;

  // generate expected result from reference (CRS) implementation
  YVector yExp("yExp", y.extent(0));
  Kokkos::deep_copy(yExp, y);
  KokkosSparse::spmv(mode, alpha, acrs, x, beta, yExp);

  // scratch space for actual value (don't modify input)
  YVector yAct("yAct", y.extent(0));
  Kokkos::deep_copy(yAct, y);

  KokkosSparse::spmv(handle, mode, alpha, a, x, beta, yAct);

  // compare yExp and yAct
  auto hyExp = Kokkos::create_mirror_view(yExp);
  auto hyAct = Kokkos::create_mirror_view(yAct);
  Kokkos::deep_copy(hyExp, yExp);
  Kokkos::deep_copy(hyAct, yAct);

  /* assume that any floating-point op may introduce eps() error
     scaling y is one op
     dot product of x is two ops per entry (mul and add)

     10x means same order of magnitude
  */
  const mag_type tolerance = KATS::eps() * KATS::abs(beta) * KATS::abs(max_y<scalar_type>()) +
                             10 * KATS::eps() * maxNnzPerRow * KATS::abs(alpha) * KATS::abs(max_a<scalar_type>()) *
                                 KATS::abs(max_x<scalar_type>());

  std::vector<ordinal_type> errIdx;

  for (ordinal_type i = 0; i < ordinal_type(hyAct.extent(0)); ++i) {
    if (KATS::abs(hyExp(i) - hyAct(i)) > tolerance) {
      errIdx.push_back(i);
    }
  }

  if (!errIdx.empty()) {
    std::string alg = KokkosSparse::get_spmv_algorithm_name(handle->get_algorithm());

    std::cerr << __FILE__ << ":" << __LINE__ << " BsrMatrix SpMV failure!" << std::endl;
    std::cerr << "alg:          " << alg << std::endl;
    std::cerr << "mode:         " << mode << std::endl;
    std::cerr << "A:            " << a.numRows() << "x" << a.numCols() << std::endl;
    std::cerr << "A blockdim:   " << a.blockDim() << std::endl;
    std::cerr << "alpha:        " << alpha << std::endl;
    std::cerr << "beta:         " << beta << std::endl;
    std::cerr << "maxNnzPerRow: " << maxNnzPerRow << std::endl;
    std::cerr << "First 100 errors:" << std::endl;
    std::cerr << "y\texp\tact\terr\ttol" << std::endl;
    std::cerr << "-\t---\t---\t---\t---" << std::endl;
    for (size_t i = 0; i < 100 && i < errIdx.size(); ++i) {
      size_t ei = errIdx[i];
      // clang-format off
      std::cerr << ei 
                << "\t" << hyExp(ei)
                << "\t" << hyAct(ei)
                << "\t" << KATS::abs(hyExp(ei) - hyAct(ei))
                << "\t" << tolerance
                << std::endl;
      // clang-format on
    }
  }

  EXPECT_TRUE(errIdx.empty());
}

template <typename Bsr>
struct VectorTypeFor {
  using type = Kokkos::View<typename Bsr::non_const_value_type *, typename Bsr::device_type>;
};

template <typename Bsr>
std::tuple<Bsr, typename VectorTypeFor<Bsr>::type, typename VectorTypeFor<Bsr>::type> spmv_corner_case_0_by_0(
    const char * /*mode*/, const int blockSize) {
  using vector_type = typename VectorTypeFor<Bsr>::type;
  Bsr a             = bsr_corner_case_0_by_0<Bsr>(blockSize);
  vector_type x("x", 0);
  vector_type y("y", 0);
  return std::make_tuple(a, x, y);
}

template <typename Bsr>
std::tuple<Bsr, typename VectorTypeFor<Bsr>::type, typename VectorTypeFor<Bsr>::type> spmv_corner_case_0_by_1(
    const char *mode, const int blockSize) {
  using vector_type     = typename VectorTypeFor<Bsr>::type;
  using execution_space = typename Bsr::execution_space;
  using scalar_type     = typename Bsr::non_const_value_type;
  Bsr a                 = bsr_corner_case_0_by_1<Bsr>(blockSize);

  size_t nx = a.numCols() * a.blockDim();
  size_t ny = a.numRows() * a.blockDim();
  if (mode_is_transpose(mode)) {
    std::swap(nx, ny);
  }
  vector_type x("x", nx);
  vector_type y("y", ny);

  Kokkos::Random_XorShift64_Pool<execution_space> random(13718);
  Kokkos::fill_random(x, random, max_x<scalar_type>());
  Kokkos::fill_random(y, random, max_y<scalar_type>());

  return std::make_tuple(a, x, y);
}

template <typename Bsr>
std::tuple<Bsr, typename VectorTypeFor<Bsr>::type, typename VectorTypeFor<Bsr>::type> spmv_corner_case_1_by_0(
    const char *mode, const int blockSize) {
  using vector_type     = typename VectorTypeFor<Bsr>::type;
  using execution_space = typename Bsr::execution_space;
  using scalar_type     = typename Bsr::non_const_value_type;
  Bsr a                 = bsr_corner_case_1_by_0<Bsr>(blockSize);

  size_t nx = a.numCols() * a.blockDim();
  size_t ny = a.numRows() * a.blockDim();
  if (mode_is_transpose(mode)) {
    std::swap(nx, ny);
  }
  vector_type x("x", nx);
  vector_type y("y", ny);

  Kokkos::Random_XorShift64_Pool<execution_space> random(13718);
  Kokkos::fill_random(x, random, max_x<scalar_type>());
  Kokkos::fill_random(y, random, max_y<scalar_type>());

  return std::make_tuple(a, x, y);
}

/*! \brief

*/
template <typename Bsr>
std::tuple<Bsr, typename VectorTypeFor<Bsr>::type, typename VectorTypeFor<Bsr>::type> spmv_random(const char *mode,
                                                                                                  const int blockSize,
                                                                                                  const int blockRows,
                                                                                                  const int blockCols) {
  using scalar_type = typename Bsr::non_const_value_type;

  // expand to Bsr matrix
  Bsr a = bsr_random<Bsr>(blockSize, blockRows, blockCols);

  // generate some random vectors
  using vector_type     = typename VectorTypeFor<Bsr>::type;
  using execution_space = typename Bsr::execution_space;

  size_t nx = a.numCols() * a.blockDim();
  size_t ny = a.numRows() * a.blockDim();
  if (mode_is_transpose(mode)) {
    std::swap(nx, ny);
  }
  vector_type x("x", nx);
  vector_type y("y", ny);

  Kokkos::Random_XorShift64_Pool<execution_space> random(13718);
  Kokkos::fill_random(x, random, max_x<scalar_type>());
  Kokkos::fill_random(y, random, max_y<scalar_type>());

  return std::make_tuple(a, x, y);
}

/*! \brief create random x and y multivectors for a given matrix and spmv mode
 */
template <typename Bsr>
auto random_vecs_for_spmv(const char *mode, const Bsr &a, const bool nans = false)
    -> std::tuple<typename VectorTypeFor<Bsr>::type, typename VectorTypeFor<Bsr>::type> {
  using scalar_type     = typename Bsr::non_const_value_type;
  using vector_type     = typename VectorTypeFor<Bsr>::type;
  using execution_space = typename Bsr::execution_space;
  using policy_type     = Kokkos::RangePolicy<typename vector_type::execution_space>;

  size_t nx = a.numCols() * a.blockDim();
  size_t ny = a.numRows() * a.blockDim();
  if (mode_is_transpose(mode)) {
    std::swap(nx, ny);
  }
  vector_type x("x", nx);
  vector_type y("y", ny);

  Kokkos::Random_XorShift64_Pool<execution_space> random(13718);
  Kokkos::fill_random(x, random, max_x<scalar_type>());
  Kokkos::fill_random(y, random, max_y<scalar_type>());

  if (nans) {
    Kokkos::parallel_for(
        policy_type(0, x.extent(0)), KOKKOS_LAMBDA(size_t i) {
          if (0 == (i % 17)) {
            x(i) = KokkosKernels::Impl::quiet_NaN<scalar_type>();
          }
        });
    Kokkos::parallel_for(
        policy_type(0, y.extent(0)), KOKKOS_LAMBDA(size_t i) {
          if (0 == (i % 17)) {
            y(i) = KokkosKernels::Impl::quiet_NaN<scalar_type>();
          }
        });
  }

  return std::make_tuple(x, y);
}

/*! \brief test all combos of the provided matrix
 */
template <typename Bsr, typename Crs>
void test_spmv_combos(const char *mode, const Bsr &a, const Crs &acrs, size_t maxNnzPerRow) {
  using namespace KokkosSparse;
  using scalar_type     = typename Bsr::non_const_value_type;
  using execution_space = typename Bsr::execution_space;

  auto [x, y]                     = random_vecs_for_spmv(mode, a);
  auto [x_with_nans, y_with_nans] = random_vecs_for_spmv(mode, a, true);

  using handle_t = SPMVHandle<execution_space, Bsr, decltype(x), decltype(y)>;

  // cover a variety of algorithms
  std::vector<std::unique_ptr<handle_t>> handles;
  for (SPMVAlgorithm algo : {SPMV_DEFAULT, SPMV_NATIVE, SPMV_BSR_V41})
    handles.push_back(std::make_unique<handle_t>(algo));

  // Tensor core algorithm temporarily disabled, fails on V100
  /*
  if constexpr (KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()) {
#if defined(KOKKOS_ENABLE_CUDA)
    if constexpr (std::is_same_v<execution_space, Kokkos::Cuda>) {
#if defined(KOKKOS_ARCH_AMPERE) || defined(KOKKOS_ARCH_VOLTA)
      handles.push_back(new handle_t(SPMV_BSR_TC));
#if defined(KOKKOS_ARCH_AMPERE)
      // Also call SPMV_BSR_TC with Precision = Double on Ampere
      handles.push_back(new handle_t(SPMV_BSR_TC));
      handles.back()->bsr_tc_precision = Experimental::Bsr_TC_Precision::Double;
#endif  // AMPERE
#endif  // AMPERE || VOLTA
    }
#endif  // CUDA
  }
  */

  for (std::unique_ptr<handle_t> &handle : handles) {
    for (scalar_type alpha : {scalar_type(0), scalar_type(1), scalar_type(-1), scalar_type(3.7)}) {
      for (scalar_type beta : {scalar_type(0), scalar_type(1), scalar_type(-1), scalar_type(-1.5)}) {
        test_spmv(handle.get(), mode, alpha, beta, a, acrs, maxNnzPerRow, x, y);
        if (beta == scalar_type(0)) {
          test_spmv(handle.get(), mode, alpha, beta, a, acrs, maxNnzPerRow, x_with_nans, y_with_nans);
        }
      }
    }
  }
}

/*! \brief test all combos of all matrices with different block sizes
 */
template <typename Scalar, typename Ordinal, typename Offset, typename Device>
void test_spmv_corner_cases() {
  using Bsr = KokkosSparse::Experimental::BsrMatrix<Scalar, Ordinal, Device, void, Offset>;
  using Crs = KokkosSparse::CrsMatrix<Scalar, Ordinal, Device, void, Offset>;
  for (auto mode : {"N", "T", "C", "H"}) {
    for (int bs : {1, 2, 5, 9}) {
      {
        auto A    = bsr_corner_case_0_by_0<Bsr>(bs);
        auto Acrs = KokkosSparse::Impl::bsr_to_crs<Crs>(A);
        test_spmv_combos(mode, A, Acrs, 0);
      }
      {
        auto A    = bsr_corner_case_0_by_1<Bsr>(bs);
        auto Acrs = KokkosSparse::Impl::bsr_to_crs<Crs>(A);
        test_spmv_combos(mode, A, Acrs, 0);
      }
      {
        auto A    = bsr_corner_case_1_by_0<Bsr>(bs);
        auto Acrs = KokkosSparse::Impl::bsr_to_crs<Crs>(A);
        test_spmv_combos(mode, A, Acrs, 0);
      }
    }
  }
}

template <typename Scalar, typename Ordinal, typename Offset, typename Device>
void test_spmv_random() {
  using Bsr = KokkosSparse::Experimental::BsrMatrix<Scalar, Ordinal, Device, void, Offset>;
  using Crs = KokkosSparse::CrsMatrix<Scalar, Ordinal, Device, void, Offset>;
  // thoroughly test smaller matrices
  std::vector<std::pair<int, int>> shapes = {{10, 10}, {10, 50}, {50, 10}};
  for (auto &shape : shapes) {
    for (int bs : {1, 2, 5, 9}) {
      auto A                   = bsr_random<Bsr>(bs, shape.first, shape.second);
      auto Acrs                = KokkosSparse::Impl::bsr_to_crs<Crs>(A);
      size_t maxNnzPerRow      = opMaxNnzPerRow(A, false);
      size_t maxNnzPerRowTrans = opMaxNnzPerRow(A, true);
      for (auto mode : {"N", "T", "C", "H"}) {
        test_spmv_combos(mode, A, Acrs, mode_is_transpose(mode) ? maxNnzPerRowTrans : maxNnzPerRow);
      }
    }
  }

  // test a tougher case on a big matrix
  {
    constexpr int blockSizePrime = 7;
    constexpr int smallPrime     = 11;
    constexpr int largePrime     = 499;
    auto A                       = bsr_random<Bsr>(blockSizePrime, smallPrime, largePrime);
    auto Acrs                    = KokkosSparse::Impl::bsr_to_crs<Crs>(A);
    size_t maxNnzPerRow          = opMaxNnzPerRow(A, false);
    size_t maxNnzPerRowTrans     = opMaxNnzPerRow(A, true);
    for (auto mode : {"N", "T"}) {
      test_spmv_combos(mode, A, Acrs, mode_is_transpose(mode) ? maxNnzPerRowTrans : maxNnzPerRow);
    }
  }
}

template <typename Scalar, typename Ordinal, typename Offset, typename Device>
void test_spmv() {
  test_spmv_corner_cases<Scalar, Ordinal, Offset, Device>();
  test_spmv_random<Scalar, Ordinal, Offset, Device>();
}

// ----------------------------------------------------------------------------
// Multivector
// ----------------------------------------------------------------------------

// Note: if mode_is_transpose(mode), then maxNnzPerRow is for A^T. Otherwise,
// it's for A.
template <typename Handle, typename Bsr, typename Crs, typename XVector, typename YVector, typename Alpha,
          typename Beta>
void test_spm_mv(Handle *handle, const char *mode, const Alpha &alpha, const Beta &beta, const Bsr &a, const Crs &acrs,
                 size_t maxNnzPerRow, const XVector &x, const YVector &y) {
  using scalar_type  = typename Bsr::non_const_value_type;
  using ordinal_type = typename Bsr::non_const_ordinal_type;
  using KATS         = Kokkos::ArithTraits<scalar_type>;
  using mag_type     = typename KATS::mag_type;

  // generate expected result from reference (CRS) implementation
  YVector yExp("yExp", y.extent(0), y.extent(1));
  Kokkos::deep_copy(yExp, y);
  KokkosSparse::spmv(mode, alpha, acrs, x, beta, yExp);

  // scratch space for actual value (don't modify input)
  YVector yAct("yAct", y.extent(0), y.extent(1));
  Kokkos::deep_copy(yAct, y);

  KokkosSparse::spmv(handle, mode, alpha, a, x, beta, yAct);

  // compare yExp and yAct
  auto hyExp = Kokkos::create_mirror_view(yExp);
  auto hyAct = Kokkos::create_mirror_view(yAct);
  Kokkos::deep_copy(hyExp, yExp);
  Kokkos::deep_copy(hyAct, yAct);

  /* assume that any floating-point op may introduce eps() error
     scaling y is one op
     dot product of x is two ops per entry (mul and add)
  */
  const mag_type tolerance = KATS::eps() * KATS::abs(beta) * KATS::abs(max_y<scalar_type>()) +
                             10 * KATS::eps() * maxNnzPerRow * KATS::abs(alpha) * KATS::abs(max_a<scalar_type>()) *
                                 KATS::abs(max_x<scalar_type>());

  std::vector<std::pair<ordinal_type, ordinal_type>> errIdx;

  for (ordinal_type i = 0; i < ordinal_type(hyAct.extent(0)); ++i) {
    for (ordinal_type j = 0; j < ordinal_type(hyAct.extent(1)); ++j) {
      if (KATS::abs(hyExp(i, j) - hyAct(i, j)) > tolerance) {
        errIdx.push_back({i, j});
      }
    }
  }

  if (!errIdx.empty()) {
    std::string alg = KokkosSparse::get_spmv_algorithm_name(handle->get_algorithm());

    std::cerr << __FILE__ << ":" << __LINE__ << " BsrMatrix SpMMV failure!" << std::endl;
    std::cerr << "alg:          " << alg << std::endl;
    std::cerr << "mode:         " << mode << std::endl;
    std::cerr << "A:            " << a.numRows() << "x" << a.numCols() << std::endl;
    std::cerr << "A blockdim:   " << a.blockDim() << std::endl;
    std::cerr << "alpha:        " << alpha << std::endl;
    std::cerr << "beta:         " << beta << std::endl;
    std::cerr << "maxNnzPerRow: " << maxNnzPerRow << std::endl;
    std::cerr << "First 100 errors:" << std::endl;
    std::cerr << "i\tj\texp\tact\terr\ttol" << std::endl;
    std::cerr << "-\t-\t---\t---\t---\t---" << std::endl;
    for (size_t e = 0; e < 100 && e < errIdx.size(); ++e) {
      auto ij = errIdx[e];
      auto i  = ij.first;
      auto j  = ij.second;
      // clang-format off
      std::cerr << i << "\t" << j 
                << "\t" << hyExp(i,j)
                << "\t" << hyAct(i,j)
                << "\t" << KATS::abs(hyExp(i,j) - hyAct(i,j))
                << "\t" << tolerance
                << std::endl;
      // clang-format on
    }
  }

  EXPECT_TRUE(errIdx.empty());
}

template <typename Layout, typename Bsr>
struct MultiVectorTypeFor {
  using type = Kokkos::View<typename Bsr::non_const_value_type **, Layout, typename Bsr::device_type>;
};

/*! \brief create random x and y multivectors for a given matrix and spmv mode
 */
template <typename Layout, typename Bsr>
auto random_multivecs_for_spm_mv(const char *mode, const Bsr &a, const size_t numVecs, const bool nans = false)
    -> std::tuple<typename MultiVectorTypeFor<Layout, Bsr>::type, typename MultiVectorTypeFor<Layout, Bsr>::type> {
  using scalar_type     = typename Bsr::non_const_value_type;
  using vector_type     = typename MultiVectorTypeFor<Layout, Bsr>::type;
  using execution_space = typename Bsr::execution_space;
  using policy_type     = Kokkos::RangePolicy<typename vector_type::execution_space>;

  size_t nx = a.numCols() * a.blockDim();
  size_t ny = a.numRows() * a.blockDim();
  if (mode_is_transpose(mode)) {
    std::swap(nx, ny);
  }
  vector_type x("x", nx, numVecs);
  vector_type y("y", ny, numVecs);

  Kokkos::Random_XorShift64_Pool<execution_space> random(13718);
  Kokkos::fill_random(x, random, max_x<scalar_type>());
  Kokkos::fill_random(y, random, max_y<scalar_type>());

  // sprinkle some "random" NaNs in
  if (nans) {
    Kokkos::parallel_for(
        policy_type(0, x.extent(0)), KOKKOS_LAMBDA(size_t i) {
          for (size_t j = 0; j < x.extent(1); ++j) {
            if (0 == ((i * x.extent(1) + j) % 13)) {
              x(i, j) = KokkosKernels::Impl::quiet_NaN<scalar_type>();
            }
          }
        });
    Kokkos::parallel_for(
        policy_type(0, y.extent(0)), KOKKOS_LAMBDA(size_t i) {
          for (size_t j = 0; j < y.extent(1); ++j) {
            if (0 == ((i * y.extent(1) + j) % 17)) {
              y(i, j) = KokkosKernels::Impl::quiet_NaN<scalar_type>();
            }
          }
        });
  }

  return std::make_tuple(x, y);
}

template <typename Layout, typename Bsr, typename Crs>
void test_spm_mv_combos(const char *mode, const Bsr &a, const Crs &acrs, size_t maxNnzPerRow) {
  using namespace KokkosSparse;
  using execution_space = typename Bsr::execution_space;
  using scalar_type     = typename Bsr::non_const_value_type;
  using multivector_t   = typename MultiVectorTypeFor<Layout, Bsr>::type;
  using handle_t        = SPMVHandle<execution_space, Bsr, multivector_t, multivector_t>;

  // cover a variety of algorithms
  std::vector<std::unique_ptr<handle_t>> handles;
  for (SPMVAlgorithm algo : {SPMV_DEFAULT, SPMV_NATIVE, SPMV_BSR_V41})
    handles.push_back(std::make_unique<handle_t>(algo));

  // Tensor core algorithm temporarily disabled, fails on V100
  /*
  if constexpr (KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()) {
#if defined(KOKKOS_ENABLE_CUDA)
    if constexpr (std::is_same_v<execution_space, Kokkos::Cuda>) {
#if defined(KOKKOS_ARCH_AMPERE) || defined(KOKKOS_ARCH_VOLTA)
      handles.push_back(new handle_t(SPMV_BSR_TC));
#if defined(KOKKOS_ARCH_AMPERE)
      // Also call SPMV_BSR_TC with Precision = Double on Ampere
      handles.push_back(new handle_t(SPMV_BSR_TC));
      handles.back()->bsr_tc_precision = Experimental::Bsr_TC_Precision::Double;
#endif  // AMPERE
#endif  // AMPERE || VOLTA
    }
#endif  // CUDA
  }
  */

  for (size_t numVecs : {1, 7}) {  // num multivecs
    auto [x, y]                     = random_multivecs_for_spm_mv<Layout>(mode, a, numVecs);
    auto [x_with_nans, y_with_nans] = random_multivecs_for_spm_mv<Layout>(mode, a, numVecs, true);
    for (std::unique_ptr<handle_t> &handle : handles) {
      for (scalar_type alpha : {scalar_type(0), scalar_type(1), scalar_type(-1), scalar_type(3.7)}) {
        for (scalar_type beta : {scalar_type(0), scalar_type(1), scalar_type(-1), scalar_type(-1.5)}) {
          test_spm_mv(handle.get(), mode, alpha, beta, a, acrs, maxNnzPerRow, x, y);
          if (beta == scalar_type(0)) {
            test_spm_mv(handle.get(), mode, alpha, beta, a, acrs, maxNnzPerRow, x_with_nans, y_with_nans);
          }
        }
      }
    }
  }
}

/*! \brief test all combos of all matrices with different block sizes
 */
template <typename Scalar, typename Ordinal, typename Offset, typename Layout, typename Device>
void test_spm_mv_corner_cases() {
  using Bsr = KokkosSparse::Experimental::BsrMatrix<Scalar, Ordinal, Device, void, Offset>;
  using Crs = KokkosSparse::CrsMatrix<Scalar, Ordinal, Device, void, Offset>;
  for (auto mode : {"N", "T", "C", "H"}) {
    for (int bs : {1, 2, 5, 9}) {
      {
        auto A    = bsr_corner_case_0_by_0<Bsr>(bs);
        auto Acrs = KokkosSparse::Impl::bsr_to_crs<Crs>(A);
        test_spm_mv_combos<Layout>(mode, A, Acrs, 0);
      }
      {
        auto A    = bsr_corner_case_0_by_1<Bsr>(bs);
        auto Acrs = KokkosSparse::Impl::bsr_to_crs<Crs>(A);
        test_spm_mv_combos<Layout>(mode, A, Acrs, 0);
      }
      {
        auto A    = bsr_corner_case_1_by_0<Bsr>(bs);
        auto Acrs = KokkosSparse::Impl::bsr_to_crs<Crs>(A);
        test_spm_mv_combos<Layout>(mode, A, Acrs, 0);
      }
    }
  }
}

template <typename Scalar, typename Ordinal, typename Offset, typename Layout, typename Device>
void test_spm_mv_random() {
  using Bsr = KokkosSparse::Experimental::BsrMatrix<Scalar, Ordinal, Device, void, Offset>;
  using Crs = KokkosSparse::CrsMatrix<Scalar, Ordinal, Device, void, Offset>;
  // thoroughly test smaller matrices
  std::vector<std::pair<int, int>> shapes = {{10, 10}, {10, 50}, {50, 10}};
  for (auto &shape : shapes) {
    for (int bs : {1, 2, 5, 9}) {
      auto A                   = bsr_random<Bsr>(bs, shape.first, shape.second);
      auto Acrs                = KokkosSparse::Impl::bsr_to_crs<Crs>(A);
      size_t maxNnzPerRow      = opMaxNnzPerRow(A, false);
      size_t maxNnzPerRowTrans = opMaxNnzPerRow(A, true);
      for (auto mode : {"N", "T", "C", "H"}) {
        test_spm_mv_combos<Layout>(mode, A, Acrs, mode_is_transpose(mode) ? maxNnzPerRowTrans : maxNnzPerRow);
      }
    }
  }

  // test a tougher case on a big matrix
  {
    constexpr int blockSizePrime = 7;
    constexpr int smallPrime     = 11;
    constexpr int largePrime     = 499;
    auto A                       = bsr_random<Bsr>(blockSizePrime, smallPrime, largePrime);
    auto Acrs                    = KokkosSparse::Impl::bsr_to_crs<Crs>(A);
    size_t maxNnzPerRow          = opMaxNnzPerRow(A, false);
    size_t maxNnzPerRowTrans     = opMaxNnzPerRow(A, true);
    for (auto mode : {"N", "T"}) {
      test_spm_mv_combos<Layout>(mode, A, Acrs, mode_is_transpose(mode) ? maxNnzPerRowTrans : maxNnzPerRow);
    }
  }
}

template <typename Scalar, typename Ordinal, typename Offset, typename Layout, typename Device>
void test_spm_mv() {
  test_spm_mv_corner_cases<Scalar, Ordinal, Offset, Layout, Device>();
  test_spm_mv_random<Scalar, Ordinal, Offset, Layout, Device>();
}

}  // namespace Test_Spmv_Bsr

//////////////////////////

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                        \
  TEST_F(TestCategory, sparse##_##bsr_spmv##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    Test_Spmv_Bsr::test_spmv<SCALAR, ORDINAL, OFFSET, DEVICE>();                           \
  }

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST

//////////////////////////

#define EXECUTE_BSR_TIMES_MVEC_TEST(SCALAR, ORDINAL, OFFSET, LAYOUT, DEVICE)                           \
  TEST_F(TestCategory, sparse##_##bsr_spmmv##_##SCALAR##_##ORDINAL##_##OFFSET##_##LAYOUT##_##DEVICE) { \
    Test_Spmv_Bsr::test_spm_mv<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>();                     \
  }

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
  EXECUTE_BSR_TIMES_MVEC_TEST(SCALAR, ORDINAL, OFFSET, LayoutLeft, TestDevice)

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST

#endif  // KOKKOSKERNELS_INST_LAYOUTLEFT

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
  EXECUTE_BSR_TIMES_MVEC_TEST(SCALAR, ORDINAL, OFFSET, LayoutRight, TestDevice)

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST

#endif  // KOKKOSKERNELS_INST_LAYOUTRIGHT

#undef EXECUTE_BSR_TIMES_MVEC_TEST
