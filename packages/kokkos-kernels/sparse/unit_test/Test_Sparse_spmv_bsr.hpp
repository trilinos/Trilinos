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

  Test the following 768 combos for at least a few matcies.

  Algorithms              Alpha     Beta     Block Sizes    Modes
  (none)                  0         0        1              N
  native              x   1      x  1     x  2           x  T
  experimental_bsr_tc     -1        -1       5              C
                          3.7       -1.5     9              H

  There are also a subset of tests on larger matrices
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
#include "KokkosKernels_Controls.hpp"
#include "KokkosKernels_default_types.hpp"

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
inline bool mode_is_transpose(const char *mode) {
  return mode[0] == 'T' || mode[0] == 'H';
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
  using Crs =
      KokkosSparse::CrsMatrix<scalar_type, ordinal_type,
                              typename Bsr::device_type, void, size_type>;
  using Graph = typename Crs::staticcrsgraph_type;

  // construct a random Crs Matrix
  Test::RandCsMatrix<scalar_type, Kokkos::LayoutLeft, typename Bsr::device_type,
                     ordinal_type, size_type>
      rcs(blockRows, blockCols, scalar_type(0), max_a<scalar_type>());

  const auto colids = Kokkos::subview(
      rcs.get_ids(), Kokkos::make_pair(size_t(0), rcs.get_nnz()));
  const auto vals = Kokkos::subview(
      rcs.get_vals(), Kokkos::make_pair(size_t(0), rcs.get_nnz()));
  Graph graph(colids, rcs.get_map());
  Crs crs("crs", blockCols, vals, graph);

  // expand to Bsr matrix
  return KokkosSparse::Impl::expand_crs_to_bsr<Bsr>(crs, blockSize);
}

/*! \brief reference SpMV is the KokkosSparse::spmv on the equivalent point
 * matrix
 */
template <typename Alpha, typename Bsr, typename XVector, typename Beta,
          typename YVector>
void reference_spmv(const char *mode, const Alpha &alpha, const Bsr &a,
                    const XVector &x, const Beta &beta, const YVector &y) {
  using Crs = KokkosSparse::CrsMatrix<
      typename Bsr::non_const_value_type, typename Bsr::non_const_ordinal_type,
      typename Bsr::device_type, void, typename Bsr::non_const_size_type>;
  const Crs crs = KokkosSparse::Impl::bsr_to_crs<Crs>(a);

  KokkosSparse::spmv(mode, alpha, crs, x, beta, y);
}

/*! \brief test a specific spmv

*/
template <typename Bsr, typename XVector, typename YVector,
          typename Alpha = typename Bsr::non_const_value_type,
          typename Beta  = typename Bsr::non_const_value_type>
void test_spmv(const char *alg, const char *mode, const Alpha &alpha,
               const Beta &beta, const Bsr &a, const XVector &x,
               const YVector &y) {
  using execution_space = typename Bsr::execution_space;
  using scalar_type     = typename Bsr::non_const_value_type;
  using ordinal_type    = typename Bsr::non_const_ordinal_type;
  using KATS            = Kokkos::ArithTraits<scalar_type>;
  using mag_type        = typename KATS::mag_type;

  // generate expected result from reference implementation
  YVector yExp("yExp", y.extent(0));
  Kokkos::deep_copy(yExp, y);
  reference_spmv(mode, alpha, a, x, beta, yExp);

  // scratch space for actual value (don't modify input)
  YVector yAct("yAct", y.extent(0));
  Kokkos::deep_copy(yAct, y);

  if (alg) {
    KokkosKernels::Experimental::Controls controls;
    controls.setParameter("algorithm", alg);
    KokkosSparse::spmv(controls, mode, alpha, a, x, beta, yAct);
  } else {
    KokkosSparse::spmv(mode, alpha, a, x, beta, yAct);
  }

  // compare yExp and yAct
  auto hyExp = Kokkos::create_mirror_view(yExp);
  auto hyAct = Kokkos::create_mirror_view(yAct);
  Kokkos::deep_copy(hyExp, yExp);
  Kokkos::deep_copy(hyAct, yAct);

  // max nnz per row is used for the tolerance
  // for a transposed computation, need to transpose the matrix before
  // seeing which rows are longest
  size_t maxNnzPerRow;
  if (mode_is_transpose(mode)) {
    auto at = KokkosSparse::Impl::transpose_bsr_matrix(a);
    maxNnzPerRow =
        at.blockDim() *
        KokkosSparse::Impl::graph_max_degree<execution_space, ordinal_type>(
            at.graph.row_map);
  } else {
    maxNnzPerRow =
        a.blockDim() *
        KokkosSparse::Impl::graph_max_degree<execution_space, ordinal_type>(
            a.graph.row_map);
  }

  /* assume that any floating-point op may introduce eps() error
     scaling y is one op
     dot product of x is two ops per entry (mul and add)

     10x means same order of magnitude
  */
  const mag_type tolerance =
      KATS::eps() * KATS::abs(beta) * KATS::abs(max_y<scalar_type>()) +
      10 * KATS::eps() * maxNnzPerRow * KATS::abs(alpha) *
          KATS::abs(max_a<scalar_type>()) * KATS::abs(max_x<scalar_type>());

  std::vector<ordinal_type> errIdx;

  for (ordinal_type i = 0; i < ordinal_type(hyAct.extent(0)); ++i) {
    if (KATS::abs(hyExp(i) - hyAct(i)) > tolerance) {
      errIdx.push_back(i);
    }
  }

  if (!errIdx.empty()) {
    std::cerr << __FILE__ << ":" << __LINE__ << " BsrMatrix SpMV failure!"
              << std::endl;
    std::cerr << "alg:          " << (alg ? alg : "<none>") << std::endl;
    std::cerr << "mode:         " << mode << std::endl;
    std::cerr << "A:            " << a.numRows() << "x" << a.numCols()
              << std::endl;
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
  using type = Kokkos::View<typename Bsr::non_const_value_type *,
                            typename Bsr::device_type>;
};

template <typename Bsr>
std::tuple<Bsr, typename VectorTypeFor<Bsr>::type,
           typename VectorTypeFor<Bsr>::type>
spmv_corner_case_0_by_0(const char * /*mode*/, const int blockSize) {
  using vector_type = typename VectorTypeFor<Bsr>::type;
  Bsr a             = bsr_corner_case_0_by_0<Bsr>(blockSize);
  vector_type x("x", 0);
  vector_type y("y", 0);
  return std::make_tuple(a, x, y);
}

template <typename Bsr>
std::tuple<Bsr, typename VectorTypeFor<Bsr>::type,
           typename VectorTypeFor<Bsr>::type>
spmv_corner_case_0_by_1(const char *mode, const int blockSize) {
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
std::tuple<Bsr, typename VectorTypeFor<Bsr>::type,
           typename VectorTypeFor<Bsr>::type>
spmv_corner_case_1_by_0(const char *mode, const int blockSize) {
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
std::tuple<Bsr, typename VectorTypeFor<Bsr>::type,
           typename VectorTypeFor<Bsr>::type>
spmv_random(const char *mode, const int blockSize, const int blockRows,
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
auto random_vecs_for_spmv(const char *mode, const Bsr &a) {
  using scalar_type     = typename Bsr::non_const_value_type;
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

  return std::make_tuple(x, y);
}

/*! \brief test all combos of the provided matrix
 */
template <typename Bsr>
void test_spmv_combos(const char *mode, const Bsr &a) {
  using scalar_type = typename Bsr::non_const_value_type;

  auto [x, y] = random_vecs_for_spmv(mode, a);

  for (auto alg : {(const char *)(nullptr), "native", "experimental_tc_bsr"}) {
    for (scalar_type alpha :
         {scalar_type(0), scalar_type(1), scalar_type(-1), scalar_type(3.7)}) {
      for (scalar_type beta : {scalar_type(0), scalar_type(1), scalar_type(-1),
                               scalar_type(-1.5)}) {
        test_spmv(alg, mode, alpha, beta, a, x, y);
      }
    }
  }
}

/*! \brief test all combos of all matrices with different block sizes
 */
template <typename Scalar, typename Ordinal, typename Offset, typename Device>
void test_spmv_corner_cases() {
  using Bsr = KokkosSparse::Experimental::BsrMatrix<Scalar, Ordinal, Device,
                                                    void, Offset>;
  for (auto mode : {"N", "T", "C", "H"}) {
    for (int bs : {1, 2, 5, 9}) {
      test_spmv_combos(mode, bsr_corner_case_0_by_0<Bsr>(bs));
      test_spmv_combos(mode, bsr_corner_case_0_by_1<Bsr>(bs));
      test_spmv_combos(mode, bsr_corner_case_1_by_0<Bsr>(bs));
    }
  }
}

template <typename Scalar, typename Ordinal, typename Offset, typename Device>
void test_spmv_random() {
  using Bsr = KokkosSparse::Experimental::BsrMatrix<Scalar, Ordinal, Device,
                                                    void, Offset>;
  for (auto mode : {"N", "T", "C", "H"}) {
    for (int bs : {1, 2, 5, 9}) {
      test_spmv_combos(mode, bsr_random<Bsr>(bs, 10, 10));
      test_spmv_combos(mode, bsr_random<Bsr>(bs, 10, 50));
      test_spmv_combos(mode, bsr_random<Bsr>(bs, 50, 10));
    }
  }

  // test a tougher case on a big matrix
  constexpr int blockSizePrime = 7;
  constexpr int smallPrime     = 11;
  constexpr int largePrime     = 499;
  for (auto mode : {"N", "T"}) {
    test_spmv_combos(mode,
                     bsr_random<Bsr>(blockSizePrime, smallPrime, largePrime));
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

template <typename Bsr, typename XVector, typename YVector, typename Alpha,
          typename Beta>
void test_spm_mv(const char *alg, const char *mode, const Alpha &alpha,
                 const Beta &beta, const Bsr &a, const XVector &x,
                 const YVector &y) {
  using execution_space = typename Bsr::execution_space;
  using scalar_type     = typename Bsr::non_const_value_type;
  using ordinal_type    = typename Bsr::non_const_ordinal_type;
  using KATS            = Kokkos::ArithTraits<scalar_type>;
  using mag_type        = typename KATS::mag_type;

  // generate expected result from reference implementation
  YVector yExp("yExp", y.extent(0), y.extent(1));
  Kokkos::deep_copy(yExp, y);
  reference_spmv(mode, alpha, a, x, beta, yExp);

  // scratch space for actual value (don't modify input)
  YVector yAct("yAct", y.extent(0), y.extent(1));
  Kokkos::deep_copy(yAct, y);

  if (alg) {
    KokkosKernels::Experimental::Controls controls;
    controls.setParameter("algorithm", alg);
    KokkosSparse::spmv(controls, mode, alpha, a, x, beta, yAct);
  } else {
    KokkosSparse::spmv(mode, alpha, a, x, beta, yAct);
  }

  // compare yExp and yAct
  auto hyExp = Kokkos::create_mirror_view(yExp);
  auto hyAct = Kokkos::create_mirror_view(yAct);
  Kokkos::deep_copy(hyExp, yExp);
  Kokkos::deep_copy(hyAct, yAct);

  // max nnz per row is used for the tolerance
  // for a transposed computation, need to transpose the matrix before
  // seeing which rows are longest
  size_t maxNnzPerRow;
  if (mode_is_transpose(mode)) {
    auto at = KokkosSparse::Impl::transpose_bsr_matrix(a);
    maxNnzPerRow =
        at.blockDim() *
        KokkosSparse::Impl::graph_max_degree<execution_space, ordinal_type>(
            at.graph.row_map);
  } else {
    maxNnzPerRow =
        a.blockDim() *
        KokkosSparse::Impl::graph_max_degree<execution_space, ordinal_type>(
            a.graph.row_map);
  }

  /* assume that any floating-point op may introduce eps() error
     scaling y is one op
     dot product of x is two ops per entry (mul and add)
  */
  const mag_type tolerance =
      KATS::eps() * KATS::abs(beta) * KATS::abs(max_y<scalar_type>()) +
      10 * KATS::eps() * maxNnzPerRow * KATS::abs(alpha) *
          KATS::abs(max_a<scalar_type>()) * KATS::abs(max_x<scalar_type>());

  std::vector<std::pair<ordinal_type, ordinal_type>> errIdx;

  for (ordinal_type i = 0; i < ordinal_type(hyAct.extent(0)); ++i) {
    for (ordinal_type j = 0; j < ordinal_type(hyAct.extent(1)); ++j) {
      if (KATS::abs(hyExp(i, j) - hyAct(i, j)) > tolerance) {
        errIdx.push_back({i, j});
      }
    }
  }

  if (!errIdx.empty()) {
    std::cerr << __FILE__ << ":" << __LINE__ << " BsrMatrix SpMMV failure!"
              << std::endl;
    std::cerr << "alg:          " << (alg ? alg : "<none>") << std::endl;
    std::cerr << "mode:         " << mode << std::endl;
    std::cerr << "A:            " << a.numRows() << "x" << a.numCols()
              << std::endl;
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
  using type = Kokkos::View<typename Bsr::non_const_value_type **, Layout,
                            typename Bsr::device_type>;
};

/*! \brief create random x and y multivectors for a given matrix and spmv mode
 */
template <typename Layout, typename Bsr>
auto random_multivecs_for_spm_mv(const char *mode, const Bsr &a,
                                 const size_t numVecs) {
  using scalar_type     = typename Bsr::non_const_value_type;
  using vector_type     = typename MultiVectorTypeFor<Layout, Bsr>::type;
  using execution_space = typename Bsr::execution_space;

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

  return std::make_tuple(x, y);
}

template <typename Layout, typename Bsr>
void test_spm_mv_combos(const char *mode, const Bsr &a) {
  using scalar_type = typename Bsr::non_const_value_type;

  for (size_t numVecs : {1, 2, 7}) {  // num multivecs
    auto [x, y] = random_multivecs_for_spm_mv<Layout>(mode, a, numVecs);
    for (auto alg :
         {(const char *)(nullptr), "native", "experimental_tc_bsr"}) {
      for (scalar_type alpha : {scalar_type(0), scalar_type(1), scalar_type(-1),
                                scalar_type(3.7)}) {
        for (scalar_type beta : {scalar_type(0), scalar_type(1),
                                 scalar_type(-1), scalar_type(-1.5)}) {
          test_spm_mv(alg, mode, alpha, beta, a, x, y);
        }
      }
    }
  }
}

/*! \brief test all combos of all matrices with different block sizes
 */
template <typename Scalar, typename Ordinal, typename Offset, typename Layout,
          typename Device>
void test_spm_mv_corner_cases() {
  using Bsr = KokkosSparse::Experimental::BsrMatrix<Scalar, Ordinal, Device,
                                                    void, Offset>;
  for (auto mode : {"N", "T", "C", "H"}) {
    for (int bs : {1, 2, 5, 9}) {
      test_spm_mv_combos<Layout>(mode, bsr_corner_case_0_by_0<Bsr>(bs));
      test_spm_mv_combos<Layout>(mode, bsr_corner_case_0_by_1<Bsr>(bs));
      test_spm_mv_combos<Layout>(mode, bsr_corner_case_1_by_0<Bsr>(bs));
    }
  }
}

template <typename Scalar, typename Ordinal, typename Offset, typename Layout,
          typename Device>
void test_spm_mv_random() {
  using Bsr = KokkosSparse::Experimental::BsrMatrix<Scalar, Ordinal, Device,
                                                    void, Offset>;
  // thoroughly test smaller matrices
  for (auto mode : {"N", "T", "C", "H"}) {
    for (int bs : {1, 2, 5, 9}) {
      test_spm_mv_combos<Layout>(mode, bsr_random<Bsr>(bs, 10, 10));
      test_spm_mv_combos<Layout>(mode, bsr_random<Bsr>(bs, 10, 50));
      test_spm_mv_combos<Layout>(mode, bsr_random<Bsr>(bs, 50, 10));
    }
  }

  // test a tougher case on a big matrix
  constexpr int blockSizePrime = 7;
  constexpr int smallPrime     = 11;
  constexpr int largePrime     = 499;
  for (auto mode : {"N", "T"}) {
    test_spm_mv_combos<Layout>(
        mode, bsr_random<Bsr>(blockSizePrime, smallPrime, largePrime));
  }
}

template <typename Scalar, typename Ordinal, typename Offset, typename Layout,
          typename Device>
void test_spm_mv() {
  test_spm_mv_corner_cases<Scalar, Ordinal, Offset, Layout, Device>();
  test_spm_mv_random<Scalar, Ordinal, Offset, Layout, Device>();
}

}  // namespace Test_Spmv_Bsr

//////////////////////////

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)          \
  TEST_F(TestCategory,                                                       \
         sparse##_##bsr_spmv##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    Test_Spmv_Bsr::test_spmv<SCALAR, ORDINAL, OFFSET, DEVICE>();             \
  }

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST

//////////////////////////

#define EXECUTE_BSR_TIMES_MVEC_TEST(SCALAR, ORDINAL, OFFSET, LAYOUT, DEVICE)          \
  TEST_F(                                                                             \
      TestCategory,                                                                   \
      sparse##_##bsr_spmmv##_##SCALAR##_##ORDINAL##_##OFFSET##_##LAYOUT##_##DEVICE) { \
    Test_Spmv_Bsr::test_spm_mv<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT,               \
                               DEVICE>();                                             \
  }

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
  EXECUTE_BSR_TIMES_MVEC_TEST(SCALAR, ORDINAL, OFFSET, LayoutLeft,  \
                              TestExecSpace)

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST

#endif  // KOKKOSKERNELS_INST_LAYOUTLEFT

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
  EXECUTE_BSR_TIMES_MVEC_TEST(SCALAR, ORDINAL, OFFSET, LayoutRight, \
                              TestExecSpace)

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST

#endif  // KOKKOSKERNELS_INST_LAYOUTRIGHT

#undef EXECUTE_BSR_TIMES_MVEC_TEST
