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

#ifndef KOKKOSKERNELS_TEST_UTILS_HPP
#define KOKKOSKERNELS_TEST_UTILS_HPP

#include <random>

#include "KokkosKernels_Utils.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "KokkosBatched_Vector.hpp"
// Make this include-able from all subdirectories
#include "../tpls/gtest/gtest/gtest.h"  //for EXPECT_**

// Simplify ETI macros
#if !defined(KOKKOSKERNELS_ETI_ONLY) && \
    !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
#define KOKKOSKERNELS_TEST_ALL_TYPES
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    defined(KOKKOSKERNELS_TEST_ALL_TYPES)
#define KOKKOSKERNELS_TEST_LAYOUTLEFT
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    defined(KOKKOSKERNELS_TEST_ALL_TYPES)
#define KOKKOSKERNELS_TEST_LAYOUTRIGHT
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || \
    defined(KOKKOSKERNELS_TEST_ALL_TYPES)
#define KOKKOSKERNELS_TEST_LAYOUTSTRIDE
#endif
#if defined(KOKKOSKERNELS_INST_FLOAT) || defined(KOKKOSKERNELS_TEST_ALL_TYPES)
#define KOKKOSKERNELS_TEST_FLOAT
#endif
#if defined(KOKKOSKERNELS_INST_DOUBLE) || defined(KOKKOSKERNELS_TEST_ALL_TYPES)
#define KOKKOSKERNELS_TEST_DOUBLE
#endif
#if defined(KOKKOSKERNELS_INST_INT) || defined(KOKKOSKERNELS_TEST_ALL_TYPES)
#define KOKKOSKERNELS_TEST_INT
#endif
#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || \
    defined(KOKKOSKERNELS_TEST_ALL_TYPES)
#define KOKKOSKERNELS_TEST_COMPLEX_FLOAT
#endif
#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    defined(KOKKOSKERNELS_TEST_ALL_TYPES)
#define KOKKOSKERNELS_TEST_COMPLEX_DOUBLE
#endif

namespace Test {
template <class ViewType,
          bool strided = std::is_same<typename ViewType::array_layout,
                                      Kokkos::LayoutStride>::value>
struct multivector_layout_adapter;

template <class ViewType>
struct multivector_layout_adapter<ViewType, true> {
  typedef typename ViewType::value_type Scalar;
  typedef typename ViewType::device_type Device;
  typedef Kokkos::View<Scalar* * [2], Kokkos::LayoutRight, Device>
      BaseTypeRight;
  typedef Kokkos::View<Scalar**, typename ViewType::array_layout, Device>
      BaseTypeDefault;
  typedef
      typename std::conditional<std::is_same<typename ViewType::array_layout,
                                             Kokkos::LayoutStride>::value,
                                BaseTypeRight, BaseTypeDefault>::type BaseType;

  static ViewType view(const BaseType& v) {
    return Kokkos::subview(v, Kokkos::ALL, Kokkos::ALL, 0);
  };
};

template <class ViewType>
struct multivector_layout_adapter<ViewType, false> {
  typedef typename ViewType::value_type Scalar;
  typedef typename ViewType::device_type Device;
  typedef Kokkos::View<Scalar* * [2], Kokkos::LayoutRight, Device>
      BaseTypeRight;
  typedef Kokkos::View<Scalar**, typename ViewType::array_layout, Device>
      BaseTypeDefault;
  typedef
      typename std::conditional<std::is_same<typename ViewType::array_layout,
                                             Kokkos::LayoutStride>::value,
                                BaseTypeRight, BaseTypeDefault>::type BaseType;

  static ViewType view(const BaseType& v) {
    return Kokkos::subview(v, Kokkos::ALL, Kokkos::ALL);
  };
};

template <class Scalar1, class Scalar2, class Scalar3>
void EXPECT_NEAR_KK(Scalar1 val1, Scalar2 val2, Scalar3 tol,
                    std::string msg = "") {
  typedef Kokkos::Details::ArithTraits<Scalar1> AT1;
  typedef Kokkos::Details::ArithTraits<Scalar3> AT3;
  EXPECT_LE((double)AT1::abs(val1 - val2), (double)AT3::abs(tol)) << msg;
}

template <class Scalar1, class Scalar2, class Scalar3>
void EXPECT_NEAR_KK_REL(Scalar1 val1, Scalar2 val2, Scalar3 tol,
                        std::string msg = "") {
  typedef typename std::remove_reference<decltype(val1)>::type hv1_type;
  typedef typename std::remove_reference<decltype(val2)>::type hv2_type;
  const auto ahv1 = Kokkos::Details::ArithTraits<hv1_type>::abs(val1);
  const auto ahv2 = Kokkos::Details::ArithTraits<hv2_type>::abs(val2);
  EXPECT_NEAR_KK(val1, val2, tol * Kokkos::max(ahv1, ahv2), msg);
}

// Special overload for accurate value by value SIMD vectors comparison
template <class Scalar, int VecLen, class Tolerance>
void EXPECT_NEAR_KK_REL(
    const KokkosBatched::Vector<KokkosBatched::SIMD<Scalar>, VecLen>& val1,
    const KokkosBatched::Vector<KokkosBatched::SIMD<Scalar>, VecLen>& val2,
    Tolerance tol, std::string msg = "") {
  for (int i = 0; i < VecLen; ++i) {
    EXPECT_NEAR_KK_REL(val1[i], val2[i], tol, msg);
  }
}

template <class ViewType1, class ViewType2, class Scalar>
void EXPECT_NEAR_KK_1DVIEW(ViewType1 v1, ViewType2 v2, Scalar tol) {
  size_t v1_size = v1.extent(0);
  size_t v2_size = v2.extent(0);
  EXPECT_EQ(v1_size, v2_size);

  typename ViewType1::HostMirror h_v1 = Kokkos::create_mirror_view(v1);
  typename ViewType2::HostMirror h_v2 = Kokkos::create_mirror_view(v2);

  KokkosKernels::Impl::safe_device_to_host_deep_copy(v1.extent(0), v1, h_v1);
  KokkosKernels::Impl::safe_device_to_host_deep_copy(v2.extent(0), v2, h_v2);

  for (size_t i = 0; i < v1_size; ++i) {
    EXPECT_NEAR_KK(h_v1(i), h_v2(i), tol);
  }
}

template <class ViewType1, class ViewType2, class Scalar>
void EXPECT_NEAR_KK_REL_1DVIEW(ViewType1 v1, ViewType2 v2, Scalar tol) {
  size_t v1_size = v1.extent(0);
  size_t v2_size = v2.extent(0);
  EXPECT_EQ(v1_size, v2_size);

  typename ViewType1::HostMirror h_v1 = Kokkos::create_mirror_view(v1);
  typename ViewType2::HostMirror h_v2 = Kokkos::create_mirror_view(v2);

  KokkosKernels::Impl::safe_device_to_host_deep_copy(v1.extent(0), v1, h_v1);
  KokkosKernels::Impl::safe_device_to_host_deep_copy(v2.extent(0), v2, h_v2);

  for (size_t i = 0; i < v1_size; ++i) {
    EXPECT_NEAR_KK_REL(h_v1(i), h_v2(i), tol);
  }
}

/// This function returns a descriptive user defined failure string for
/// insertion into gtest macros such as FAIL() and EXPECT_LE(). \param file The
/// filename where the failure originated \param func The function where the
/// failure originated \param line The line number where the failure originated
/// \return a new string containing: "  > from file:func:line\n    > "
static inline const std::string kk_failure_str(std::string file,
                                               std::string func,
                                               const int line) {
  std::string failure_msg = "  > from ";
  failure_msg += (file + ":" + func + ":" + std::to_string(line) + "\n    > ");
  return std::string(failure_msg);
}

#if defined(KOKKOS_HALF_T_IS_FLOAT)
using halfScalarType = Kokkos::Experimental::half_t;
#endif  // KOKKOS_HALF_T_IS_FLOAT

#if defined(KOKKOS_BHALF_T_IS_FLOAT)
using bhalfScalarType = Kokkos::Experimental::bhalf_t;
#endif  // KOKKOS_BHALF_T_IS_FLOAT

template <class ViewTypeA, class ViewTypeB, class ViewTypeC,
          class ExecutionSpace>
struct SharedVanillaGEMM {
  bool A_t, B_t, A_c, B_c;
  int C_rows, C_cols, A_cols;
  ViewTypeA A;
  ViewTypeB B;
  ViewTypeC C;

  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;
  typedef typename ViewTypeC::value_type ScalarC;
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride,
                       typename ViewTypeA::device_type>
      SubviewTypeA;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutStride,
                       typename ViewTypeB::device_type>
      SubviewTypeB;
  typedef Kokkos::Details::ArithTraits<ScalarC> APT;
  typedef typename APT::mag_type mag_type;
  ScalarA alpha;
  ScalarC beta;

  KOKKOS_INLINE_FUNCTION
  void operator()(
      const typename Kokkos::TeamPolicy<ExecutionSpace>::member_type& team)
      const {
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, C_rows), [&](const int& i) {
          // Give each kokkos thread a vector of A
          SubviewTypeA a_vec;
          if (A_t)
            a_vec = Kokkos::subview(A, Kokkos::ALL(), i);
          else
            a_vec = Kokkos::subview(A, i, Kokkos::ALL());

          // Have all vector lanes perform the dot product
          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(team, C_cols), [&](const int& j) {
                SubviewTypeB b_vec;
                if (B_t)
                  b_vec = Kokkos::subview(B, j, Kokkos::ALL());
                else
                  b_vec = Kokkos::subview(B, Kokkos::ALL(), j);
                ScalarC ab = ScalarC(0);
                for (int k = 0; k < A_cols; k++) {
                  auto a = A_c ? APT::conj(a_vec(k)) : a_vec(k);
                  auto b = B_c ? APT::conj(b_vec(k)) : b_vec(k);
                  ab += a * b;
                }
                C(i, j) = beta * C(i, j) + alpha * ab;
              });
        });
  }
};
// C(i,:,:) = alpha * (A(i,:,:) * B(i,:,:)) + beta * C(i,:,:)
template <class ViewTypeA, class ViewTypeB, class ViewTypeC,
          class ExecutionSpace>
struct Functor_BatchedVanillaGEMM {
  bool A_t, B_t, A_c, B_c, batch_size_last_dim = false;
  ViewTypeA A;
  ViewTypeB B;
  ViewTypeC C;

  using ScalarA      = typename ViewTypeA::value_type;
  using ScalarB      = typename ViewTypeB::value_type;
  using ScalarC      = typename ViewTypeC::value_type;
  using SubviewTypeA = typename Kokkos::View<ScalarA**, Kokkos::LayoutStride,
                                             typename ViewTypeA::device_type>;
  using SubviewTypeB = typename Kokkos::View<ScalarB**, Kokkos::LayoutStride,
                                             typename ViewTypeA::device_type>;
  using SubviewTypeC = typename Kokkos::View<ScalarC**, Kokkos::LayoutStride,
                                             typename ViewTypeA::device_type>;

  ScalarA alpha;
  ScalarC beta;

  KOKKOS_INLINE_FUNCTION
  void operator()(
      const typename Kokkos::TeamPolicy<ExecutionSpace>::member_type& team)
      const {
    int i = team.league_rank();
    SubviewTypeA _A;
    SubviewTypeB _B;
    SubviewTypeC _C;

    if (batch_size_last_dim) {
      _A = Kokkos::subview(A, Kokkos::ALL(), Kokkos::ALL(), i);
      _B = Kokkos::subview(B, Kokkos::ALL(), Kokkos::ALL(), i);
      _C = Kokkos::subview(C, Kokkos::ALL(), Kokkos::ALL(), i);
    } else {
      _A = Kokkos::subview(A, i, Kokkos::ALL(), Kokkos::ALL());
      _B = Kokkos::subview(B, i, Kokkos::ALL(), Kokkos::ALL());
      _C = Kokkos::subview(C, i, Kokkos::ALL(), Kokkos::ALL());
    }
    struct SharedVanillaGEMM<SubviewTypeA, SubviewTypeB, SubviewTypeC,
                             ExecutionSpace>
        vgemm;
    vgemm.A_t    = A_t;
    vgemm.B_t    = B_t;
    vgemm.A_c    = A_c;
    vgemm.B_c    = B_c;
    vgemm.C_rows = batch_size_last_dim ? C.extent(0) : C.extent(1);
    vgemm.C_cols = batch_size_last_dim ? C.extent(1) : C.extent(2);
    vgemm.A_cols = batch_size_last_dim ? (A_t ? A.extent(0) : A.extent(1))
                                       : (A_t ? A.extent(1) : A.extent(2));
    vgemm.A     = _A;
    vgemm.B     = _B;
    vgemm.C     = _C;
    vgemm.alpha = alpha;
    vgemm.beta  = beta;
    vgemm(team);
  }

  inline void run() {
    Kokkos::parallel_for(
        "Test::VanillaGEMM",
        Kokkos::TeamPolicy<ExecutionSpace>(
            batch_size_last_dim ? C.extent(2) : C.extent(0), Kokkos::AUTO,
            KokkosKernels::Impl::kk_get_max_vector_size<ExecutionSpace>()),
        *this);
  }
};

// Compute C := alpha * AB + beta * C
template <class ViewTypeA, class ViewTypeB, class ViewTypeC>
void vanillaGEMM(typename ViewTypeC::non_const_value_type alpha,
                 const ViewTypeA& A, const ViewTypeB& B,
                 typename ViewTypeC::non_const_value_type beta,
                 const ViewTypeC& C) {
  using value_type = typename ViewTypeC::non_const_value_type;
  using KAT        = Kokkos::ArithTraits<value_type>;
  int m            = A.extent(0);
  int k            = A.extent(1);
  int n            = B.extent(1);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      value_type sum = KAT::zero();
      for (int ii = 0; ii < k; ii++) {
        sum += A(i, ii) * B(ii, j);
      }
      C(i, j) = alpha * sum + beta * C(i, j);
    }
  }
}

template <class AlphaType, class ViewTypeA, class ViewTypeX, class BetaType,
          class ViewTypeY>
KOKKOS_INLINE_FUNCTION void vanillaGEMV(char mode, AlphaType alpha,
                                        const ViewTypeA& A, const ViewTypeX& x,
                                        BetaType beta, const ViewTypeY& y) {
  using ScalarY = typename ViewTypeY::non_const_value_type;
  using KAT_A   = Kokkos::ArithTraits<typename ViewTypeA::non_const_value_type>;
  const bool transposed = mode == 'T' || mode == 'C';
  const bool conjugated = mode == 'C';
  const bool has_beta   = beta != Kokkos::ArithTraits<BetaType>::zero();
  int M                 = A.extent(transposed ? 1 : 0);
  int N                 = A.extent(transposed ? 0 : 1);
  for (int i = 0; i < M; i++) {
    ScalarY y_i{};
    if (has_beta) y_i = beta * y(i);
    for (int j = 0; j < N; j++) {
      const auto a   = transposed ? A(j, i) : A(i, j);
      const auto Aij = conjugated ? KAT_A::conj(a) : a;
      y_i += alpha * Aij * x(j);
    }
    y(i) = y_i;
  }
}

template <class T>
class epsilon {
 public:
  constexpr static double value = std::numeric_limits<T>::epsilon();
};

// explicit epsilon specializations
#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
template <>
class epsilon<Kokkos::Experimental::half_t> {
 public:
  constexpr static double value = 0.0009765625F;
};
#endif  // KOKKOS_HALF_T_IS_FLOAT

// explicit epsilon specializations
#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
template <>
class epsilon<Kokkos::Experimental::bhalf_t> {
 public:
  constexpr static double value = 0.0078125F;
};
#endif  // KOKKOS_HALF_T_IS_FLOAT

using KokkosKernels::Impl::getRandomBounds;

template <typename scalar_t, typename lno_t, typename size_type,
          typename device, typename crsMat_t>
crsMat_t symmetrize(crsMat_t A) {
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  auto host_rowmap =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
  auto host_entries =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);
  auto host_values =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.values);
  lno_t numRows = A.numRows();
  // symmetrize as input_mat + input_mat^T, to still have a diagonally dominant
  // matrix
  typedef std::map<lno_t, scalar_t> Row;
  std::vector<Row> symRows(numRows);
  for (lno_t r = 0; r < numRows; r++) {
    auto& row = symRows[r];
    for (size_type i = host_rowmap(r); i < host_rowmap(r + 1); i++) {
      lno_t c   = host_entries(i);
      auto& col = symRows[c];
      auto it   = row.find(c);
      if (it == row.end())
        row[c] = host_values(i);
      else
        row[c] += host_values(i);
      it = col.find(r);
      if (it == col.end())
        col[r] = host_values(i);
      else
        col[r] += host_values(i);
    }
  }
  // Count entries
  Kokkos::View<size_type*, Kokkos::LayoutLeft, Kokkos::HostSpace>
      new_host_rowmap("Rowmap", numRows + 1);
  size_t accum = 0;
  for (lno_t r = 0; r <= numRows; r++) {
    new_host_rowmap(r) = accum;
    if (r < numRows) accum += symRows[r].size();
  }
  // Allocate new entries/values
  Kokkos::View<lno_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> new_host_entries(
      "Entries", accum);
  Kokkos::View<scalar_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>
      new_host_values("Values", accum);
  for (lno_t r = 0; r < numRows; r++) {
    auto rowIt = symRows[r].begin();
    for (size_type i = new_host_rowmap(r); i < new_host_rowmap(r + 1); i++) {
      new_host_entries(i) = rowIt->first;
      new_host_values(i)  = rowIt->second;
      rowIt++;
    }
  }
  lno_view_t new_rowmap("Rowmap", numRows + 1);
  lno_nnz_view_t new_entries("Entries", accum);
  scalar_view_t new_values("Values", accum);
  Kokkos::deep_copy(new_rowmap, new_host_rowmap);
  Kokkos::deep_copy(new_entries, new_host_entries);
  Kokkos::deep_copy(new_values, new_host_values);
  return crsMat_t("SymA", numRows, numRows, accum, new_values, new_rowmap,
                  new_entries);
}

// create_random_x_vector and create_random_y_vector can be used together to
// generate a random linear system Ax = y.
template <typename vec_t>
vec_t create_random_x_vector(vec_t& kok_x, double max_value = 10.0) {
  typedef typename vec_t::value_type scalar_t;
  auto h_x = Kokkos::create_mirror_view(kok_x);
  for (size_t j = 0; j < h_x.extent(1); ++j) {
    for (size_t i = 0; i < h_x.extent(0); ++i) {
      scalar_t r = static_cast<scalar_t>(rand()) /
                   static_cast<scalar_t>(RAND_MAX / max_value);
      h_x.access(i, j) = r;
    }
  }
  Kokkos::deep_copy(kok_x, h_x);
  return kok_x;
}

/// \brief SharedParamTag class used to specify how to invoke templates within
///                       batched unit tests
/// \var TA Indicates which transpose operation to apply to the A matrix
/// \var TB Indicates which transpose operation to apply to the B matrix
/// \var BL Indicates whether the batch size is in the leftmost or rightmost
///         dimension
template <typename TA, typename TB, typename BL>
struct SharedParamTag {
  using transA      = TA;
  using transB      = TB;
  using batchLayout = BL;
};

/// \brief value_type_name returns a string with the value type name
template <typename T>
std::string value_type_name() {
  return "::UnknownValueType";
}

template <>
std::string value_type_name<float>() {
  return "::Float";
}

template <>
std::string value_type_name<double>() {
  return "::Double";
}

template <>
std::string value_type_name<int>() {
  return "::Int";
}

template <>
std::string value_type_name<Kokkos::complex<float>>() {
  return "::ComplexFloat";
}

template <>
std::string value_type_name<Kokkos::complex<double>>() {
  return "::ComplexDouble";
}

int string_compare_no_case(const char* str1, const char* str2) {
  std::string str1_s(str1);
  std::string str2_s(str2);
  for (size_t i = 0; i < str1_s.size(); i++)
    str1_s[i] = std::tolower(str1_s[i]);
  for (size_t i = 0; i < str2_s.size(); i++)
    str2_s[i] = std::tolower(str2_s[i]);
  return strcmp(str1_s.c_str(), str2_s.c_str());
}

/// /brief Cs (Compressed Sparse) matrix class for testing purposes.
/// This class is for testing purposes only and will generate a random
/// Crs / Ccs matrix when instantiated. The class is intentionally written
/// without the use of "row" and "column" member names.
/// dim1 refers to either rows for Crs matrix or columns for a Ccs matrix.
/// dim2 refers to either columns for a Crs matrix or rows for a Ccs matrix.
/// \tparam ScalarType
/// \tparam LayoutType
/// \tparam ExeSpaceType
template <class ScalarType, class LayoutType, class ExeSpaceType>
class RandCsMatrix {
 private:
  using ValViewTypeD = Kokkos::View<ScalarType*, LayoutType, ExeSpaceType>;
  using IdViewTypeD  = Kokkos::View<int64_t*, LayoutType, ExeSpaceType>;
  using MapViewTypeD = Kokkos::View<int64_t*, LayoutType, ExeSpaceType>;
  int64_t __dim2;
  int64_t __dim1;
  int64_t __nnz = 0;
  MapViewTypeD __map_d;
  IdViewTypeD __ids_d;
  ValViewTypeD __vals_d;
  using MapViewTypeH = typename MapViewTypeD::HostMirror;
  using IdViewTypeH  = typename IdViewTypeD::HostMirror;
  using ValViewTypeH = typename ValViewTypeD::HostMirror;
  MapViewTypeH __map;
  IdViewTypeH __ids;
  ValViewTypeH __vals;
  bool __fully_sparse;

  /// Generates a random map where (using Ccs terminology):
  ///  1. __map(i) is in [__ids.data(), &row_ids.data()[nnz - 1]
  ///  2. __map(i) > col_map(i - 1) for i > 1
  ///  3. __map(i) == col_map(j) iff __map(i) == col_map(j) == nullptr
  ///  4. __map(i) - col_map(i - 1) is in [0, m]
  void __populate_random_cs_mat(uint64_t ticks) {
    std::srand(ticks);
    for (int64_t col_idx = 0; col_idx < __dim1; col_idx++) {
      int64_t r = std::rand() % (__dim2 + 1);
      if (r == 0 || __fully_sparse) {  // 100% sparse vector
        __map(col_idx) = __nnz;
      } else {  // sparse vector with r elements
        // Populate r row ids
        std::vector<int64_t> v(r);

        for (int64_t i = 0; i < r; i++) v.at(i) = i;

        std::shuffle(v.begin(), v.end(), std::mt19937(std::random_device()()));

        for (int64_t i = 0; i < r; i++) __ids(i + __nnz) = v.at(i);

        // Point to new column and accumulate number of non zeros
        __map(col_idx) = __nnz;
        __nnz += r;
      }
    }

    // last entry in map points to end of id list
    __map(__dim1) = __nnz;

    // Copy to device
    Kokkos::deep_copy(__map_d, __map);
    Kokkos::deep_copy(__ids_d, __ids);
    ExeSpaceType().fence();
  }

  template <class T>
  T __getter_copy_helper(T src) {
    T dst(std::string("RandCsMatrix.") + typeid(T).name() + " copy",
          src.extent(0));
    Kokkos::deep_copy(dst, src);
    ExeSpaceType().fence();
    return dst;
  }

 public:
  std::string info;
  /// Constructs a random cs matrix.
  /// \param dim1 The first dimension: rows for Crs or columns for Ccs
  /// \param dim2 The second dimension: columns for Crs or rows for Ccs
  /// \param min_val The minimum scalar value in the matrix.
  /// \param max_val The maximum scalar value in the matrix.
  RandCsMatrix(int64_t dim1, int64_t dim2, ScalarType min_val,
               ScalarType max_val, bool fully_sparse = false) {
    __dim1         = dim1;
    __dim2         = dim2;
    __fully_sparse = fully_sparse;
    __map_d        = MapViewTypeD("RandCsMatrix.ColMapViewType", __dim1 + 1);
    __map          = Kokkos::create_mirror_view(__map_d);
    __ids_d        = IdViewTypeD("RandCsMatrix.RowIdViewType",
                          dim2 * dim1 + 1);  // over-allocated
    __ids          = Kokkos::create_mirror_view(__ids_d);

    uint64_t ticks =
        std::chrono::high_resolution_clock::now().time_since_epoch().count() %
        UINT32_MAX;

    info = std::string(
        std::string("RandCsMatrix<") + typeid(ScalarType).name() + ", " +
        typeid(LayoutType).name() + ", " + typeid(ExeSpaceType).name() + ">(" +
        std::to_string(dim2) + ", " + std::to_string(dim1) +
        "...): rand seed: " + std::to_string(ticks) +
        ", fully sparse: " + (__fully_sparse ? "true" : "false") + "\n");
    Kokkos::Random_XorShift64_Pool<Kokkos::HostSpace> random(ticks);
    __populate_random_cs_mat(ticks);

    __vals_d = ValViewTypeD("RandCsMatrix.ValViewType", __nnz + 1);
    __vals   = Kokkos::create_mirror_view(__vals_d);
    Kokkos::fill_random(__vals, random, min_val, max_val);  // random scalars
    Kokkos::fence();
    __vals(__nnz) = ScalarType(0);

    // Copy to device
    Kokkos::deep_copy(__vals_d, __vals);
    ExeSpaceType().fence();
  }

  // O(c), where c is a constant.
  ScalarType operator()(int64_t idx) { return __vals(idx); }
  int64_t get_nnz() { return __nnz; }
  // dimension2: This is either columns for a Crs matrix or rows for a Ccs
  // matrix.
  int64_t get_dim2() { return __dim2; }
  // dimension1: This is either rows for Crs matrix or columns for a Ccs matrix.
  int64_t get_dim1() { return __dim1; }
  ValViewTypeD get_vals() { return __getter_copy_helper(__vals_d); }
  IdViewTypeD get_ids() { return __getter_copy_helper(__ids_d); }
  MapViewTypeD get_map() { return __getter_copy_helper(__map_d); }
};

/// \brief Randomly shuffle the entries in each row (col) of a Crs (Ccs) matrix.
template <typename Rowptrs, typename Entries, typename Values>
void shuffleMatrixEntries(Rowptrs rowptrs, Entries entries, Values values) {
  using size_type    = typename Rowptrs::non_const_value_type;
  using ordinal_type = typename Entries::value_type;
  auto rowptrsHost =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rowptrs);
  auto entriesHost =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), entries);
  auto valuesHost =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), values);
  ordinal_type numRows =
      rowptrsHost.extent(0) ? (rowptrsHost.extent(0) - 1) : 0;
  for (ordinal_type i = 0; i < numRows; i++) {
    size_type rowBegin = rowptrsHost(i);
    size_type rowEnd   = rowptrsHost(i + 1);
    for (size_type j = rowBegin; j < rowEnd - 1; j++) {
      ordinal_type swapRange = rowEnd - j;
      size_type swapOffset   = j + (rand() % swapRange);
      std::swap(entriesHost(j), entriesHost(swapOffset));
      std::swap(valuesHost(j), valuesHost(swapOffset));
    }
  }
  Kokkos::deep_copy(entries, entriesHost);
  Kokkos::deep_copy(values, valuesHost);
}

}  // namespace Test
#endif
