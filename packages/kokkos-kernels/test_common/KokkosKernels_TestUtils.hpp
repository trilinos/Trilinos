/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSKERNELS_TEST_UTILS_HPP
#define KOKKOSKERNELS_TEST_UTILS_HPP

#include "KokkosKernels_Utils.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "KokkosSparse_spmv.hpp"
// Make this include-able from all subdirectories
#include "../tpls/gtest/gtest/gtest.h"  //for EXPECT_**

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

template <class ViewTypeA, class ViewTypeX, class ViewTypeY>
void vanillaGEMV(char mode, typename ViewTypeA::non_const_value_type alpha,
                 const ViewTypeA& A, const ViewTypeX& x,
                 typename ViewTypeY::non_const_value_type beta,
                 const ViewTypeY& y) {
  using ScalarY = typename ViewTypeY::non_const_value_type;
  using KAT_A   = Kokkos::ArithTraits<typename ViewTypeA::non_const_value_type>;
  using KAT_Y   = Kokkos::ArithTraits<ScalarY>;
  int M         = A.extent(0);
  int N         = A.extent(1);
  if (beta == KAT_Y::zero()) Kokkos::deep_copy(y, KAT_Y::zero());
  if (mode == 'N') {
    for (int i = 0; i < M; i++) {
      ScalarY y_i = beta * y(i);
      for (int j = 0; j < N; j++) {
        y_i += alpha * A(i, j) * x(j);
      }
      y(i) = y_i;
    }
  } else if (mode == 'T') {
    for (int j = 0; j < N; j++) {
      ScalarY y_j = beta * y(j);
      for (int i = 0; i < M; i++) {
        y_j += alpha * A(i, j) * x(i);
      }
      y(j) = y_j;
    }
  } else if (mode == 'C') {
    for (int j = 0; j < N; j++) {
      ScalarY y_j = beta * y(j);
      for (int i = 0; i < M; i++) {
        y_j += alpha * KAT_A::conj(A(i, j)) * x(i);
      }
      y(j) = y_j;
    }
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

template <typename crsMat_t, typename vector_t>
vector_t create_random_y_vector(crsMat_t crsMat, vector_t x_vector) {
  vector_t y_vector(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Y VECTOR"),
                    crsMat.numRows());
  KokkosSparse::spmv("N", 1, crsMat, x_vector, 0, y_vector);
  return y_vector;
}

template <typename crsMat_t, typename vector_t>
vector_t create_random_y_vector_mv(crsMat_t crsMat, vector_t x_vector) {
  vector_t y_vector(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Y VECTOR"),
                    crsMat.numRows(), x_vector.extent(1));
  KokkosSparse::spmv("N", 1, crsMat, x_vector, 0, y_vector);
  return y_vector;
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

}  // namespace Test
#endif
