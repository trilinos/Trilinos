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
#ifndef TEST_BLAS2_GEMV_UTIL_HPP
#define TEST_BLAS2_GEMV_UTIL_HPP

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Vector_SIMD.hpp>

namespace Test {

template <typename ValueType, int length = KokkosBatched::DefaultVectorLength<ValueType, TestDevice>::value>
using simd_vector = KokkosBatched::Vector<KokkosBatched::SIMD<ValueType>, length>;

template <class AType, class XType, class YType, class ScalarType>
struct GemvOpBase {
  GemvOpBase(char trans_, ScalarType alpha_, AType A_, XType x_, ScalarType beta_, YType y_)
      : trans(trans_), alpha(alpha_), beta(beta_), A(A_), x(x_), y(y_) {}

 protected:
  // parameters
  char trans;
  ScalarType alpha;
  ScalarType beta;
  // data
  AType A;
  XType x;
  YType y;
};

// Note: vanillaGEMV is called on device here - alternatively one can move
//       _strided_ data using safe_device_to_host_deep_copy() etc.
template <class AType, class XType, class YType, class ScalarType>
struct RefGEMVOp : public GemvOpBase<AType, XType, YType, ScalarType> {
  using params = GemvOpBase<AType, XType, YType, ScalarType>;

  RefGEMVOp(char trans_, ScalarType alpha_, AType A_, XType x_, ScalarType beta_, YType y_)
      : params(trans_, alpha_, A_, x_, beta_, y_) {}

  template <typename TeamMember>
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMember & /* member */) const {
    vanillaGEMV(params::trans, params::alpha, params::A, params::x, params::beta, params::y);
  }
};  // RefGEMVOp

// fill regular view with random values
template <class ViewType, class PoolType, class ScalarType = typename ViewType::non_const_value_type>
typename std::enable_if<!KokkosBatched::is_vector<ScalarType>::value>::type fill_random_view(
    ViewType A, PoolType &rand_pool, const ScalarType max_val = 10.0) {
  Kokkos::fill_random(A, rand_pool, max_val);
  Kokkos::fence();
}

// fill rank-1 view of SIMD vectors with random values
template <class ValueType, int VecLength, class Layout, class... Props, class PoolType>
void fill_random_view(
    Kokkos::View<KokkosBatched::Vector<KokkosBatched::SIMD<ValueType>, VecLength> *, Layout, Props...> x,
    PoolType &rand_pool, const ValueType max_val = 10.0) {
  // the view can be strided and have Vector<SIMD> values, so randoms
  // are generated in a plain, linear view first and then copied
  using device_type = typename decltype(x)::device_type;
  Kokkos::View<ValueType *, device_type> rnd("random_vals", x.extent(0) * VecLength);
  Kokkos::fill_random(rnd, rand_pool, max_val);
  using size_type = decltype(x.extent(0));
  for (size_type i = 0; i < x.extent(0); ++i) {
    x(i).loadUnaligned(&rnd(i * VecLength));
  }
}

// fill rank-2 view of SIMD vectors with random values
template <class ValueType, int VecLength, class Layout, class... Props, class PoolType>
static void fill_random_view(
    Kokkos::View<KokkosBatched::Vector<KokkosBatched::SIMD<ValueType>, VecLength> **, Layout, Props...> A,
    PoolType &rand_pool, const ValueType max_val = 10.0) {
  // the view can be strided and have Vector<SIMD> values, so randoms
  // are generated in a plain, linear view first and then copied
  using device_type = typename decltype(A)::device_type;
  Kokkos::View<ValueType *, device_type> rnd("random_vals", A.extent(0) * A.extent(1) * VecLength);
  Kokkos::fill_random(rnd, rand_pool, max_val);
  using size_type = decltype(A.extent(0));
  size_type idx   = 0;
  for (size_type i = 0; i < A.extent(0); ++i) {
    for (size_type j = 0; j < A.extent(1); ++j) {
      A(i, j).loadUnaligned(&rnd(idx));
      idx += VecLength;
    }
  }
}

template <class GemvFunc, class ScalarA, class ScalarX, class ScalarY, class Device, class ScalarCoef = void>
struct GEMVTest {
  static void run(const char *mode) { run_algorithms<0, typename GemvFunc::algorithms>(mode); }

 private:
  // ScalarCoef==void default behavior is to derive alpha/beta scalar types
  // from A and X scalar types
  using ScalarType = typename std::conditional<!std::is_void<ScalarCoef>::value, ScalarCoef,
                                               typename std::common_type<ScalarA, ScalarX>::type>::type;

  template <int Idx, class AlgorithmsTuple>
  static std::enable_if_t<Idx == std::tuple_size<AlgorithmsTuple>::value> run_algorithms(const char * /*mode*/) {}

  template <int Idx, class AlgorithmsTuple>
  static typename std::enable_if<(Idx < std::tuple_size<AlgorithmsTuple>::value)>::type run_algorithms(
      const char *mode) {
    run_layouts<typename std::tuple_element<Idx, AlgorithmsTuple>::type>(mode);
    run_algorithms<Idx + 1, AlgorithmsTuple>(mode);
  }

  // Note: all layouts listed here are subview'ed to test Kokkos::LayoutStride
  template <class AlgoTag>
  static void run_layouts(const char *mode) {
#ifdef KOKKOSKERNELS_TEST_LAYOUTLEFT
    run_view_types<AlgoTag, Kokkos::LayoutLeft>(mode);
#endif
#ifdef KOKKOSKERNELS_TEST_LAYOUTRIGHT
    run_view_types<AlgoTag, Kokkos::LayoutRight>(mode);
#endif
#if defined(KOKKOSKERNELS_TEST_LAYOUTLEFT) && defined(KOKKOSKERNELS_TEST_LAYOUTRIGHT)
    using A_t = typename Kokkos::View<ScalarA **, Kokkos::LayoutRight, Device>;
    using x_t = typename Kokkos::View<ScalarX *, Kokkos::LayoutLeft, Device>;
    using y_t = typename Kokkos::View<ScalarY *, Kokkos::LayoutRight, Device>;
    run_sizes<AlgoTag, A_t, x_t, y_t>(mode);
#endif
  }

  template <class AlgoTag, class Layout>
  static void run_view_types(const char *mode) {
    typedef Kokkos::View<ScalarA **, Layout, Device> view_type_A;
    typedef Kokkos::View<ScalarX *, Layout, Device> view_type_x;
    typedef Kokkos::View<ScalarY *, Layout, Device> view_type_y;
    run_sizes<AlgoTag, view_type_A, view_type_x, view_type_y>(mode);
  }

  template <class AlgoTag, class ViewAType, class ViewXType, class ViewYType>
  static void run_sizes(const char *mode) {
    // zero cases
    run_size<AlgoTag, ViewAType, ViewXType, ViewYType>(mode, 0, 0);
    run_size<AlgoTag, ViewAType, ViewXType, ViewYType>(mode, 0, 4);
    run_size<AlgoTag, ViewAType, ViewXType, ViewYType>(mode, 4, 0);
    // small block sizes
    for (int n = 1; n <= 16; ++n) {
      run_size<AlgoTag, ViewAType, ViewXType, ViewYType>(mode, n, n);
    }
    // other cases
    run_size<AlgoTag, ViewAType, ViewXType, ViewYType>(mode, 1024, 1);
    run_size<AlgoTag, ViewAType, ViewXType, ViewYType>(mode, 1024, 13);
    run_size<AlgoTag, ViewAType, ViewXType, ViewYType>(mode, 1024, 124);
  }

  template <class AlgoTag, class ViewTypeA, class ViewTypeX, class ViewTypeY>
  static void run_size(const char *mode, int N, int M) {
    using A_layout = typename ViewTypeA::array_layout;
    using x_layout = typename ViewTypeX::array_layout;
    using y_layout = typename ViewTypeY::array_layout;
    static_assert(!std::is_same<A_layout, Kokkos::LayoutStride>::value, "");
    static_assert(!std::is_same<x_layout, Kokkos::LayoutStride>::value, "");
    static_assert(!std::is_same<y_layout, Kokkos::LayoutStride>::value, "");

    const auto trans      = mode[0];
    const bool transposed = trans == (char)'T' || trans == (char)'C';
    const auto Nt         = transposed ? M : N;
    const auto Mt         = transposed ? N : M;

    // 1. run on regular (non-strided) views
    ViewTypeA A1("A1", Nt, Mt);
    ViewTypeX x1("X1", M);
    ViewTypeY y1("Y1", N);
    run_views<AlgoTag>(trans, A1, x1, y1);

    // 2. run on strided subviews (enforced by adding extra rank on both sides)
    // Note: strided views are not supported by MKL routines
    if (!std::is_same<AlgoTag, KokkosBlas::Algo::Gemv::CompactMKL>::value) {
      typedef Kokkos::View<ScalarA ****, A_layout, Device> BaseTypeA;
      typedef Kokkos::View<ScalarX ***, x_layout, Device> BaseTypeX;
      typedef Kokkos::View<ScalarY ***, y_layout, Device> BaseTypeY;

      BaseTypeA b_A("A", 2, Nt, Mt, 2);
      BaseTypeX b_x("X", 2, M, 2);
      BaseTypeY b_y("Y", 2, N, 2);
      auto A = Kokkos::subview(b_A, 0, Kokkos::ALL(), Kokkos::ALL(), 0);
      auto x = Kokkos::subview(b_x, 0, Kokkos::ALL(), 0);
      auto y = Kokkos::subview(b_y, 0, Kokkos::ALL(), 0);

      // make sure it's actually LayoutStride there
      static_assert(std::is_same<typename decltype(A)::array_layout, Kokkos::LayoutStride>::value, "");
      static_assert(std::is_same<typename decltype(x)::array_layout, Kokkos::LayoutStride>::value, "");
      static_assert(std::is_same<typename decltype(y)::array_layout, Kokkos::LayoutStride>::value, "");
      run_views<AlgoTag>(trans, A, x, y);
    }
  }

  template <class AlgoTag, class ViewTypeA, class ViewTypeX, class ViewTypeY>
  static void run_views(const char trans, ViewTypeA A, ViewTypeX x, ViewTypeY y) {
    Kokkos::TeamPolicy<typename Device::execution_space> teams(1, 1);  // just run on device
    fill_inputs(A, x, y);
    ScalarType alpha = 3;  // TODO: test also with zero alpha/beta ?
    ScalarType beta  = 5;

    // get reference results
    Kokkos::View<ScalarY *, Device> y_ref("Y_ref", y.extent(0));
    Kokkos::deep_copy(y_ref, y);
    RefGEMVOp<ViewTypeA, ViewTypeX, decltype(y_ref), ScalarType> gemv_ref(trans, alpha, A, x, beta, y_ref);
    Kokkos::parallel_for(teams, gemv_ref);

    // 1. check non-consts
    run_case<AlgoTag>(trans, alpha, A, x, beta, y, y_ref);

    // 2. check const x
    typename ViewTypeX::const_type c_x = x;
    run_case<AlgoTag>(trans, alpha, A, c_x, beta, y, y_ref);

    // 3. check const A and x
    typename ViewTypeA::const_type c_A = A;
    run_case<AlgoTag>(trans, alpha, c_A, c_x, beta, y, y_ref);
  }

  template <class AlgoTag, class ViewTypeA, class ViewTypeX, class ViewTypeY, class ViewTypeYRef, class ScalarType>
  static void run_case(const char trans, ScalarType alpha, ViewTypeA A, ViewTypeX x, ScalarType beta, ViewTypeY y,
                       ViewTypeYRef y_ref) {
    // run on original y view (not to alter the test)
    // but backup it and restore, so it can be reused
    Kokkos::View<ScalarY *, Device> y_backup("Y2", y.extent(0));
    Kokkos::deep_copy(y_backup, y);

    // fetch GEMV functor from the factory
    using op_type =
        typename GemvFunc::template functor_type<AlgoTag, ViewTypeA, ViewTypeX, ViewTypeY, Device, ScalarType>;

    op_type gemv_op(trans, alpha, A, x, beta, y);
    Kokkos::parallel_for(Kokkos::TeamPolicy<typename Device::execution_space>(1, 1), gemv_op);

    const double eps = epsilon(ScalarY{});
    EXPECT_NEAR_KK_REL_1DVIEW(y, y_ref, eps);
    Kokkos::deep_copy(y, y_backup);
  }

  //----- utilities -----//

  // GEMV tolerance for scalar types
  static double epsilon(float) { return 2 * 1e-5; }
  static double epsilon(double) { return 1e-7; }
  static double epsilon(int) { return 0; }
  // tolerance for derived types
  template <class ScalarType>
  static double epsilon(Kokkos::complex<ScalarType>) {
    return epsilon(ScalarType{});
  }
  template <class ScalarType, int VecLen>
  static double epsilon(simd_vector<ScalarType, VecLen>) {
    return epsilon(ScalarType{});
  }

  template <class ViewTypeA, class ViewTypeX, class ViewTypeY>
  static void fill_inputs(ViewTypeA A, ViewTypeX x, ViewTypeY y) {
    using exec_space = typename Device::execution_space;
    Kokkos::Random_XorShift64_Pool<exec_space> rand_pool(13718);
    fill_random_view(A, rand_pool);
    fill_random_view(x, rand_pool);
    fill_random_view(y, rand_pool);
  }
};  // struct GEMVTest

}  // namespace Test

#define TEST_CASE4(PREFIX, FACTORY, NAME, SCALAR_A, SCALAR_X, SCALAR_Y, SCALAR_COEF)            \
  using PREFIX##_##NAME##_gemv_test =                                                           \
      ::Test::GEMVTest<::Test::FACTORY, SCALAR_A, SCALAR_X, SCALAR_Y, TestDevice, SCALAR_COEF>; \
  TEST_F(TestCategory, PREFIX##_gemv_nt_##NAME) { PREFIX##_##NAME##_gemv_test::run("N"); }      \
  TEST_F(TestCategory, PREFIX##_gemv_t_##NAME) { PREFIX##_##NAME##_gemv_test::run("T"); }       \
  TEST_F(TestCategory, PREFIX##_gemv_ct_##NAME) { PREFIX##_##NAME##_gemv_test::run("C"); }

#define TEST_CASE2(PREFIX, FACTORY, NAME, SCALAR, SCALAR_COEF) \
  TEST_CASE4(PREFIX, FACTORY, NAME, SCALAR, SCALAR, SCALAR, SCALAR_COEF)
#define TEST_CASE(PREFIX, FACTORY, NAME, SCALAR) TEST_CASE2(PREFIX, FACTORY, NAME, SCALAR, SCALAR)

#endif  // TEST_BLAS2_GEMV_UTIL_HPP
