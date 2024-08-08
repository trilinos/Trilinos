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

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "KokkosBlas1_set.hpp"
#include "KokkosBlas1_scal.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBlas;

namespace Test {

enum : int { BlasSet = 0, BlasScale = 1 };

struct KokkosKernelTag {};
struct NaiveTag {};

template <typename DeviceType, typename ViewType, typename ScalarType, typename AlgoTagType, int TestID>
struct Functor_TestBlasSerialMatUtil {
  ScalarType _alpha;
  ViewType _a;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBlasSerialMatUtil(const ScalarType alpha, const ViewType &a) : _alpha(alpha), _a(a) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const KokkosKernelTag &, const int i) const {
    auto A = Kokkos::subview(_a, i, Kokkos::ALL(), Kokkos::ALL());
    switch (TestID) {
      case BlasSet: KokkosBlas::SerialSet::invoke(_alpha, A); break;
      case BlasScale: KokkosBlas::SerialScale::invoke(_alpha, A); break;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const NaiveTag &, const int k) const {
    // MD Note: changing because of the error with -werror
    auto A      = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
    const int m = A.extent(0), n = A.extent(1);
    switch (TestID) {
      case BlasSet: {
        for (int i = 0; i < m; ++i)
          for (int j = 0; j < n; ++j) A(i, j) = _alpha;
        break;
      }
      case BlasScale: {
        for (int i = 0; i < m; ++i)
          for (int j = 0; j < n; ++j) A(i, j) *= _alpha;
        break;
      }
    }
  }

  inline int run() {
    typedef typename ViewType::value_type value_type;
    std::string name_region("KokkosBlas::Test::SerialMatUtil");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name_work_tag         = (std::is_same<AlgoTagType, KokkosKernelTag>::value ? "::KokkosBlas"
                                         : std::is_same<AlgoTagType, NaiveTag>::value      ? "::Naive"
                                                                                           : "::UnknownWorkTag");
    std::string name_test_id          = (TestID == BlasSet ? "Set" : TestID == BlasScale ? "Scale" : "UnknownTest");
    std::string name                  = name_region + name_value_type + name_work_tag + name_test_id;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<typename DeviceType::execution_space, AlgoTagType> policy(0, _a.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
    return 0;
  }
};

template <typename DeviceType, typename ViewType, typename ScalarType, int TestID>
void impl_test_blas_matutil(const int N, const int BlkSize) {
  /// typedefs
  typedef typename ViewType::value_type value_type;
  typedef Kokkos::ArithTraits<value_type> ats;

  /// radomized input testing views
  const ScalarType alpha = 11.1;
  ViewType a("a", N, BlkSize, BlkSize);
  ViewType b("b", N, BlkSize, BlkSize);

  Kokkos::Random_XorShift64_Pool<typename DeviceType::execution_space> random(13718);
  Kokkos::fill_random(a, random, value_type(1.0));

  Kokkos::fence();

  Kokkos::deep_copy(b, a);

  /// test body
  Functor_TestBlasSerialMatUtil<DeviceType, ViewType, ScalarType, NaiveTag, TestID>(alpha, a).run();
  Functor_TestBlasSerialMatUtil<DeviceType, ViewType, ScalarType, KokkosKernelTag, TestID>(alpha, b).run();

  Kokkos::fence();

  /// for comparison send it to host
  typename ViewType::HostMirror a_host = Kokkos::create_mirror_view(a);
  typename ViewType::HostMirror b_host = Kokkos::create_mirror_view(b);

  Kokkos::deep_copy(a_host, a);
  Kokkos::deep_copy(b_host, b);

  /// check a = b
  typename ats::mag_type eps = 100 * std::numeric_limits<typename ats::mag_type>::epsilon();
  for (int k = 0; k < N; ++k)
    for (int i = 0; i < BlkSize; ++i)
      for (int j = 0; j < BlkSize; ++j) EXPECT_NEAR_KK(b_host(k, i, j), a_host(k, i, j), eps);
}
}  // namespace Test

template <typename DeviceType, typename ValueType, typename ScalarType, int TestID>
int test_blas_matutil() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutLeft, DeviceType> ViewType;
    Test::impl_test_blas_matutil<DeviceType, ViewType, ScalarType, TestID>(0, 10);
    Test::impl_test_blas_matutil<DeviceType, ViewType, ScalarType, TestID>(10, 15);
    Test::impl_test_blas_matutil<DeviceType, ViewType, ScalarType, TestID>(1024, 9);
    Test::impl_test_blas_matutil<DeviceType, ViewType, ScalarType, TestID>(132231, 3);
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutRight, DeviceType> ViewType;
    Test::impl_test_blas_matutil<DeviceType, ViewType, ScalarType, TestID>(0, 10);
    Test::impl_test_blas_matutil<DeviceType, ViewType, ScalarType, TestID>(10, 15);
    Test::impl_test_blas_matutil<DeviceType, ViewType, ScalarType, TestID>(1024, 9);
    Test::impl_test_blas_matutil<DeviceType, ViewType, ScalarType, TestID>(132231, 3);
  }
#endif

  return 0;
}

// Real test cases

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, blas_scalar_serial_set_float_float) {
  test_blas_matutil<TestDevice, float, float, ::Test::BlasSet>();
}
TEST_F(TestCategory, blas_scalar_serial_scale_float_float) {
  test_blas_matutil<TestDevice, float, float, ::Test::BlasScale>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, blas_scalar_serial_set_double_double) {
  test_blas_matutil<TestDevice, double, double, ::Test::BlasSet>();
}
TEST_F(TestCategory, blas_scalar_serial_scale_double_double) {
  test_blas_matutil<TestDevice, double, double, ::Test::BlasScale>();
}
#endif

// Complex test cases

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, blas_scalar_serial_set_dcomplex_dcomplex) {
  test_blas_matutil<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, ::Test::BlasSet>();
}
TEST_F(TestCategory, blas_scalar_serial_scale_dcomplex_dcomplex) {
  test_blas_matutil<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, ::Test::BlasScale>();
}
TEST_F(TestCategory, blas_scalar_serial_set_dcomplex_double) {
  test_blas_matutil<TestDevice, Kokkos::complex<double>, double, ::Test::BlasSet>();
}
TEST_F(TestCategory, blas_scalar_serial_scale_dcomplex_double) {
  test_blas_matutil<TestDevice, Kokkos::complex<double>, double, ::Test::BlasScale>();
}
#endif
