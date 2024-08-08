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

#ifndef TEST_BLAS_SERIAL_NRM2_HPP_
#define TEST_BLAS_SERIAL_NRM2_HPP_

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"
#include "KokkosKernels_TestUtils.hpp"

#include "KokkosBlas1_nrm2.hpp"

namespace Test {

template <typename DeviceType, typename ViewType, typename AlgoTagType>
struct Functor_TestBlasSerialNrm2 {
  using execution_space = typename DeviceType::execution_space;
  using value_type      = typename ViewType::non_const_value_type;
  using IPT             = Kokkos::Details::InnerProductSpaceTraits<value_type>;
  using norm_type       = typename IPT::mag_type;
  using norm_view_type  = Kokkos::View<norm_type *, execution_space>;

  ViewType _x;
  norm_view_type _nrm;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBlasSerialNrm2(const ViewType &x, const norm_view_type &nrm) : _x(x), _nrm(nrm) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const KokkosKernelTag &, const int i) const {
    auto X  = Kokkos::subview(_x, i, Kokkos::ALL());
    _nrm(i) = KokkosBlas::serial_nrm2(X);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const NaiveTag &, const int k) const {
    auto X  = Kokkos::subview(_x, k, Kokkos::ALL());
    _nrm(k) = Kokkos::ArithTraits<norm_type>::zero();
    for (int i = 0; i < X.extent_int(0); ++i) {
      _nrm(k) += IPT::norm(IPT::dot(X(i), X(i)));
    }

    _nrm(k) = Kokkos::ArithTraits<norm_type>::sqrt(_nrm(k));
  }

  inline void run() {
    std::string name_region("KokkosBlas::Test::SerialNrm2");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name_work_tag         = (std::is_same<AlgoTagType, KokkosKernelTag>::value ? "::KokkosBlas"
                                         : std::is_same<AlgoTagType, NaiveTag>::value      ? "::Naive"
                                                                                           : "::UnknownWorkTag");
    std::string name_test_id          = "Nrm2";
    std::string name                  = name_region + name_value_type + name_work_tag + name_test_id;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space, AlgoTagType> policy(0, _x.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
    return;
  }
};

template <typename DeviceType, typename ViewType, typename AlgoTagType>
struct Functor_TestBlasSerialNrm2MV {
  using execution_space = typename DeviceType::execution_space;
  using value_type      = typename ViewType::non_const_value_type;
  using IPT             = Kokkos::Details::InnerProductSpaceTraits<value_type>;
  using norm_type       = typename IPT::mag_type;
  using norm_view_type  = Kokkos::View<norm_type **, execution_space>;

  ViewType _x;
  norm_view_type _nrm;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBlasSerialNrm2MV(const ViewType &x, const norm_view_type &nrm) : _x(x), _nrm(nrm) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const KokkosKernelTag &, const int i) const {
    auto X = Kokkos::subview(_x, i, Kokkos::ALL(), Kokkos::ALL());
    auto R = Kokkos::subview(_nrm, i, Kokkos::ALL());
    KokkosBlas::serial_nrm2(X, R);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const NaiveTag &, const int k) const {
    auto X = Kokkos::subview(_x, k, Kokkos::ALL(), Kokkos::ALL());
    auto R = Kokkos::subview(_nrm, k, Kokkos::ALL());

    for (int colIdx = 0; colIdx < X.extent_int(1); ++colIdx) {
      R(colIdx) = Kokkos::ArithTraits<norm_type>::zero();
      for (int rowIdx = 0; rowIdx < X.extent_int(0); ++rowIdx) {
        R(colIdx) += IPT::norm(IPT::dot(X(rowIdx, colIdx), X(rowIdx, colIdx)));
      }
      R(colIdx) = Kokkos::ArithTraits<norm_type>::sqrt(R(colIdx));
    }
  }

  inline void run() {
    std::string name_region("KokkosBlas::Test::SerialNrm2MV");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name_work_tag         = (std::is_same<AlgoTagType, KokkosKernelTag>::value ? "::KokkosBlas"
                                         : std::is_same<AlgoTagType, NaiveTag>::value      ? "::Naive"
                                                                                           : "::UnknownWorkTag");
    std::string name_test_id          = "Nrm2";
    std::string name                  = name_region + name_value_type + name_work_tag + name_test_id;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space, AlgoTagType> policy(0, _x.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
    return;
  }
};

template <typename DeviceType, typename ViewType>
void impl_test_blas_serial_nrm2(const int N, const int BlkSize) {
  /// typedefs
  using execution_space = typename DeviceType::execution_space;
  using value_type      = typename ViewType::non_const_value_type;
  using ats             = Kokkos::ArithTraits<value_type>;
  using IPT             = Kokkos::Details::InnerProductSpaceTraits<value_type>;
  using norm_type       = typename IPT::mag_type;
  using norm_view_type  = Kokkos::View<norm_type *, execution_space>;

  /// radomized input testing views
  ViewType X("X", N, BlkSize);
  Kokkos::Random_XorShift64_Pool<execution_space> random(13718);
  Kokkos::fill_random(X, random, ats::one());
  Kokkos::fence();

  norm_view_type norms("norms", N);
  norm_view_type norms_ref("ref norms", N);

  /// test body
  Functor_TestBlasSerialNrm2<DeviceType, ViewType, NaiveTag>(X, norms).run();
  Functor_TestBlasSerialNrm2<DeviceType, ViewType, KokkosKernelTag>(X, norms_ref).run();

  Kokkos::fence();

  /// for comparison send it to host
  typename norm_view_type::HostMirror norms_host     = Kokkos::create_mirror_view(norms);
  typename norm_view_type::HostMirror norms_ref_host = Kokkos::create_mirror_view(norms_ref);

  Kokkos::deep_copy(norms_host, norms);
  Kokkos::deep_copy(norms_ref_host, norms_ref);

  /// check a = b
  typename ats::mag_type eps = 100 * std::numeric_limits<typename ats::mag_type>::epsilon();
  for (int k = 0; k < N; ++k) EXPECT_NEAR_KK(norms_host(k), norms_ref_host(k), eps);
}

template <typename DeviceType, typename ViewType>
void impl_test_blas_serial_nrm2mv(const int N, const int vecLength, const int numVecs) {
  /// typedefs
  using execution_space = typename DeviceType::execution_space;
  using value_type      = typename ViewType::non_const_value_type;
  using ats             = Kokkos::ArithTraits<value_type>;
  using IPT             = Kokkos::Details::InnerProductSpaceTraits<value_type>;
  using norm_type       = typename IPT::mag_type;
  using norm_view_type  = Kokkos::View<norm_type **, execution_space>;

  /// radomized input testing views
  ViewType X("X", N, vecLength, numVecs);
  Kokkos::Random_XorShift64_Pool<execution_space> random(13718);
  Kokkos::fill_random(X, random, ats::one());
  Kokkos::fence();

  norm_view_type norms("norms", N, numVecs);
  norm_view_type norms_ref("ref norms", N, numVecs);

  /// test body
  Functor_TestBlasSerialNrm2MV<DeviceType, ViewType, NaiveTag>(X, norms).run();
  Functor_TestBlasSerialNrm2MV<DeviceType, ViewType, KokkosKernelTag>(X, norms_ref).run();

  Kokkos::fence();

  /// for comparison send it to host
  typename norm_view_type::HostMirror norms_host     = Kokkos::create_mirror_view(norms);
  typename norm_view_type::HostMirror norms_ref_host = Kokkos::create_mirror_view(norms_ref);

  Kokkos::deep_copy(norms_host, norms);
  Kokkos::deep_copy(norms_ref_host, norms_ref);

  /// check a = b
  typename ats::mag_type eps = 100 * std::numeric_limits<typename ats::mag_type>::epsilon();
  for (int k = 0; k < N; ++k)
    for (int vecIdx = 0; vecIdx < numVecs; ++vecIdx)
      EXPECT_NEAR_KK(norms_host(k, vecIdx), norms_ref_host(k, vecIdx), eps);
}

}  // namespace Test

template <typename DeviceType, typename ValueType>
int test_blas_serial_nrm2() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using ViewType = Kokkos::View<ValueType **, Kokkos::LayoutLeft, DeviceType>;
    Test::impl_test_blas_serial_nrm2<DeviceType, ViewType>(0, 10);
    Test::impl_test_blas_serial_nrm2<DeviceType, ViewType>(10, 15);
    Test::impl_test_blas_serial_nrm2<DeviceType, ViewType>(1024, 9);
    Test::impl_test_blas_serial_nrm2<DeviceType, ViewType>(132231, 3);

    using MVViewType = Kokkos::View<ValueType ***, Kokkos::LayoutLeft, DeviceType>;
    Test::impl_test_blas_serial_nrm2mv<DeviceType, MVViewType>(0, 10, 5);
    Test::impl_test_blas_serial_nrm2mv<DeviceType, MVViewType>(10, 15, 7);
    Test::impl_test_blas_serial_nrm2mv<DeviceType, MVViewType>(1024, 9, 5);
    Test::impl_test_blas_serial_nrm2mv<DeviceType, MVViewType>(132231, 3, 3);
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using ViewType = Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType>;
    Test::impl_test_blas_serial_nrm2<DeviceType, ViewType>(0, 10);
    Test::impl_test_blas_serial_nrm2<DeviceType, ViewType>(10, 15);
    Test::impl_test_blas_serial_nrm2<DeviceType, ViewType>(1024, 9);
    Test::impl_test_blas_serial_nrm2<DeviceType, ViewType>(132231, 3);

    using MVViewType = Kokkos::View<ValueType ***, Kokkos::LayoutLeft, DeviceType>;
    Test::impl_test_blas_serial_nrm2mv<DeviceType, MVViewType>(0, 10, 5);
    Test::impl_test_blas_serial_nrm2mv<DeviceType, MVViewType>(10, 15, 5);
    Test::impl_test_blas_serial_nrm2mv<DeviceType, MVViewType>(1024, 9, 5);
    Test::impl_test_blas_serial_nrm2mv<DeviceType, MVViewType>(132231, 3, 3);
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, serial_nrm2_float_float) { test_blas_serial_nrm2<TestDevice, float>(); }
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, serial_nrm2_double_double) { test_blas_serial_nrm2<TestDevice, double>(); }
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, serial_nrm2_fcomplex_float) { test_blas_serial_nrm2<TestDevice, Kokkos::complex<float> >(); }
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, serial_nrm2_dcomplex_dcomplex) { test_blas_serial_nrm2<TestDevice, Kokkos::complex<double> >(); }
#endif

#endif  // TEST_BLAS_SERIAL_NRM2_HPP_
