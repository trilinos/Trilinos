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
/// \author Kim Liegeois (knliege@sandia.gov)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "KokkosBatched_Axpy.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched;

namespace Test {
namespace Axpy {

template <typename DeviceType, typename ViewType, typename alphaViewType>
struct Functor_TestBatchedSerialAxpy {
  using execution_space = typename DeviceType::execution_space;
  const alphaViewType _alpha;
  const ViewType _X;
  const ViewType _Y;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedSerialAxpy(const alphaViewType &alpha, const ViewType &X, const ViewType &Y)
      : _alpha(alpha), _X(X), _Y(Y) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    auto alpha = Kokkos::subview(_alpha, Kokkos::make_pair(k, k + 1));
    auto x     = Kokkos::subview(_X, Kokkos::make_pair(k, k + 1), Kokkos::ALL);
    auto y     = Kokkos::subview(_Y, Kokkos::make_pair(k, k + 1), Kokkos::ALL);

    KokkosBatched::SerialAxpy::invoke(alpha, x, y);
  }

  inline void run() {
    typedef typename ViewType::value_type value_type;
    std::string name_region("KokkosBatched::Test::SerialAxpy");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space> policy(0, _X.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename ViewType, typename alphaViewType>
void impl_test_batched_axpy(const int N, const int BlkSize) {
  typedef typename ViewType::value_type value_type;
  typedef typename ViewType::const_value_type const_value_type;
  typedef typename alphaViewType::const_value_type alpha_const_value_type;
  typedef Kokkos::ArithTraits<value_type> ats;

  ViewType X0("x0", N, BlkSize), X1("x1", N, BlkSize), Y0("y0", N, BlkSize), Y1("y1", N, BlkSize);

  alphaViewType alpha("alpha", N);

  Kokkos::Random_XorShift64_Pool<typename DeviceType::execution_space> random(13718);
  Kokkos::fill_random(X0, random, const_value_type(1.0));
  Kokkos::fill_random(Y0, random, const_value_type(1.0));
  Kokkos::fill_random(alpha, random, alpha_const_value_type(1.0));

  Kokkos::fence();

  Kokkos::deep_copy(X1, X0);
  Kokkos::deep_copy(Y1, Y0);

  /// test body
  auto alpha_host = Kokkos::create_mirror_view(alpha);
  auto X0_host    = Kokkos::create_mirror_view(X0);
  auto Y0_host    = Kokkos::create_mirror_view(Y0);

  Kokkos::deep_copy(alpha_host, alpha);
  Kokkos::deep_copy(X0_host, X0);
  Kokkos::deep_copy(Y0_host, Y0);

  for (int l = 0; l < N; ++l)
    for (int i = 0; i < BlkSize; ++i) Y0_host(l, i) += alpha_host(l) * X0_host(l, i);

  Functor_TestBatchedSerialAxpy<DeviceType, ViewType, alphaViewType>(alpha, X1, Y1).run();

  Kokkos::fence();

  /// for comparison send it to host
  auto Y1_host = Kokkos::create_mirror_view(Y1);

  Kokkos::deep_copy(Y1_host, Y1);

  /// check c0 = c1 ; this eps is about 10^-14
  typedef typename ats::mag_type mag_type;
  mag_type sum(1), diff(0);
  const mag_type eps = 1.0e3 * ats::epsilon();

  for (int l = 0; l < N; ++l)
    for (int i = 0; i < BlkSize; ++i) {
      sum += ats::abs(Y0_host(l, i));
      diff += ats::abs(Y0_host(l, i) - Y1_host(l, i));
    }
  EXPECT_NEAR_KK(diff / sum, 0, eps);
}
}  // namespace Axpy
}  // namespace Test

template <typename DeviceType, typename ValueType, typename ScalarType>
int test_batched_axpy() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    typedef Kokkos::View<ValueType **, Kokkos::LayoutLeft, DeviceType> ViewType;
    typedef Kokkos::View<ScalarType *, Kokkos::LayoutLeft, DeviceType> alphaViewType;

    for (int i = 3; i < 10; ++i) {
      Test::Axpy::impl_test_batched_axpy<DeviceType, ViewType, alphaViewType>(1024, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    typedef Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType> ViewType;
    typedef Kokkos::View<ScalarType *, Kokkos::LayoutRight, DeviceType> alphaViewType;

    for (int i = 3; i < 10; ++i) {
      Test::Axpy::impl_test_batched_axpy<DeviceType, ViewType, alphaViewType>(1024, i);
    }
  }
#endif

  return 0;
}
