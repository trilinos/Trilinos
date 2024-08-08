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
#include <KokkosBlas1_rotg.hpp>

namespace Test {
template <class Device, class Scalar>
void test_rotg_impl(typename Device::execution_space const& space, Scalar const a_in, Scalar const b_in) {
  using magnitude_type = typename Kokkos::ArithTraits<Scalar>::mag_type;
  using SViewType      = Kokkos::View<Scalar, Device>;
  using MViewType      = Kokkos::View<magnitude_type, Device>;

  // const magnitude_type eps = Kokkos::ArithTraits<Scalar>::eps();
  // const Scalar zero        = Kokkos::ArithTraits<Scalar>::zero();

  // Initialize inputs/outputs
  SViewType a("a");
  Kokkos::deep_copy(a, a_in);
  SViewType b("b");
  Kokkos::deep_copy(b, b_in);
  MViewType c("c");
  SViewType s("s");

  KokkosBlas::rotg(space, a, b, c, s);

  // Check that a*c - b*s == 0
  // and a == sqrt(a*a + b*b)
  // EXPECT_NEAR_KK(a_in * s - b_in * c, zero, 10 * eps);
  // EXPECT_NEAR_KK(Kokkos::sqrt(a_in * a_in + b_in * b_in), a, 10 * eps);
}
}  // namespace Test

template <class Scalar, class Device>
int test_rotg() {
  const Scalar zero = Kokkos::ArithTraits<Scalar>::zero();
  const Scalar one  = Kokkos::ArithTraits<Scalar>::one();
  const Scalar two  = one + one;

  typename Device::execution_space space{};

  Test::test_rotg_impl<Device, Scalar>(space, one, zero);
  Test::test_rotg_impl<Device, Scalar>(space, one / two, one / two);
  Test::test_rotg_impl<Device, Scalar>(space, 2.1 * one, 1.3 * one);

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, rotg_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::rotg");
  test_rotg<float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, rotg_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::rotg");
  test_rotg<double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, rotg_complex_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::rotg");
  test_rotg<Kokkos::complex<float>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, rotg_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::rotg");
  test_rotg<Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif
