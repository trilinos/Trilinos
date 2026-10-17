// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#include <KokkosBlas1_rotg.hpp>

namespace Test {
template <class Device, class Scalar>
void test_rotg_impl(typename Device::execution_space const& space, Scalar const a_in, Scalar const b_in) {
  using magnitude_type = typename KokkosKernels::ArithTraits<Scalar>::mag_type;
  using SViewType      = Kokkos::View<Scalar, Device>;
  using MViewType      = Kokkos::View<magnitude_type, Device>;

  const magnitude_type eps = 10 * KokkosKernels::ArithTraits<Scalar>::eps();
  const Scalar zero        = KokkosKernels::ArithTraits<Scalar>::zero();

  // Initialize inputs/outputs
  SViewType a("a");
  Kokkos::deep_copy(a, a_in);
  SViewType b("b");
  Kokkos::deep_copy(b, b_in);
  MViewType c("c");
  SViewType s("s");

  KokkosBlas::rotg(space, a, b, c, s);
  space.fence();

  auto h_a = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, a);
  auto h_c = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, c);
  auto h_s = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, s);

  // The returned r must agree with applying the returned rotation.
  EXPECT_NEAR_KK(h_c() * a_in + h_s() * b_in, h_a(), eps);
  EXPECT_NEAR_KK(h_c() * b_in - KokkosKernels::ArithTraits<Scalar>::conj(h_s()) * a_in, zero, eps);
  EXPECT_NEAR_KK(h_c() * h_c() + Kokkos::abs(h_s()) * Kokkos::abs(h_s()), magnitude_type(1), eps);
}
}  // namespace Test

template <class Scalar, class Device>
int test_rotg() {
  const Scalar zero = KokkosKernels::ArithTraits<Scalar>::zero();
  const Scalar one  = KokkosKernels::ArithTraits<Scalar>::one();
  const Scalar two  = one + one;

  typename Device::execution_space space{};

  Test::test_rotg_impl<Device, Scalar>(space, one, zero);
  Test::test_rotg_impl<Device, Scalar>(space, zero, -one);
  Test::test_rotg_impl<Device, Scalar>(space, -zero, -one);
  Test::test_rotg_impl<Device, Scalar>(space, zero, one);
  Test::test_rotg_impl<Device, Scalar>(space, -one, zero);
  Test::test_rotg_impl<Device, Scalar>(space, zero, zero);
  Test::test_rotg_impl<Device, Scalar>(space, -one, two);
  Test::test_rotg_impl<Device, Scalar>(space, two, -one);
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
