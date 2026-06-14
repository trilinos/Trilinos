// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Rotg.hpp>

namespace Test {
namespace Rotg {

template <typename DeviceType, typename SViewType, typename MViewType>
struct Functor_BatchedRotg {
  using execution_space = typename DeviceType::execution_space;
  SViewType m_a;
  SViewType m_b;
  MViewType m_c;
  SViewType m_s;

  Functor_BatchedRotg(const SViewType &a, const SViewType &b, const MViewType &c, const SViewType &s)
      : m_a(a), m_b(b), m_c(c), m_s(s) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    auto sub_a = Kokkos::subview(m_a, k);
    auto sub_b = Kokkos::subview(m_b, k);
    auto sub_c = Kokkos::subview(m_c, k);
    auto sub_s = Kokkos::subview(m_s, k);
    KokkosBatched::Rotg::invoke(sub_a, sub_b, sub_c, sub_s);
  }

  inline void run() {
    using value_type = typename SViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::Rotg");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space> policy(0, m_a.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename LayoutType, typename SType, typename MType>
void impl_test_batched_rotg_analytical(const std::size_t Nb, SType a_in, SType b_in, MType c_ref, SType s_ref) {
  using ats       = typename KokkosKernels::ArithTraits<SType>;
  using SViewType = Kokkos::View<SType *, LayoutType, DeviceType>;
  using MViewType = Kokkos::View<MType *, LayoutType, DeviceType>;

  SViewType a("a", Nb), b("b", Nb), s("s", Nb);
  MViewType c("c", Nb);

  Kokkos::deep_copy(a, a_in);
  Kokkos::deep_copy(b, b_in);

  Functor_BatchedRotg<DeviceType, SViewType, MViewType>(a, b, c, s).run();

  auto h_a = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, a);
  auto h_b = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, b);
  auto h_c = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, c);
  auto h_s = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, s);

  typename ats::mag_type eps = 1.0e1 * ats::epsilon();
  for (std::size_t i = 0; i < Nb; i++) {
    // Check if the computed c and s satisfy the Givens rotation properties
    EXPECT_NEAR_KK(h_c(i), c_ref, eps);
    EXPECT_NEAR_KK(h_s(i), s_ref, eps);
    auto rotation_norm = h_c(i) * h_c(i) + Kokkos::abs(h_s(i)) * Kokkos::abs(h_s(i));
    EXPECT_NEAR_KK(rotation_norm, 1.0, eps);

    // Check if the rotated values satisfy the expected relationships
    // For the given a_in and b_in, the expected rotated values are:
    // r = c*a_in + s*b_in should be approximately equal to sqrt(a_in^2 + b_in^2)
    // z = c*b_in - conj(s)*a should be approximately zero
    using Op = std::conditional_t<KokkosKernels::ArithTraits<SType>::is_complex, KokkosBlas::Impl::OpConj,
                                  KokkosBlas::Impl::OpID>;
    Op op;
    auto r = h_c(i) * a_in + h_s(i) * b_in;
    auto z = h_c(i) * b_in - op(h_s(i)) * a_in;

    SType r_ref = Kokkos::sqrt(Kokkos::abs(a_in) * Kokkos::abs(a_in) + Kokkos::abs(b_in) * Kokkos::abs(b_in));
    if constexpr (KokkosKernels::ArithTraits<SType>::is_complex) {
      if (Kokkos::abs(a_in) == KokkosKernels::ArithTraits<SType>::zero()) {
        r_ref = Kokkos::abs(b_in);
      } else {
        r_ref *= a_in / Kokkos::abs(a_in);
      }
    } else {
      auto sign = (Kokkos::abs(a_in) > Kokkos::abs(b_in)) ? a_in : b_in;
      r_ref     = Kokkos::copysign(r_ref, sign);
    }
    EXPECT_NEAR_KK(r, r_ref, eps);
    EXPECT_NEAR_KK(z, 0.0, eps);

    auto a_ref = r_ref;
    auto b_ref = b_in;

    // b is updated only for the real case
    if constexpr (!KokkosKernels::ArithTraits<SType>::is_complex) {
      if (Kokkos::abs(a_in) > Kokkos::abs(b_in)) {
        b_ref = h_s(i);
      } else if (h_c(i) != 0) {
        b_ref = KokkosKernels::ArithTraits<SType>::one() / h_c(i);
      } else {
        b_ref = KokkosKernels::ArithTraits<SType>::one();
      }
    }

    EXPECT_NEAR_KK(h_a(i), a_ref, eps);
    EXPECT_NEAR_KK(h_b(i), b_ref, eps);
  }
}

template <typename DeviceType, typename LayoutType, typename SType, typename MType>
void impl_test_batched_rotg_random(const std::size_t Nb) {
  using ats       = typename KokkosKernels::ArithTraits<SType>;
  using SViewType = Kokkos::View<SType *, LayoutType, DeviceType>;
  using MViewType = Kokkos::View<MType *, LayoutType, DeviceType>;

  SViewType a("a", Nb), b("b", Nb), s("s", Nb), a_init("a_init", Nb), b_init("b_init", Nb);
  MViewType c("c", Nb);

  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  SType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(a_init, rand_pool, randStart, randEnd);
  Kokkos::fill_random(b_init, rand_pool, randStart, randEnd);

  Kokkos::deep_copy(a, a_init);
  Kokkos::deep_copy(b, b_init);

  Functor_BatchedRotg<DeviceType, SViewType, MViewType>(a, b, c, s).run();

  auto h_a      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, a);
  auto h_b      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, b);
  auto h_c      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, c);
  auto h_s      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, s);
  auto h_a_init = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, a_init);
  auto h_b_init = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, b_init);

  typename ats::mag_type eps = 1.0e1 * ats::epsilon();
  for (std::size_t i = 0; i < Nb; i++) {
    // Check if the computed c and s satisfy the Givens rotation properties
    auto rotation_norm = h_c(i) * h_c(i) + Kokkos::abs(h_s(i)) * Kokkos::abs(h_s(i));
    EXPECT_NEAR_KK(rotation_norm, 1.0, eps);

    // Check if the rotated values satisfy the expected relationships
    // For the given a_in and b_in, the expected rotated values are:
    // r = c*a_in + s*b_in should be approximately equal to sqrt(a_in^2 + b_in^2)
    // z = c*b_in - conj(s)*a should be approximately zero
    using Op = std::conditional_t<KokkosKernels::ArithTraits<SType>::is_complex, KokkosBlas::Impl::OpConj,
                                  KokkosBlas::Impl::OpID>;
    Op op;
    auto r      = h_c(i) * h_a_init(i) + h_s(i) * h_b_init(i);
    auto z      = h_c(i) * h_b_init(i) - op(h_s(i)) * h_a_init(i);
    SType r_ref = Kokkos::sqrt(Kokkos::abs(h_a_init(i)) * Kokkos::abs(h_a_init(i)) +
                               Kokkos::abs(h_b_init(i)) * Kokkos::abs(h_b_init(i)));
    if constexpr (KokkosKernels::ArithTraits<SType>::is_complex) {
      if (Kokkos::abs(h_a_init(i)) == KokkosKernels::ArithTraits<SType>::zero()) {
        r_ref = Kokkos::abs(h_b_init(i));
      } else {
        r_ref *= h_a_init(i) / Kokkos::abs(h_a_init(i));
      }
    } else {
      auto sign = (Kokkos::abs(h_a_init(i)) > Kokkos::abs(h_b_init(i))) ? h_a_init(i) : h_b_init(i);
      r_ref     = Kokkos::copysign(r_ref, sign);
    }
    EXPECT_NEAR_KK(r, r_ref, eps);
    EXPECT_NEAR_KK(z, 0.0, eps);

    auto a_ref = r_ref;
    auto b_ref = h_b_init(i);

    // b is updated only for the real case
    if constexpr (!KokkosKernels::ArithTraits<SType>::is_complex) {
      if (Kokkos::abs(h_a_init(i)) > Kokkos::abs(h_b_init(i))) {
        b_ref = h_s(i);
      } else if (h_c(i) != 0) {
        b_ref = KokkosKernels::ArithTraits<SType>::one() / h_c(i);
      } else {
        b_ref = KokkosKernels::ArithTraits<SType>::one();
      }
    }

    EXPECT_NEAR_KK(h_a(i), a_ref, eps);
    EXPECT_NEAR_KK(h_b(i), b_ref, eps);
  }
}

template <typename DeviceType, typename SType, typename LayoutType>
void impl_test_batched_rotg(const std::size_t Nb) {
  using MType = typename KokkosKernels::ArithTraits<SType>::mag_type;
  // Test with analytical values
  constexpr auto pi = Kokkos::numbers::pi_v<MType>;
  MType c_ref       = Kokkos::cos(pi / 6.0);
  SType s_ref       = Kokkos::sin(pi / 6.0);
  SType r           = 2.5;  // scaling factor
  if constexpr (KokkosKernels::ArithTraits<SType>::is_complex) {
    r = SType(2.5, 2.0);
  }
  SType a_in = r * c_ref;
  SType b_in = r * s_ref;
  impl_test_batched_rotg_analytical<DeviceType, LayoutType, SType, MType>(Nb, a_in, b_in, c_ref, s_ref);

  // Test with a == 0
  a_in  = SType(0);
  b_in  = 3.0;
  c_ref = 0.0;
  s_ref = 1.0;
  if constexpr (KokkosKernels::ArithTraits<SType>::is_complex) {
    b_in  = SType(3.0, 4.0);
    s_ref = SType(0.6, -0.8);  // s_ref should be conj(b_in) / |b_in| when a == 0
  }
  impl_test_batched_rotg_analytical<DeviceType, LayoutType, SType, MType>(Nb, a_in, b_in, c_ref, s_ref);

  // Test with b == 0
  a_in  = 4.0;
  b_in  = 0.0;
  c_ref = 1.0;
  s_ref = 0.0;
  if constexpr (KokkosKernels::ArithTraits<SType>::is_complex) {
    a_in = SType(4.0, 5.0);
  }
  impl_test_batched_rotg_analytical<DeviceType, LayoutType, SType, MType>(Nb, a_in, b_in, c_ref, s_ref);

  // Test with random values
  impl_test_batched_rotg_random<DeviceType, LayoutType, SType, MType>(Nb);
}

}  // namespace Rotg
}  // namespace Test

template <typename DeviceType, typename ScalarType>
int test_batched_rotg() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    for (int i = 1; i < 3; i++) {
      Test::Rotg::impl_test_batched_rotg<DeviceType, ScalarType, LayoutType>(i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    for (int i = 1; i < 3; i++) {
      Test::Rotg::impl_test_batched_rotg<DeviceType, ScalarType, LayoutType>(i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_rotg_float) { test_batched_rotg<TestDevice, float>(); }
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_rotg_double) { test_batched_rotg<TestDevice, double>(); }
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, test_batched_rotg_fcomplex) { test_batched_rotg<TestDevice, Kokkos::complex<float>>(); }
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_rotg_dcomplex) { test_batched_rotg<TestDevice, Kokkos::complex<double>>(); }
#endif
