// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Rotmg.hpp>

namespace Test {
namespace Rotmg {

template <typename DeviceType, typename DXViewType, typename YViewType, typename PViewType>
struct Functor_BatchedRotmg {
  using execution_space = typename DeviceType::execution_space;
  DXViewType m_d1;
  DXViewType m_d2;
  DXViewType m_x1;
  YViewType m_y1;
  PViewType m_param;

  Functor_BatchedRotmg(const DXViewType &d1, const DXViewType &d2, const DXViewType &x1, const YViewType &y1,
                       const PViewType &param)
      : m_d1(d1), m_d2(d2), m_x1(x1), m_y1(y1), m_param(param) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    auto sub_d1    = Kokkos::subview(m_d1, k);
    auto sub_d2    = Kokkos::subview(m_d2, k);
    auto sub_x1    = Kokkos::subview(m_x1, k);
    auto sub_y1    = Kokkos::subview(m_y1, k);
    auto sub_param = Kokkos::subview(m_param, k, Kokkos::ALL());
    KokkosBatched::Rotmg::invoke(sub_d1, sub_d2, sub_x1, sub_y1, sub_param);
  }

  inline void run() {
    using value_type = typename DXViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::Rotmg");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space> policy(0, m_d1.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename LayoutType, typename DXType, typename YType, typename PType>
void impl_test_batched_rotmg_analytical(const std::size_t Nb, DXType d1_in, DXType d2_in, DXType x1_in, YType y1_in,
                                        const DXType ref[8]) {
  using ats        = typename KokkosKernels::ArithTraits<DXType>;
  using DXViewType = Kokkos::View<DXType *, LayoutType, DeviceType>;
  using YViewType  = Kokkos::View<YType *, LayoutType, DeviceType>;
  using PViewType  = Kokkos::View<PType **, LayoutType, DeviceType>;

  DXViewType d1("d1", Nb), d2("d2", Nb), x1("x1", Nb);
  YViewType y1("y1", Nb);
  PViewType param("param", Nb, 5);

  Kokkos::deep_copy(d1, d1_in);
  Kokkos::deep_copy(d2, d2_in);
  Kokkos::deep_copy(x1, x1_in);
  Kokkos::deep_copy(y1, y1_in);

  Functor_BatchedRotmg<DeviceType, DXViewType, YViewType, PViewType>(d1, d2, x1, y1, param).run();

  auto h_d1    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, d1);
  auto h_d2    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, d2);
  auto h_x1    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x1);
  auto h_y1    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, y1);
  auto h_param = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, param);

  typename ats::mag_type eps = 1.0e1 * ats::epsilon();
  for (std::size_t i = 0; i < Nb; i++) {
    EXPECT_NEAR_KK_REL(h_d1(i), ref[0], eps, "rotmg: d1 does not match reference");
    EXPECT_NEAR_KK_REL(h_d2(i), ref[1], eps, "rotmg: d2 does not match reference");
    EXPECT_NEAR_KK_REL(h_x1(i), ref[2], eps, "rotmg: x1 does not match reference");
    EXPECT_NEAR_KK_REL(h_param(i, 0), ref[3], eps, "rotmg: param[0] (flag) does not match reference");
    EXPECT_NEAR_KK_REL(h_param(i, 1), ref[4], eps, "rotmg: param[1] (h11) does not match reference");
    EXPECT_NEAR_KK_REL(h_param(i, 2), ref[5], eps, "rotmg: param[2] (h21) does not match reference");
    EXPECT_NEAR_KK_REL(h_param(i, 3), ref[6], eps, "rotmg: param[3] (h12) does not match reference");
    EXPECT_NEAR_KK_REL(h_param(i, 4), ref[7], eps, "rotmg: param[4] (h22) does not match reference");
  }
}

template <typename DeviceType, typename SType, typename LayoutType>
void impl_test_batched_rotmg(const std::size_t Nb) {
  // Test with analytical values
  SType d1_in, d2_in, x1_in, y1_in;
  SType ref[8];  // d1, d2, x1, flag, h11, h21, h12, h22

  // flag == -2 case (p2 == 0 branch)
  // where p2 = d2 * y1^2
  d1_in  = SType(1.0);
  d2_in  = SType(2.0);
  x1_in  = SType(3.0);
  y1_in  = SType(0.0);
  ref[0] = SType(1.0);   // d1
  ref[1] = SType(2.0);   // d2
  ref[2] = SType(3.0);   // x1
  ref[3] = SType(-2.0);  // flag
  ref[4] = SType(0.0);   // h11
  ref[5] = SType(0.0);   // h21
  ref[6] = SType(0.0);   // h12
  ref[7] = SType(0.0);   // h22
  impl_test_batched_rotmg_analytical<DeviceType, LayoutType, SType, SType, SType>(Nb, d1_in, d2_in, x1_in, y1_in, ref);

  // flag == -1 case (d1 < 0 branch)
  d1_in  = SType(-1.0);
  d2_in  = SType(2.0);
  x1_in  = SType(3.0);
  y1_in  = SType(4.0);
  ref[0] = SType(0.0);   // d1
  ref[1] = SType(0.0);   // d2
  ref[2] = SType(0.0);   // x1
  ref[3] = SType(-1.0);  // flag
  ref[4] = SType(0.0);   // h11
  ref[5] = SType(0.0);   // h21
  ref[6] = SType(0.0);   // h12
  ref[7] = SType(0.0);   // h22
  impl_test_batched_rotmg_analytical<DeviceType, LayoutType, SType, SType, SType>(Nb, d1_in, d2_in, x1_in, y1_in, ref);

  // flag == 0 case (abs(q1) > abs(q2) branch)
  // where q1 = d1 * x1^2 and q2 = d2 * y1^2
  d1_in = SType(10.0);
  d2_in = SType(1.0);
  x1_in = SType(5.0);
  y1_in = SType(1.0);

  auto p1  = d1_in * x1_in;
  auto p2  = d2_in * y1_in;
  auto h21 = -y1_in / x1_in;
  auto h12 = p2 / p1;
  auto u   = 1.0 - h12 * h21;
  ref[0]   = d1_in / u;   // d1
  ref[1]   = d2_in / u;   // d2
  ref[2]   = x1_in * u;   // x1
  ref[3]   = SType(0.0);  // flag
  ref[4]   = SType(0.0);  // h11
  ref[5]   = h21;         // h21
  ref[6]   = h12;         // h12
  ref[7]   = SType(0.0);  // h22
  impl_test_batched_rotmg_analytical<DeviceType, LayoutType, SType, SType, SType>(Nb, d1_in, d2_in, x1_in, y1_in, ref);

  // flag == 1 case (abs(q1) <= abs(q2) branch)
  // where q1 = d1 * x1^2 and q2 = d2 * y1^2
  d1_in    = SType(1.0);
  d2_in    = SType(10.0);
  x1_in    = SType(1.0);
  y1_in    = SType(5.0);
  p1       = d1_in * x1_in;
  p2       = d2_in * y1_in;
  auto h11 = p1 / p2;
  auto h22 = x1_in / y1_in;
  u        = 1.0 + h11 * h22;
  ref[0]   = d2_in / u;   // d1
  ref[1]   = d1_in / u;   // d2
  ref[2]   = y1_in * u;   // x1
  ref[3]   = SType(1.0);  // flag
  ref[4]   = h11;         // h11
  ref[5]   = SType(0.0);  // h21
  ref[6]   = SType(0.0);  // h12
  ref[7]   = h22;         // h22
  impl_test_batched_rotmg_analytical<DeviceType, LayoutType, SType, SType, SType>(Nb, d1_in, d2_in, x1_in, y1_in, ref);

  // small d1 case (|q1| > |q2| with d1 scaled up by gam^2)
  constexpr SType gam = SType(4096.0);
  d1_in               = SType(6.0e-10);
  d2_in               = SType(2.0e-2);
  x1_in               = SType(1.0e5);
  y1_in               = SType(10.0);
  p1                  = d1_in * x1_in;
  p2                  = d2_in * y1_in;
  h21                 = -1.0e-4;                  // -y1_in / x1_in
  h12                 = 1.0 / 3.0 * 1.0e4;        // (d2_in * y1_in) / (d1_in * x1_in)
  u                   = 4.0 / 3.0;                // 1.0 - h12 * h21
  ref[0]              = 4.5e-10 * gam * gam;      // (d1_in / u) * (gam * gam);
  ref[1]              = 1.5e-2;                   // (d2_in / u)
  ref[2]              = 4.0 / 3.0 * 1.0e5 / gam;  // (x1 * u / gam)
  ref[3]              = SType(-1.0);              // flag
  ref[4]              = 1.0 / gam;                // h11
  ref[5]              = h21;                      // h21
  ref[6]              = 1.0e4 / (3.0 * gam);      // h12 / gam
  ref[7]              = SType(1.0);               // h22

  impl_test_batched_rotmg_analytical<DeviceType, LayoutType, SType, SType, SType>(Nb, d1_in, d2_in, x1_in, y1_in, ref);

  // big d1 case (|q1| > |q2| with d1 scaled down by gam^2)
  d1_in  = SType(6.0e10);
  d2_in  = SType(2.0e-2);
  x1_in  = SType(1.0e-5);
  y1_in  = SType(10.0);
  p1     = d1_in * x1_in;
  p2     = d2_in * y1_in;
  h21    = -1.0e6;                    // -y1_in / x1_in
  h12    = 1.0 / 3.0 * 1.0e-6;        // (d2_in * y1_in) / (d1_in * x1_in)
  u      = 4.0 / 3.0;                 // 1.0 - h12 * h21
  ref[0] = 4.5e10 / (gam * gam);      // (d1_in / u) / (gam * gam);
  ref[1] = 1.5e-2;                    // (d2_in / u)
  ref[2] = 4.0 / 3.0 * 1.0e-5 * gam;  // (x1 * u * gam)
  ref[3] = SType(-1.0);               // flag
  ref[4] = 1.0 * gam;                 // h11
  ref[5] = h21;                       // h21
  ref[6] = h12 * gam;                 // h12 * gam
  ref[7] = SType(1.0);                // h22

  impl_test_batched_rotmg_analytical<DeviceType, LayoutType, SType, SType, SType>(Nb, d1_in, d2_in, x1_in, y1_in, ref);

  // small d1 case (|q1| <= |q2| with d1 scaled up by gam^2)
  d1_in  = SType(2.0e-10);
  d2_in  = SType(4.0e-2);
  x1_in  = SType(1.0e5);
  y1_in  = SType(10.0);
  p1     = d1_in * x1_in;
  p2     = d2_in * y1_in;
  h11    = 5.0e-5;                             // (d1_in * x1_in) / (d2_in * y1_in)
  h22    = 1.0e4;                              // (x1_in / y1_in)
  u      = 1.5;                                // 1.0 + h11 * h22
  ref[0] = 8.0 / 3.0 * 1.0e-2;                 // (d2_in / u)
  ref[1] = 4.0 / 3.0 * 1.0e-10 * (gam * gam);  // (d1_in / u) * (gam * gam);
  ref[2] = 15;                                 // (x1 * u)
  ref[3] = SType(-1.0);                        // flag
  ref[4] = h11;                                // h11
  ref[5] = SType(-1.0) / gam;                  // 1.0 / gam
  ref[6] = SType(1.0);                         // h12
  ref[7] = h22 / gam;                          // h22 / gam

  impl_test_batched_rotmg_analytical<DeviceType, LayoutType, SType, SType, SType>(Nb, d1_in, d2_in, x1_in, y1_in, ref);

  // big d1 case (|q1| <= |q2| with d1 scaled down by gam^2)
  d1_in  = SType(2.0e10);
  d2_in  = SType(4.0e-2);
  x1_in  = SType(1.0e-5);
  y1_in  = SType(10.0);
  p1     = d1_in * x1_in;
  p2     = d2_in * y1_in;
  h11    = 5.0e5;                             // (d1_in * x1_in) / (d2_in * y1_in)
  h22    = 1.0e-6;                            // (x1_in / y1_in)
  u      = 1.5;                               // 1.0 + h11 * h22
  ref[0] = 8.0 / 3.0 * 1.0e-2;                // (d2_in / u)
  ref[1] = 4.0 / 3.0 * 1.0e10 / (gam * gam);  // (d1_in / u) / (gam * gam);
  ref[2] = 15;                                // (x1 * u)
  ref[3] = SType(-1.0);                       // flag
  ref[4] = h11;                               // h11
  ref[5] = -gam;                              // h21
  ref[6] = SType(1.0);                        // h12
  ref[7] = h22 * gam;                         // h22 / gam

  impl_test_batched_rotmg_analytical<DeviceType, LayoutType, SType, SType, SType>(Nb, d1_in, d2_in, x1_in, y1_in, ref);
}

}  // namespace Rotmg
}  // namespace Test

template <typename DeviceType, typename ScalarType>
int test_batched_rotmg() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    for (int i = 1; i < 3; i++) {
      Test::Rotmg::impl_test_batched_rotmg<DeviceType, ScalarType, LayoutType>(i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    for (int i = 1; i < 3; i++) {
      Test::Rotmg::impl_test_batched_rotmg<DeviceType, ScalarType, LayoutType>(i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_rotmg_float) { test_batched_rotmg<TestDevice, float>(); }
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_rotmg_double) { test_batched_rotmg<TestDevice, double>(); }
#endif
