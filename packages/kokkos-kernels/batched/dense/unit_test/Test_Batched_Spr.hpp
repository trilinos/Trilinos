// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include <concepts>
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Spr.hpp>
#include <KokkosBatched_Syr.hpp>
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Spr {

template <typename U, typename T, typename M>
struct ParamTag {
  using uplo  = U;
  using trans = T;
  using mode  = M;
};

template <typename DeviceType, typename ScalarType, typename XViewType, typename APViewType, typename ParamTagType>
struct Functor_BatchedSpr {
  using execution_space = typename DeviceType::execution_space;
  using member_type     = typename Kokkos::TeamPolicy<execution_space>::member_type;
  using ArgMode         = typename ParamTagType::mode;
  using ArgUplo         = typename ParamTagType::uplo;
  using ArgTrans        = typename ParamTagType::trans;
  XViewType m_x;
  APViewType m_ap;
  ScalarType m_alpha;

  Functor_BatchedSpr(const ScalarType alpha, const XViewType &x, const APViewType &ap)
      : m_x(x), m_ap(ap), m_alpha(alpha) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type &member, int &info) const {
    const int k = member.league_rank();

    auto sub_x  = Kokkos::subview(m_x, k, Kokkos::ALL());
    auto sub_ap = Kokkos::subview(m_ap, k, Kokkos::ALL());

    if constexpr (std::same_as<ArgMode, KokkosBatched::Mode::Serial>) {
      Kokkos::single(Kokkos::PerTeam(member),
                     [&]() { info += KokkosBatched::SerialSpr<ArgUplo, ArgTrans>::invoke(m_alpha, sub_x, sub_ap); });
    } else if constexpr (std::same_as<ArgMode, KokkosBatched::Mode::Team>) {
      info += KokkosBatched::TeamSpr<member_type, ArgUplo, ArgTrans>::invoke(member, m_alpha, sub_x, sub_ap);
    } else if constexpr (std::same_as<ArgMode, KokkosBatched::Mode::TeamVector>) {
      info += KokkosBatched::TeamVectorSpr<member_type, ArgUplo, ArgTrans>::invoke(member, m_alpha, sub_x, sub_ap);
    }
  }

  inline int run() {
    using value_type        = typename APViewType::non_const_value_type;
    std::string name_region = std::same_as<ArgMode, KokkosBatched::Mode::Serial> ? "KokkosBatched::Test::SerialSpr"
                              : std::same_as<ArgMode, KokkosBatched::Mode::Team> ? "KokkosBatched::Test::TeamSpr"
                                                                                 : "KokkosBatched::Test::TeamVectorSpr";
    const std::string name_value_type = TestUtils::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    int info_sum                      = 0;
    Kokkos::Profiling::pushRegion(name.c_str());
    const int league_size = m_ap.extent_int(0);
    Kokkos::TeamPolicy<execution_space> policy(league_size, Kokkos::AUTO);
    Kokkos::parallel_reduce(name.c_str(), policy, *this, info_sum);
    Kokkos::Profiling::popRegion();
    return info_sum;
  }
};

template <typename DeviceType, typename ScalarType, typename XViewType, typename AViewType, typename ParamTagType>
struct Functor_BatchedSyr {
  using execution_space = typename DeviceType::execution_space;
  XViewType m_x;
  AViewType m_A;
  ScalarType m_alpha;

  Functor_BatchedSyr(const ScalarType alpha, const XViewType &x, const AViewType &A) : m_x(x), m_A(A), m_alpha(alpha) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    auto sub_x = Kokkos::subview(m_x, k, Kokkos::ALL());
    auto sub_A = Kokkos::subview(m_A, k, Kokkos::ALL(), Kokkos::ALL());

    KokkosBatched::SerialSyr<typename ParamTagType::uplo, typename ParamTagType::trans>::invoke(m_alpha, sub_x, sub_A);
  }

  inline void run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialSyr");
    const std::string name_value_type = TestUtils::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space> policy(0, m_A.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

/// \brief Implementation details of batched spr analytical test
///        to confirm A:= x*x**T + A is computed correctly
///        alpha = 1.5
///        4x4 matrix (upper)
///        U: [[1, -3, -2,  0],
///            [0,  3, -1, -2],
///            [0,  0,  9,  5],
///            [0,  0,  0, 27]]
///        U_packed: [1, -3, 3, -2, -1, 9, 0, -2, 5, 27]
///        x: [1, 2, 3, 4]
///        Ref: [2.5,  0.,  9.,  2.5,  8.,  22.5,  6.,  10.,  23.,  51.]
///
///        4x4 matrix (lower)
///        L: [[ 1  0  0  0]
///            [ 3  3  0  0]
///            [ 2  1  9  0]
///            [ 0  2 -5 27]]
///        L_packed:  [1, 3, 2, 0, 3, 1, 2, 9, -5, 27]
///        x: [1, 2, 3, 4]
///        Ref: [2.5, 6., 6.5, 6.,  9.,  10.,  14.,  22.5, 13.,  51.]
///
/// \param Nb [in] Batch size of matrices
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType>
void impl_test_batched_spr_analytical(const std::size_t Nb) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;
  using ArgUplo           = typename ParamTagType::uplo;

  const std::size_t BlkSize     = 4;
  const std::size_t packed_size = BlkSize * (BlkSize + 1) / 2;
  View2DType ap("ap", Nb, packed_size), ap_s("ap_s", Nb, packed_size), ap_ref("ap_ref", Nb, packed_size),
      x("x", Nb, BlkSize);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout{Nb, incx, BlkSize, Nb * incx};
  StridedView2DType x_s("x_s", layout);

  // Only filling x2, A2 and deep_copy from its subview
  auto h_ap     = Kokkos::create_mirror_view(ap);
  auto h_ap_ref = Kokkos::create_mirror_view(ap_ref);
  auto h_x      = Kokkos::create_mirror_view(x);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    if (std::is_same_v<ArgUplo, KokkosBatched::Uplo::Upper>) {
      h_ap(ib, 0) = 1;
      h_ap(ib, 1) = -3;
      h_ap(ib, 2) = 3;
      h_ap(ib, 3) = -2;
      h_ap(ib, 4) = -1;
      h_ap(ib, 5) = 9;
      h_ap(ib, 6) = 0;
      h_ap(ib, 7) = -2;
      h_ap(ib, 8) = 5;
      h_ap(ib, 9) = 27;

      h_ap_ref(ib, 0) = 2.5;
      h_ap_ref(ib, 1) = 0;
      h_ap_ref(ib, 2) = 9;
      h_ap_ref(ib, 3) = 2.5;
      h_ap_ref(ib, 4) = 8;
      h_ap_ref(ib, 5) = 22.5;
      h_ap_ref(ib, 6) = 6;
      h_ap_ref(ib, 7) = 10;
      h_ap_ref(ib, 8) = 23;
      h_ap_ref(ib, 9) = 51;
    } else {
      h_ap(ib, 0) = 1;
      h_ap(ib, 1) = 3;
      h_ap(ib, 2) = 2;
      h_ap(ib, 3) = 0;
      h_ap(ib, 4) = 3;
      h_ap(ib, 5) = 1;
      h_ap(ib, 6) = 2;
      h_ap(ib, 7) = 9;
      h_ap(ib, 8) = -5;
      h_ap(ib, 9) = 27;

      h_ap_ref(ib, 0) = 2.5;
      h_ap_ref(ib, 1) = 6;
      h_ap_ref(ib, 2) = 6.5;
      h_ap_ref(ib, 3) = 6;
      h_ap_ref(ib, 4) = 9;
      h_ap_ref(ib, 5) = 10;
      h_ap_ref(ib, 6) = 14;
      h_ap_ref(ib, 7) = 22.5;
      h_ap_ref(ib, 8) = 13;
      h_ap_ref(ib, 9) = 51;
    }

    for (std::size_t j = 0; j < BlkSize; j++) {
      h_x(ib, j) = static_cast<ScalarType>(j + 1);
    }
  }

  Kokkos::deep_copy(ap, h_ap);
  Kokkos::deep_copy(ap_s, ap);
  Kokkos::deep_copy(x, h_x);
  Kokkos::deep_copy(x_s, x);

  const ScalarType alpha = 1.5;

  auto info = Functor_BatchedSpr<DeviceType, ScalarType, View2DType, View2DType, ParamTagType>(alpha, x, ap).run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  // With strided views
  info =
      Functor_BatchedSpr<DeviceType, ScalarType, StridedView2DType, View2DType, ParamTagType>(alpha, x_s, ap_s).run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  RealType eps = 1.0e1 * ats::epsilon();
  Kokkos::deep_copy(h_ap, ap);
  auto h_ap_s = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, ap_s);

  // Check if A:= alpha * x * y**T + A
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < packed_size; i++) {
      EXPECT_NEAR_KK(h_ap(ib, i), h_ap_ref(ib, i), eps);
      EXPECT_NEAR_KK(h_ap_s(ib, i), h_ap_ref(ib, i), eps);
    }
  }
}

/// \brief Implementation details of batched spr test
///
/// \param[in] Nb Batch size of matrices
/// \param[in] BlkSize Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType>
void impl_test_batched_spr(const std::size_t Nb, const std::size_t BlkSize) {
  using ats        = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType   = typename ats::mag_type;
  using View2DType = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using View3DType = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using ArgUplo    = typename ParamTagType::uplo;
  using ArgTrans   = typename ParamTagType::trans;

  // Dense Matrix
  View3DType A("A", Nb, BlkSize, BlkSize), A0("A0", Nb, BlkSize, BlkSize), A_ref("A_ref", Nb, BlkSize, BlkSize),
      A0_ref("A0_ref", Nb, BlkSize, BlkSize);

  // Packed Matrix
  const std::size_t packed_size = BlkSize * (BlkSize + 1) / 2;
  View2DType ap("ap", Nb, packed_size), ap0("ap0", Nb, packed_size), ap_ref("ap_ref", Nb, packed_size),
      ap0_ref("ap0_ref", Nb, packed_size);

  // Vector x
  View2DType x("x", Nb, BlkSize);

  // Create a random matrix A and x
  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);
  Kokkos::fill_random(x, rand_pool, randStart, randEnd);

  // Upper or lower triangular part of A
  dense_to_symmetric<ArgUplo, ArgTrans>(A, A_ref);
  dense_to_packed<ArgUplo>(A_ref, ap);

  auto h_A_ref2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, A_ref);
  auto h_ap2    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, ap);

  // Cleanup A
  Kokkos::deep_copy(A, 0);

  // When A0 is zero (ap0 is already zero-initialized)
  const ScalarType alpha = 1.5;
  auto info0 = Functor_BatchedSpr<DeviceType, ScalarType, View2DType, View2DType, ParamTagType>(alpha, x, ap0).run();

  // When A is a random matrix
  auto info1 = Functor_BatchedSpr<DeviceType, ScalarType, View2DType, View2DType, ParamTagType>(alpha, x, ap).run();

  Kokkos::fence();
  EXPECT_EQ(info0, 0);
  EXPECT_EQ(info1, 0);

  // Make references with Syr
  Functor_BatchedSyr<DeviceType, ScalarType, View2DType, View3DType, ParamTagType>(alpha, x, A0_ref).run();
  Functor_BatchedSyr<DeviceType, ScalarType, View2DType, View3DType, ParamTagType>(alpha, x, A_ref).run();

  // Unpack
  packed_to_dense<ArgUplo, ArgTrans>(ap0, A0);
  packed_to_dense<ArgUplo, ArgTrans>(ap, A);

  RealType eps  = 1.0e1 * ats::epsilon();
  auto h_A      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, A);
  auto h_A0     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, A0);
  auto h_A_ref  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, A_ref);
  auto h_A0_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, A0_ref);

  // Check if A:= alpha * x * y**T + A or A:= alpha * x * y**H + A
  // Note: only the upper or lower part of the matrix is updated depending on ArgUplo by Syr
  for (std::size_t ib = 0; ib < Nb; ib++) {
    if (std::is_same_v<ArgUplo, KokkosBatched::Uplo::Upper>) {
      for (std::size_t i = 0; i < BlkSize; i++) {
        for (std::size_t j = 0; j < i + 1; j++) {
          EXPECT_NEAR_KK(h_A(ib, j, i), h_A_ref(ib, j, i), eps);
          EXPECT_NEAR_KK(h_A0(ib, j, i), h_A0_ref(ib, j, i), eps);
        }
      }
    } else {
      for (std::size_t i = 0; i < BlkSize; i++) {
        for (std::size_t j = i; j < BlkSize; j++) {
          EXPECT_NEAR_KK(h_A(ib, j, i), h_A_ref(ib, j, i), eps);
          EXPECT_NEAR_KK(h_A0(ib, j, i), h_A0_ref(ib, j, i), eps);
        }
      }
    }
  }
}

}  // namespace Spr
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ParamTagType>
int test_batched_spr() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Spr::impl_test_batched_spr_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(1);
    Test::Spr::impl_test_batched_spr_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Spr::impl_test_batched_spr<DeviceType, ScalarType, LayoutType, ParamTagType>(1, i);
      Test::Spr::impl_test_batched_spr<DeviceType, ScalarType, LayoutType, ParamTagType>(2, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Spr::impl_test_batched_spr_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(1);
    Test::Spr::impl_test_batched_spr_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Spr::impl_test_batched_spr<DeviceType, ScalarType, LayoutType, ParamTagType>(1, i);
      Test::Spr::impl_test_batched_spr<DeviceType, ScalarType, LayoutType, ParamTagType>(2, i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_serial_spr_l_t_float) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::Transpose, Mode::Serial>;
  test_batched_spr<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_spr_l_c_float) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::ConjTranspose, Mode::Serial>;
  test_batched_spr<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_spr_u_t_float) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::Transpose, Mode::Serial>;
  test_batched_spr<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_spr_u_c_float) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::ConjTranspose, Mode::Serial>;
  test_batched_spr<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_spr_l_t_float) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::Transpose, Mode::Team>;
  test_batched_spr<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_spr_l_c_float) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::ConjTranspose, Mode::Team>;
  test_batched_spr<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_spr_u_t_float) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::Transpose, Mode::Team>;
  test_batched_spr<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_spr_u_c_float) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::ConjTranspose, Mode::Team>;
  test_batched_spr<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_spr_l_t_float) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::Transpose, Mode::TeamVector>;
  test_batched_spr<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_spr_l_c_float) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::ConjTranspose, Mode::TeamVector>;
  test_batched_spr<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_spr_u_t_float) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::Transpose, Mode::TeamVector>;
  test_batched_spr<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_spr_u_c_float) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::ConjTranspose, Mode::TeamVector>;
  test_batched_spr<TestDevice, float, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_serial_spr_l_t_double) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::Transpose, Mode::Serial>;
  test_batched_spr<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_spr_l_c_double) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::ConjTranspose, Mode::Serial>;
  test_batched_spr<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_spr_u_t_double) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::Transpose, Mode::Serial>;
  test_batched_spr<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_spr_u_c_double) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::ConjTranspose, Mode::Serial>;
  test_batched_spr<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_spr_l_t_double) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::Transpose, Mode::Team>;
  test_batched_spr<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_spr_l_c_double) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::ConjTranspose, Mode::Team>;
  test_batched_spr<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_spr_u_t_double) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::Transpose, Mode::Team>;
  test_batched_spr<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_spr_u_c_double) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::ConjTranspose, Mode::Team>;
  test_batched_spr<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_spr_l_t_double) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::Transpose, Mode::TeamVector>;
  test_batched_spr<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_spr_l_c_double) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::ConjTranspose, Mode::TeamVector>;
  test_batched_spr<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_spr_u_t_double) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::Transpose, Mode::TeamVector>;
  test_batched_spr<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_spr_u_c_double) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::ConjTranspose, Mode::TeamVector>;
  test_batched_spr<TestDevice, double, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, test_batched_serial_spr_l_t_fcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::Transpose, Mode::Serial>;
  test_batched_spr<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_spr_l_c_fcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::ConjTranspose, Mode::Serial>;
  test_batched_spr<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_spr_u_t_fcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::Transpose, Mode::Serial>;
  test_batched_spr<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_spr_u_c_fcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::ConjTranspose, Mode::Serial>;
  test_batched_spr<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_spr_l_t_fcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::Transpose, Mode::Team>;
  test_batched_spr<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_spr_l_c_fcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::ConjTranspose, Mode::Team>;
  test_batched_spr<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_spr_u_t_fcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::Transpose, Mode::Team>;
  test_batched_spr<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_spr_u_c_fcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::ConjTranspose, Mode::Team>;
  test_batched_spr<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_spr_l_t_fcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::Transpose, Mode::TeamVector>;
  test_batched_spr<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_spr_l_c_fcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::ConjTranspose, Mode::TeamVector>;
  test_batched_spr<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_spr_u_t_fcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::Transpose, Mode::TeamVector>;
  test_batched_spr<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_spr_u_c_fcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::ConjTranspose, Mode::TeamVector>;
  test_batched_spr<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_serial_spr_l_t_dcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::Transpose, Mode::Serial>;
  test_batched_spr<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_spr_l_c_dcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::ConjTranspose, Mode::Serial>;
  test_batched_spr<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_spr_u_t_dcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::Transpose, Mode::Serial>;
  test_batched_spr<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_spr_u_c_dcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::ConjTranspose, Mode::Serial>;
  test_batched_spr<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_spr_l_t_dcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::Transpose, Mode::Team>;
  test_batched_spr<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_spr_l_c_dcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::ConjTranspose, Mode::Team>;
  test_batched_spr<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_spr_u_t_dcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::Transpose, Mode::Team>;
  test_batched_spr<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_spr_u_c_dcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::ConjTranspose, Mode::Team>;
  test_batched_spr<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_spr_l_t_dcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::Transpose, Mode::TeamVector>;
  test_batched_spr<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_spr_l_c_dcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Lower, Trans::ConjTranspose, Mode::TeamVector>;
  test_batched_spr<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_spr_u_t_dcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::Transpose, Mode::TeamVector>;
  test_batched_spr<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_spr_u_c_dcomplex) {
  using param_tag_type = ::Test::Spr::ParamTag<Uplo::Upper, Trans::ConjTranspose, Mode::TeamVector>;
  test_batched_spr<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
#endif
