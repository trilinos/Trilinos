// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "KokkosBatched_Copy_Decl.hpp"
#include "KokkosKernels_TestUtils.hpp"

namespace Test {
namespace Copy {

template <typename DeviceType, typename AViewType, typename BViewType, typename ArgTrans>
struct Functor_TestBatchedSerialCopy {
  using execution_space = typename DeviceType::execution_space;
  const AViewType m_A;
  const BViewType m_B;

  Functor_TestBatchedSerialCopy(const AViewType &A, const BViewType &B) : m_A(A), m_B(B) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k, int &info) const {
    if constexpr (AViewType::rank() == 3) {
      auto sub_A = Kokkos::subview(m_A, k, Kokkos::ALL, Kokkos::ALL);
      auto sub_B = Kokkos::subview(m_B, k, Kokkos::ALL, Kokkos::ALL);
      info += KokkosBatched::SerialCopy<ArgTrans>::invoke(sub_A, sub_B);
    } else {
      auto sub_A = Kokkos::subview(m_A, k, Kokkos::ALL);
      auto sub_B = Kokkos::subview(m_B, k, Kokkos::ALL);
      info += KokkosBatched::SerialCopy<ArgTrans>::invoke(sub_A, sub_B);
    }
  }

  inline int run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialCopy");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    int info_sum                      = 0;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space> policy(0, m_A.extent(0));
    Kokkos::parallel_reduce(name.c_str(), policy, *this, info_sum);
    Kokkos::Profiling::popRegion();
    return info_sum;
  }
};

template <typename DeviceType, typename AViewType, typename BViewType, typename ArgTrans, typename ArgMode>
struct Functor_TestBatchedTeamCopy {
  using execution_space = typename DeviceType::execution_space;
  const AViewType m_A;
  const BViewType m_B;

  Functor_TestBatchedTeamCopy(const AViewType &A, const BViewType &B) : m_A(A), m_B(B) {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    const int k = member.league_rank();
    if constexpr (AViewType::rank() == 3) {
      auto sub_A = Kokkos::subview(m_A, k, Kokkos::ALL, Kokkos::ALL);
      auto sub_B = Kokkos::subview(m_B, k, Kokkos::ALL, Kokkos::ALL);

      if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Team>) {
        KokkosBatched::TeamCopy<MemberType, ArgTrans>::invoke(member, sub_A, sub_B);
      } else if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::TeamVector>) {
        KokkosBatched::TeamVectorCopy<MemberType, ArgTrans>::invoke(member, sub_A, sub_B);
      }
    } else {
      auto sub_A = Kokkos::subview(m_A, k, Kokkos::ALL);
      auto sub_B = Kokkos::subview(m_B, k, Kokkos::ALL);
      if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Team>) {
        KokkosBatched::TeamCopy<MemberType, ArgTrans>::invoke(member, sub_A, sub_B);
      } else if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::TeamVector>) {
        KokkosBatched::TeamVectorCopy<MemberType, ArgTrans>::invoke(member, sub_A, sub_B);
      }
    }
  }

  inline void run() {
    using value_type                  = typename AViewType::non_const_value_type;
    std::string name_region           = std::is_same_v<ArgMode, KokkosBatched::Mode::Team>
                                            ? "KokkosBatched::Test::TeamCopy"
                                            : "KokkosBatched::Test::TeamVectorCopy";
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    const int league_size = m_A.extent_int(0);
    Kokkos::TeamPolicy<execution_space> policy(league_size, Kokkos::AUTO);
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

/// \brief Implementation details of batched copy analytical test
///        Confirm B == A or B == A^T or A^H
///
/// A0: [[1],[2],[3]]
/// B0: [[1],[2],[3]] if no transpose
/// B0: [[1,2,3]] if transpose or conjugate transpose
///
/// A1: [[1,2,3],[4,5,6]]
/// B1: [[1,2,3],[4,5,6]] if no transpose
/// B1: [[1,4],[2,5],[3,6]] if transpose or conjugate transpose
/// \tparam DeviceType Kokkos device type
/// \tparam ScalarType Value type of input/output views
/// \tparam LayoutType Layout type of input/output views
/// \tparam ArgTrans Type indicating whether the transpose (Trans::Transpose) or conjugate transpose
/// (Trans::ConjTranspose) of A is used
/// \tparam ArgMode: one of Mode::Serial, Mode::Team, Mode::TeamVector
/// \param[in] N Batch size of matrices
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ArgTrans, typename ArgMode>
void impl_test_batched_copy_analytical(const std::size_t N) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using View3DType        = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;
  using StridedView3DType = Kokkos::View<ScalarType ***, Kokkos::LayoutStride, DeviceType>;

  const std::size_t m = std::is_same_v<ArgTrans, KokkosBatched::Trans::NoTranspose> ? 2 : 3;
  const std::size_t n = std::is_same_v<ArgTrans, KokkosBatched::Trans::NoTranspose> ? 3 : 2;
  View2DType A0("A0", N, 3), B0("B0", N, 3), Ref_B0("Ref_B0", N, 3);
  View3DType A1("A1", N, 2, 3), B1("B1", N, m, n), Ref_B1("Ref_B1", N, m, n);

  // Testing with strided views
  const std::size_t stride = 2;
  Kokkos::LayoutStride layout0{N, stride, 3, N * stride};
  StridedView2DType A0_s("A0_s", layout0), B0_s("B0_s", layout0);

  Kokkos::LayoutStride layout1{N, stride, 2, N * stride, 3, N * stride * 2},
      layout2{N, stride, m, N * stride, n, N * stride * m};

  StridedView3DType A1_s("A1_s", layout1), B1_s("B1_s", layout2);

  // Initialize A0 and A1
  auto h_A0     = Kokkos::create_mirror_view(A0);
  auto h_A1     = Kokkos::create_mirror_view(A1);
  auto h_Ref_B1 = Kokkos::create_mirror_view(Ref_B1);

  for (std::size_t ib = 0; ib < N; ib++) {
    h_A0(ib, 0) = 1.0;
    h_A0(ib, 1) = 2.0;
    h_A0(ib, 2) = 3.0;

    h_A1(ib, 0, 0) = 1.0;
    h_A1(ib, 0, 1) = 2.0;
    h_A1(ib, 0, 2) = 3.0;
    h_A1(ib, 1, 0) = 4.0;
    h_A1(ib, 1, 1) = 5.0;
    h_A1(ib, 1, 2) = 6.0;

    for (std::size_t i = 0; i < h_Ref_B1.extent(1); i++) {
      for (std::size_t j = 0; j < h_Ref_B1.extent(2); j++) {
        h_Ref_B1(ib, i, j) =
            std::is_same_v<ArgTrans, KokkosBatched::Trans::NoTranspose> ? h_A1(ib, i, j) : h_A1(ib, j, i);
      }
    }
  }
  Kokkos::deep_copy(A0, h_A0);
  Kokkos::deep_copy(A1, h_A1);
  Kokkos::deep_copy(Ref_B0, A0);

  // Strided views can be copied only on the same device
  Kokkos::deep_copy(A0_s, A0);
  Kokkos::deep_copy(A1_s, A1);

  if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Serial>) {
    auto info0 = Functor_TestBatchedSerialCopy<DeviceType, View2DType, View2DType, ArgTrans>(A0, B0).run();
    auto info1 = Functor_TestBatchedSerialCopy<DeviceType, View3DType, View3DType, ArgTrans>(A1, B1).run();

    EXPECT_EQ(info0, 0);
    EXPECT_EQ(info1, 0);

    // For strided views
    info0 = Functor_TestBatchedSerialCopy<DeviceType, StridedView2DType, StridedView2DType, ArgTrans>(A0_s, B0_s).run();
    info1 = Functor_TestBatchedSerialCopy<DeviceType, StridedView3DType, StridedView3DType, ArgTrans>(A1_s, B1_s).run();

    EXPECT_EQ(info0, 0);
    EXPECT_EQ(info1, 0);
  } else {
    // Team or TeamVector
    Functor_TestBatchedTeamCopy<DeviceType, View2DType, View2DType, ArgTrans, ArgMode>(A0, B0).run();
    Functor_TestBatchedTeamCopy<DeviceType, View3DType, View3DType, ArgTrans, ArgMode>(A1, B1).run();
    Functor_TestBatchedTeamCopy<DeviceType, StridedView2DType, StridedView2DType, ArgTrans, ArgMode>(A0_s, B0_s).run();
    Functor_TestBatchedTeamCopy<DeviceType, StridedView3DType, StridedView3DType, ArgTrans, ArgMode>(A1_s, B1_s).run();
  }

  RealType eps  = 1.0e1 * ats::epsilon();
  auto h_B0     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, B0);
  auto h_B1     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, B1);
  auto h_Ref_B0 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, Ref_B0);

  // Check if B0 and B1 are same as Ref_B0 and Ref_B1
  for (std::size_t ib = 0; ib < h_B0.extent(0); ib++) {
    for (std::size_t j = 0; j < h_B0.extent(1); j++) {
      EXPECT_NEAR_KK(h_B0(ib, j), h_Ref_B0(ib, j), eps);
    }
    for (std::size_t i = 0; i < h_B1.extent(1); i++) {
      for (std::size_t j = 0; j < h_B1.extent(2); j++) {
        EXPECT_NEAR_KK(h_B1(ib, i, j), h_Ref_B1(ib, i, j), eps);
      }
    }
  }

  // Testing for strided views, reusing B0 and B1
  Kokkos::deep_copy(B0, B0_s);
  Kokkos::deep_copy(B1, B1_s);
  Kokkos::deep_copy(h_B0, B0);
  Kokkos::deep_copy(h_B1, B1);
  for (std::size_t ib = 0; ib < h_B0.extent(0); ib++) {
    for (std::size_t j = 0; j < h_B0.extent(1); j++) {
      EXPECT_NEAR_KK(h_B0(ib, j), h_Ref_B0(ib, j), eps);
    }
    for (std::size_t i = 0; i < h_B1.extent(1); i++) {
      for (std::size_t j = 0; j < h_B1.extent(2); j++) {
        EXPECT_NEAR_KK(h_B1(ib, i, j), h_Ref_B1(ib, i, j), eps);
      }
    }
  }
}

/// \brief Implementation details of batched copy test
///        Confirm B == A or B == A^T or A^H
///
/// \tparam DeviceType Kokkos device type
/// \tparam ScalarType Value type of input/output views
/// \tparam LayoutType Layout type of input/output views
/// \tparam ArgTrans Type indicating whether the transpose (Trans::Transpose) or conjugate transpose
/// (Trans::ConjTranspose) of A is used
/// \tparam ArgMode: one of Mode::Serial, Mode::Team, Mode::TeamVector
/// \param[in] N Batch size of matrices
/// \param[in] m_A Number of rows of matrix A
/// \param[in] n_A Number of columns of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ArgTrans, typename ArgMode>
void impl_test_batched_copy(const std::size_t N, const std::size_t m_A, const std::size_t n_A) {
  using ats        = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType   = typename ats::mag_type;
  using View2DType = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using View3DType = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;

  const std::size_t m_B = std::is_same_v<ArgTrans, KokkosBatched::Trans::NoTranspose> ? m_A : n_A;
  const std::size_t n_B = std::is_same_v<ArgTrans, KokkosBatched::Trans::NoTranspose> ? n_A : m_A;
  View2DType A0("A0", N, m_A), B0("B0", N, m_A), Ref_B0("Ref_B0", N, m_A);
  View3DType A1("A1", N, m_A, n_A), B1("B1", N, m_B, n_B), Ref_B1("Ref_B1", N, m_B, n_B);

  // Initialize A0 and A1 with random numbers
  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A0, rand_pool, randStart, randEnd);
  Kokkos::fill_random(A1, rand_pool, randStart, randEnd);

  using Op = std::conditional_t<std::is_same_v<ArgTrans, KokkosBatched::Trans::ConjTranspose>, KokkosBlas::Impl::OpConj,
                                KokkosBlas::Impl::OpID>;
  Op op;

  auto h_A0     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, A0);
  auto h_A1     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, A1);
  auto h_Ref_B0 = Kokkos::create_mirror_view(Ref_B0);
  auto h_Ref_B1 = Kokkos::create_mirror_view(Ref_B1);

  for (std::size_t ib = 0; ib < N; ib++) {
    for (std::size_t j = 0; j < h_Ref_B0.extent(1); j++) {
      h_Ref_B0(ib, j) = op(h_A0(ib, j));
    }
    for (std::size_t i = 0; i < h_Ref_B1.extent(1); i++) {
      for (std::size_t j = 0; j < h_Ref_B1.extent(2); j++) {
        h_Ref_B1(ib, i, j) =
            std::is_same_v<ArgTrans, KokkosBatched::Trans::NoTranspose> ? h_A1(ib, i, j) : op(h_A1(ib, j, i));
      }
    }
  }

  // Copy operations
  if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Serial>) {
    auto info0 = Functor_TestBatchedSerialCopy<DeviceType, View2DType, View2DType, ArgTrans>(A0, B0).run();
    auto info1 = Functor_TestBatchedSerialCopy<DeviceType, View3DType, View3DType, ArgTrans>(A1, B1).run();
    EXPECT_EQ(info0, 0);
    EXPECT_EQ(info1, 0);
  } else {
    Functor_TestBatchedTeamCopy<DeviceType, View2DType, View2DType, ArgTrans, ArgMode>(A0, B0).run();
    Functor_TestBatchedTeamCopy<DeviceType, View3DType, View3DType, ArgTrans, ArgMode>(A1, B1).run();
  }

  RealType eps = 1.0e1 * ats::epsilon();
  auto h_B0    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, B0);
  auto h_B1    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, B1);

  // Check if B0 and B1 are same as Ref_B0 and Ref_B1
  for (std::size_t ib = 0; ib < h_B0.extent(0); ib++) {
    for (std::size_t j = 0; j < h_B0.extent(1); j++) {
      EXPECT_NEAR_KK(h_B0(ib, j), h_Ref_B0(ib, j), eps);
    }
    for (std::size_t i = 0; i < h_B1.extent(1); i++) {
      for (std::size_t j = 0; j < h_B1.extent(2); j++) {
        EXPECT_NEAR_KK(h_B1(ib, i, j), h_Ref_B1(ib, i, j), eps);
      }
    }
  }
}

}  // namespace Copy
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ArgTrans, typename ArgMode>
int test_batched_copy() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Copy::impl_test_batched_copy_analytical<DeviceType, ScalarType, LayoutType, ArgTrans, ArgMode>(1);
    Test::Copy::impl_test_batched_copy_analytical<DeviceType, ScalarType, LayoutType, ArgTrans, ArgMode>(10);

    for (std::size_t ib = 0; ib < 3; ib++) {
      for (std::size_t m = 0; m <= 3; m++) {
        for (std::size_t n = 0; n <= 3; n++) {
          Test::Copy::impl_test_batched_copy<DeviceType, ScalarType, LayoutType, ArgTrans, ArgMode>(ib, m, n);
        }
      }
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Copy::impl_test_batched_copy_analytical<DeviceType, ScalarType, LayoutType, ArgTrans, ArgMode>(1);
    Test::Copy::impl_test_batched_copy_analytical<DeviceType, ScalarType, LayoutType, ArgTrans, ArgMode>(10);

    for (std::size_t ib = 0; ib < 3; ib++) {
      for (std::size_t m = 0; m <= 3; m++) {
        for (std::size_t n = 0; n <= 3; n++) {
          Test::Copy::impl_test_batched_copy<DeviceType, ScalarType, LayoutType, ArgTrans, ArgMode>(ib, m, n);
        }
      }
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
// Serial
TEST_F(TestCategory, test_batched_serial_copy_nt_float) {
  test_batched_copy<TestDevice, float, KokkosBatched::Trans::NoTranspose, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_copy_t_float) {
  test_batched_copy<TestDevice, float, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_copy_c_float) {
  test_batched_copy<TestDevice, float, KokkosBatched::Trans::ConjTranspose, KokkosBatched::Mode::Serial>();
}

// Team
TEST_F(TestCategory, test_batched_team_copy_nt_float) {
  test_batched_copy<TestDevice, float, KokkosBatched::Trans::NoTranspose, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_copy_t_float) {
  test_batched_copy<TestDevice, float, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_copy_c_float) {
  test_batched_copy<TestDevice, float, KokkosBatched::Trans::ConjTranspose, KokkosBatched::Mode::Team>();
}

// Team Vector
TEST_F(TestCategory, test_batched_teamvector_copy_nt_float) {
  test_batched_copy<TestDevice, float, KokkosBatched::Trans::NoTranspose, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_copy_t_float) {
  test_batched_copy<TestDevice, float, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_copy_c_float) {
  test_batched_copy<TestDevice, float, KokkosBatched::Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
// Serial
TEST_F(TestCategory, test_batched_serial_copy_nt_double) {
  test_batched_copy<TestDevice, double, KokkosBatched::Trans::NoTranspose, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_copy_t_double) {
  test_batched_copy<TestDevice, double, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_copy_c_double) {
  test_batched_copy<TestDevice, double, KokkosBatched::Trans::ConjTranspose, KokkosBatched::Mode::Serial>();
}

// Team
TEST_F(TestCategory, test_batched_team_copy_nt_double) {
  test_batched_copy<TestDevice, double, KokkosBatched::Trans::NoTranspose, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_copy_t_double) {
  test_batched_copy<TestDevice, double, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_copy_c_double) {
  test_batched_copy<TestDevice, double, KokkosBatched::Trans::ConjTranspose, KokkosBatched::Mode::Team>();
}

// Team Vector
TEST_F(TestCategory, test_batched_teamvector_copy_nt_double) {
  test_batched_copy<TestDevice, double, KokkosBatched::Trans::NoTranspose, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_copy_t_double) {
  test_batched_copy<TestDevice, double, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_copy_c_double) {
  test_batched_copy<TestDevice, double, KokkosBatched::Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
// Serial
TEST_F(TestCategory, test_batched_serial_copy_nt_fcomplex) {
  test_batched_copy<TestDevice, Kokkos::complex<float>, KokkosBatched::Trans::NoTranspose,
                    KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_copy_t_fcomplex) {
  test_batched_copy<TestDevice, Kokkos::complex<float>, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_copy_c_fcomplex) {
  test_batched_copy<TestDevice, Kokkos::complex<float>, KokkosBatched::Trans::ConjTranspose,
                    KokkosBatched::Mode::Serial>();
}

// Team
TEST_F(TestCategory, test_batched_team_copy_nt_fcomplex) {
  test_batched_copy<TestDevice, Kokkos::complex<float>, KokkosBatched::Trans::NoTranspose, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_copy_t_fcomplex) {
  test_batched_copy<TestDevice, Kokkos::complex<float>, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_copy_c_fcomplex) {
  test_batched_copy<TestDevice, Kokkos::complex<float>, KokkosBatched::Trans::ConjTranspose,
                    KokkosBatched::Mode::Team>();
}

// Team Vector
TEST_F(TestCategory, test_batched_teamvector_copy_nt_fcomplex) {
  test_batched_copy<TestDevice, Kokkos::complex<float>, KokkosBatched::Trans::NoTranspose,
                    KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_copy_t_fcomplex) {
  test_batched_copy<TestDevice, Kokkos::complex<float>, KokkosBatched::Trans::Transpose,
                    KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_copy_c_fcomplex) {
  test_batched_copy<TestDevice, Kokkos::complex<float>, KokkosBatched::Trans::ConjTranspose,
                    KokkosBatched::Mode::TeamVector>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
// Serial
TEST_F(TestCategory, test_batched_serial_copy_nt_dcomplex) {
  test_batched_copy<TestDevice, Kokkos::complex<double>, KokkosBatched::Trans::NoTranspose,
                    KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_copy_t_dcomplex) {
  test_batched_copy<TestDevice, Kokkos::complex<double>, KokkosBatched::Trans::Transpose,
                    KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_copy_c_dcomplex) {
  test_batched_copy<TestDevice, Kokkos::complex<double>, KokkosBatched::Trans::ConjTranspose,
                    KokkosBatched::Mode::Serial>();
}

// Team
TEST_F(TestCategory, test_batched_team_copy_nt_dcomplex) {
  test_batched_copy<TestDevice, Kokkos::complex<double>, KokkosBatched::Trans::NoTranspose,
                    KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_copy_t_dcomplex) {
  test_batched_copy<TestDevice, Kokkos::complex<double>, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_copy_c_dcomplex) {
  test_batched_copy<TestDevice, Kokkos::complex<double>, KokkosBatched::Trans::ConjTranspose,
                    KokkosBatched::Mode::Team>();
}

// Team Vector
TEST_F(TestCategory, test_batched_teamvector_copy_nt_dcomplex) {
  test_batched_copy<TestDevice, Kokkos::complex<double>, KokkosBatched::Trans::NoTranspose,
                    KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_copy_t_dcomplex) {
  test_batched_copy<TestDevice, Kokkos::complex<double>, KokkosBatched::Trans::Transpose,
                    KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_copy_c_dcomplex) {
  test_batched_copy<TestDevice, Kokkos::complex<double>, KokkosBatched::Trans::ConjTranspose,
                    KokkosBatched::Mode::TeamVector>();
}
#endif
