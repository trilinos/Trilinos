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
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Tbsv.hpp"
#include "Test_Batched_DenseUtils.hpp"

using namespace KokkosBatched;

namespace Test {
namespace Tbsv {

template <typename U, typename T, typename D>
struct ParamTag {
  using uplo  = U;
  using trans = T;
  using diag  = D;
};

template <typename DeviceType, typename AViewType, typename BViewType, typename ScalarType, typename ParamTagType,
          typename AlgoTagType>
struct Functor_BatchedSerialTrsv {
  using execution_space = typename DeviceType::execution_space;
  AViewType _a;
  BViewType _b;

  ScalarType _alpha;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialTrsv(const ScalarType alpha, const AViewType &a, const BViewType &b)
      : _a(a), _b(b), _alpha(alpha) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const ParamTagType &, const int k) const {
    auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb = Kokkos::subview(_b, k, Kokkos::ALL());

    KokkosBatched::SerialTrsv<typename ParamTagType::uplo, typename ParamTagType::trans, typename ParamTagType::diag,
                              AlgoTagType>::invoke(_alpha, aa, bb);
  }

  inline void run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialTbsv");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space, ParamTagType> policy(0, _b.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

template <typename DeviceType, typename AViewType, typename BViewType, typename ParamTagType, typename AlgoTagType>
struct Functor_BatchedSerialTbsv {
  using execution_space = typename DeviceType::execution_space;
  AViewType _a;
  BViewType _b;
  int _k;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialTbsv(const AViewType &a, const BViewType &b, const int k) : _a(a), _b(b), _k(k) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const ParamTagType &, const int k) const {
    auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb = Kokkos::subview(_b, k, Kokkos::ALL());

    KokkosBatched::SerialTbsv<typename ParamTagType::uplo, typename ParamTagType::trans, typename ParamTagType::diag,
                              AlgoTagType>::invoke(aa, bb, _k);
  }

  inline void run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialTbsv");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space, ParamTagType> policy(0, _b.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType, typename AlgoTagType>
/// \brief Implementation details of batched tbsv test
///
/// \param N [in] Batch size of RHS (banded matrix can also be batched matrix)
/// \param k [in] Number of superdiagonals or subdiagonals of matrix A
/// \param BlkSize [in] Block size of matrix A
void impl_test_batched_tbsv(const int N, const int k, const int BlkSize) {
  using execution_space = typename DeviceType::execution_space;
  using View2DType      = Kokkos::View<ScalarType **, LayoutType, execution_space>;
  using View3DType      = Kokkos::View<ScalarType ***, LayoutType, execution_space>;

  // Reference is created by trsv from triangular matrix
  View3DType A("A", N, BlkSize, BlkSize), Ref("Ref", N, BlkSize, BlkSize);
  View3DType Ab("Ab", N, k + 1, BlkSize);                 // Banded matrix
  View2DType x0("x0", N, BlkSize), x1("x1", N, BlkSize);  // Solutions

  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;
  Test::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(Ref, rand_pool, randStart, randEnd);
  Kokkos::fill_random(x0, rand_pool, randStart, randEnd);

  Kokkos::deep_copy(x1, x0);

  // Create triangluar or banded matrix
  create_banded_triangular_matrix<View3DType, View3DType, typename ParamTagType::uplo>(Ref, A, k, false);
  create_banded_triangular_matrix<View3DType, View3DType, typename ParamTagType::uplo>(Ref, Ab, k, true);

  // Reference trsv
  Functor_BatchedSerialTrsv<DeviceType, View3DType, View2DType, ScalarType, ParamTagType, Algo::Trsv::Unblocked>(1.0, A,
                                                                                                                 x0)
      .run();

  // tbsv
  Functor_BatchedSerialTbsv<DeviceType, View3DType, View2DType, ParamTagType, AlgoTagType>(Ab, x1, k).run();

  Kokkos::fence();

  // this eps is about 10^-14
  using ats      = typename Kokkos::ArithTraits<ScalarType>;
  using mag_type = typename ats::mag_type;
  mag_type eps   = 1.0e3 * ats::epsilon();

  // Check x0 = x1
  auto h_x0 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x0);
  auto h_x1 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x1);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < BlkSize; j++) {
      EXPECT_NEAR_KK(h_x0(i, j), h_x1(i, j), eps);
    }
  }
}

template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType, typename AlgoTagType>
/// \brief Implementation details of batched tbsv test
///
/// \param N [in] Batch size of RHS (banded matrix can also be batched matrix)
void impl_test_batched_tbsv_analytical(const std::size_t N) {
  using execution_space   = typename DeviceType::execution_space;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, execution_space>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, execution_space>;
  using View3DType        = Kokkos::View<ScalarType ***, LayoutType, execution_space>;

  // Reference is created by trsv from triangular matrix
  constexpr std::size_t BlkSize = 3, k = 2, incx = 2;

  View3DType A("A", N, BlkSize, BlkSize), ref("Ref", N, BlkSize, BlkSize);
  View3DType Ab("Ab", N, k + 1, BlkSize);                       // Banded matrix
  View2DType x0("x0", N, BlkSize), x_ref("x_ref", N, BlkSize);  // Solutions

  // Testing incx argument with strided Views
  Kokkos::LayoutStride layout{N, incx, BlkSize, N * incx};
  StridedView2DType x1("x1", layout);  // Solutions

  Kokkos::RangePolicy<execution_space> policy(0, N);
  Kokkos::parallel_for(
      "KokkosBatched::Test::SerialTbsv::Initialize", policy, KOKKOS_LAMBDA(const std::size_t ib) {
        for (std::size_t i = 0; i < BlkSize; i++) {
          for (std::size_t j = 0; j < BlkSize; j++) {
            ref(ib, i, j) = i + 1;
          }
        }
        for (std::size_t j = 0; j < BlkSize; j++) {
          x0(ib, j) = 1;
          x1(ib, j) = 1;
        }

        if (std::is_same_v<typename ParamTagType::uplo, KokkosBatched::Uplo::Upper>) {
          if (std::is_same_v<typename ParamTagType::trans, Trans::NoTranspose>) {
            if (std::is_same_v<typename ParamTagType::diag, Diag::NonUnit>) {
              x_ref(ib, 0) = 1.0 / 2.0;
              x_ref(ib, 1) = 1.0 / 6.0;
              x_ref(ib, 2) = 1.0 / 3.0;
            } else {
              x_ref(ib, 0) = 1.0;
              x_ref(ib, 1) = -1.0;
              x_ref(ib, 2) = 1.0;
            }
          } else {
            if (std::is_same_v<typename ParamTagType::diag, Diag::NonUnit>) {
              x_ref(ib, 0) = 1.0;
              x_ref(ib, 1) = 0.0;
              x_ref(ib, 2) = 0.0;
            } else {
              x_ref(ib, 0) = 1.0;
              x_ref(ib, 1) = 0.0;
              x_ref(ib, 2) = 0.0;
            }
          }
        } else {
          if (std::is_same_v<typename ParamTagType::trans, Trans::NoTranspose>) {
            if (std::is_same_v<typename ParamTagType::diag, Diag::NonUnit>) {
              x_ref(ib, 0) = 1.0;
              x_ref(ib, 1) = -1.0 / 2.0;
              x_ref(ib, 2) = -1.0 / 6.0;
            } else {
              x_ref(ib, 0) = 1.0;
              x_ref(ib, 1) = -1.0;
              x_ref(ib, 2) = 1.0;
            }
          } else {
            if (std::is_same_v<typename ParamTagType::diag, Diag::NonUnit>) {
              x_ref(ib, 0) = 0.0;
              x_ref(ib, 1) = 0.0;
              x_ref(ib, 2) = 1.0 / 3.0;
            } else {
              x_ref(ib, 0) = 2.0;
              x_ref(ib, 1) = -2.0;
              x_ref(ib, 2) = 1.0;
            }
          }
        }
      });

  Kokkos::fence();

  // Create triangluar or banded matrix
  create_banded_triangular_matrix<View3DType, View3DType, typename ParamTagType::uplo>(ref, A, k, false);
  create_banded_triangular_matrix<View3DType, View3DType, typename ParamTagType::uplo>(ref, Ab, k, true);

  // tbsv
  Functor_BatchedSerialTbsv<DeviceType, View3DType, View2DType, ParamTagType, AlgoTagType>(Ab, x0, k).run();

  // tbsv with incx == 2
  Functor_BatchedSerialTbsv<DeviceType, View3DType, StridedView2DType, ParamTagType, AlgoTagType>(Ab, x1, k).run();

  Kokkos::fence();

  // Check x0 = x_ref and x1 = x_ref
  // Firstly, prepare contiguous views on host
  auto h_x0 = Kokkos::create_mirror_view(x0);
  auto h_x1 = Kokkos::create_mirror_view(x0);

  Kokkos::deep_copy(h_x0, x0);

  // Pack x1 into x0 for contiguous storage
  Kokkos::parallel_for(
      "KokkosBatched::Test::SerialTbsv::Copy", policy, KOKKOS_LAMBDA(const std::size_t ib) {
        for (std::size_t j = 0; j < BlkSize; j++) {
          x0(ib, j) = x1(ib, j);
        }
      });

  Kokkos::fence();
  Kokkos::deep_copy(h_x1, x0);

  // this eps is about 10^-14
  using ats      = typename Kokkos::ArithTraits<ScalarType>;
  using mag_type = typename ats::mag_type;
  mag_type eps   = 1.0e3 * ats::epsilon();

  auto h_x_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x_ref);
  for (std::size_t ib = 0; ib < N; ib++) {
    for (std::size_t j = 0; j < BlkSize; j++) {
      // Check x0 = x_ref
      EXPECT_NEAR_KK(h_x0(ib, j), h_x_ref(ib, j), eps);

      // Check x1 = x_ref
      EXPECT_NEAR_KK(h_x1(ib, j), h_x_ref(ib, j), eps);
    }
  }
}

}  // namespace Tbsv
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
int test_batched_tbsv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Tbsv::impl_test_batched_tbsv_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(0);
    Test::Tbsv::impl_test_batched_tbsv_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1);
    Test::Tbsv::impl_test_batched_tbsv<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(0, 1, 10);
    for (int i = 0; i < 10; i++) {
      Test::Tbsv::impl_test_batched_tbsv<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1, 1, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Tbsv::impl_test_batched_tbsv_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(0);
    Test::Tbsv::impl_test_batched_tbsv_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1);
    Test::Tbsv::impl_test_batched_tbsv<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(0, 1, 10);
    for (int i = 0; i < 10; i++) {
      Test::Tbsv::impl_test_batched_tbsv<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1, 1, i);
    }
  }
#endif

  return 0;
}
