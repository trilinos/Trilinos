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
#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "KokkosBatched_Trtri_Decl.hpp"
#include "KokkosBatched_Trtri_Serial_Impl.hpp"

#include "KokkosKernels_TestUtils.hpp"

#include <chrono>

#define PRINT_MAT 0

using namespace KokkosBatched;

namespace Test {
namespace Trtri {

template <class ViewTypeA, class ExecutionSpace>
struct UnitDiagTRTRI {
  ViewTypeA A_;
  using ScalarA = typename ViewTypeA::value_type;

  UnitDiagTRTRI(const ViewTypeA& A) : A_(A) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i) const { A_(i, i) = ScalarA(1); }
};
template <class ViewTypeA, class ExecutionSpace>
struct NonUnitDiagTRTRI {
  ViewTypeA A_;
  using ScalarA = typename ViewTypeA::value_type;

  NonUnitDiagTRTRI(const ViewTypeA& A) : A_(A) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i) const { A_(i, i) = A_(i, i) + 10; }
};
template <class ViewTypeA, class ViewTypeB, class ViewTypeC, class ExecutionSpace>
struct VanillaGEMM {
  bool A_t, B_t, A_c, B_c;
  int N, K;
  ViewTypeA A;
  ViewTypeB B;
  ViewTypeC C;

  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;
  typedef typename ViewTypeC::value_type ScalarC;
  typedef Kokkos::ArithTraits<ScalarC> APT;
  typedef typename APT::mag_type mag_type;
  ScalarA alpha;
  ScalarC beta;

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<ExecutionSpace>::member_type& team) const {
// GNU COMPILER BUG WORKAROUND
#if defined(KOKKOS_COMPILER_GNU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
    int i = team.league_rank();
#else
    const int i = team.league_rank();
#endif
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int& j) {
      ScalarC C_ij = 0.0;

      // GNU 5.3, 5.4 and 6.1 (and maybe more) crash with another nested lambda
      // here

#if defined(KOKKOS_COMPILER_GNU) && !defined(KOKKOS_COMPILER_NVCC)
      for (int k = 0; k < K; k++) {
        ScalarA A_ik = A_t ? (A_c ? APT::conj(A(k, i)) : A(k, i)) : A(i, k);
        ScalarB B_kj = B_t ? (B_c ? APT::conj(B(j, k)) : B(j, k)) : B(k, j);
        C_ij += A_ik * B_kj;
      }
#else
        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,K), [&] (const int& k, ScalarC& lsum) {
           ScalarA A_ik = A_t?(A_c?APT::conj(A(k,i)):A(k,i)):A(i,k);
           ScalarB B_kj = B_t?(B_c?APT::conj(B(j,k)):B(j,k)):B(k,j);
           lsum += A_ik*B_kj;
        },C_ij);
#endif

      C(i, j) = beta * C(i, j) + alpha * C_ij;
    });
  }
};

template <typename U, typename D>
struct ParamTag {
  typedef U uplo;
  typedef D diag;
};

template <typename DeviceType, typename ViewType, typename ParamTagType, typename AlgoTagType>
struct Functor_TestBatchedSerialTrtri {
  using execution_space = typename DeviceType::execution_space;
  ViewType _a;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedSerialTrtri(const ViewType& a) : _a(a) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const ParamTagType&, const int k) const {
    auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());

    SerialTrtri<typename ParamTagType::uplo, typename ParamTagType::diag, AlgoTagType>::invoke(aa);
  }

  inline void run() {
    typedef typename ViewType::value_type value_type;
    std::string name_region("KokkosBatched::Test::SerialTrtri");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space, ParamTagType> policy(0, _a.extent(0));
    Kokkos::parallel_for("Functor_TestBatchedSerialTrtri", policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename ViewType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
void impl_test_batched_trtri(const int N, const int K) {
  typedef typename ViewType::value_type value_type;
  typedef typename DeviceType::execution_space execution_space;
  typedef Kokkos::ArithTraits<value_type> ats;

  ScalarType alpha(1.0);
  ScalarType beta(0.0);

  // eps is ~ 10^-13 for double
  typedef typename ats::mag_type mag_type;
  const mag_type eps = 1.0e8 * ats::epsilon();
  bool fail_flag     = false;
  ScalarType cur_check_val;  // Either 1 or 0, to check A_I

  const bool is_A_lower = std::is_same<typename ParamTagType::uplo, Uplo::Lower>::value;
  ViewType A("A", N, K, K);
  ViewType A_original("A_original", N, K, K);
  ViewType A_I("A_I", N, K, K);

  typename ViewType::HostMirror I_host = Kokkos::create_mirror_view(A_I);
  typename ViewType::HostMirror A_host = Kokkos::create_mirror_view(A);

  uint64_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();

  using ViewTypeSubA = decltype(Kokkos::subview(A, 0, Kokkos::ALL(), Kokkos::ALL()));

  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);

  if (std::is_same<typename ParamTagType::diag, Diag::Unit>::value) {
    // Initialize A with deterministic random numbers
    Kokkos::fill_random(A, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarType>::max());
    using functor_type = UnitDiagTRTRI<ViewTypeSubA, execution_space>;
    for (int k = 0; k < N; ++k) {
      functor_type udtrtri(Kokkos::subview(A, k, Kokkos::ALL(), Kokkos::ALL()));
      // Initialize As diag with 1s
      Kokkos::parallel_for("KokkosBlas::Test::UnitDiagTRTRI", Kokkos::RangePolicy<execution_space>(0, K), udtrtri);
    }
  } else {  //(diag[0]=='N')||(diag[0]=='n')
    // Initialize A with random numbers
    Kokkos::fill_random(A, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarType>::max());
    using functor_type = NonUnitDiagTRTRI<ViewTypeSubA, execution_space>;
    for (int k = 0; k < N; ++k) {
      functor_type nudtrtri(Kokkos::subview(A, k, Kokkos::ALL(), Kokkos::ALL()));
      // Initialize As diag with A(i,i)+10
      Kokkos::parallel_for("KokkosBlas::Test::NonUnitDiagTRTRI", Kokkos::RangePolicy<execution_space>(0, K), nudtrtri);
    }
  }
  Kokkos::fence();

  Kokkos::deep_copy(A_host, A);
  // Make A_host a lower triangle
  for (int k = 0; k < N; k++) {
    if (is_A_lower) {
      for (int i = 0; i < K - 1; i++)
        for (int j = i + 1; j < K; j++) A_host(k, i, j) = ScalarType(0);
    } else {
      // Make A_host a upper triangle
      for (int i = 1; i < K; i++)
        for (int j = 0; j < i; j++) A_host(k, i, j) = ScalarType(0);
    }
  }
  Kokkos::deep_copy(A, A_host);
  Kokkos::deep_copy(A_original, A);
  Kokkos::fence();

#if PRINT_MAT
  printf("A_original:\n");
  for (int k = 0; k < N; ++k) {
    for (int i = 0; i < K; i++) {
      for (int j = 0; j < K; j++) {
        printf("%*.13lf ", 20, A_original(k, i, j));
      }
      printf("\n");
    }
  }
#endif

#if PRINT_MAT
  printf("A:\n");
  for (int k = 0; k < N; ++k) {
    for (int i = 0; i < K; i++) {
      for (int j = 0; j < K; j++) {
        printf("%*.13lf ", 20, A(k, i, j));
      }
      printf("\n");
    }
  }
#endif

  Functor_TestBatchedSerialTrtri<DeviceType, ViewType, ParamTagType, Algo::Trtri::Unblocked>(A).run();

#if PRINT_MAT
  printf("A_original:\n");
  for (int k = 0; k < N; ++k) {
    for (int i = 0; i < K; i++) {
      for (int j = 0; j < K; j++) {
        printf("%*.13lf ", 20, A_original(k, i, j));
      }
      printf("\n");
    }
  }
#endif

#if PRINT_MAT
  printf("A:\n");
  for (int k = 0; k < N; ++k) {
    for (int i = 0; i < K; i++) {
      for (int j = 0; j < K; j++) {
        printf("%*.13lf ", 20, A(k, i, j));
      }
      printf("\n");
    }
  }
#endif

  Kokkos::fence();

  struct VanillaGEMM<ViewTypeSubA, ViewTypeSubA, ViewTypeSubA, execution_space> vgemm;
  vgemm.A_t   = false;
  vgemm.B_t   = false;
  vgemm.A_c   = false;
  vgemm.B_c   = false;
  vgemm.N     = K;
  vgemm.K     = K;
  vgemm.alpha = alpha;
  vgemm.beta  = beta;
  for (int i = 0; i < N; i++) {
    vgemm.A = Kokkos::subview(A, i, Kokkos::ALL(), Kokkos::ALL());
    vgemm.B = Kokkos::subview(A_original, i, Kokkos::ALL(), Kokkos::ALL());
    ;
    vgemm.C = Kokkos::subview(A_I, i, Kokkos::ALL(), Kokkos::ALL());
    ;
    Kokkos::parallel_for("KokkosBlas::Test::VanillaGEMM", Kokkos::TeamPolicy<execution_space>(K, Kokkos::AUTO, 16),
                         vgemm);
  }

  Kokkos::fence();
  Kokkos::deep_copy(I_host, A_I);
  Kokkos::fence();

#if PRINT_MAT
  printf("I_host:\n");
  for (int k = 0; k < N; ++k) {
    for (int i = 0; i < K; i++) {
      for (int j = 0; j < K; j++) {
        printf("%*.13lf ", 20, I_host(k, i, j));
      }
      printf("\n");
    }
  }
#endif

  for (int k = 0; k < N; ++k) {
    for (int i = 0; i < K; ++i) {
      for (int j = 0; j < K; ++j) {
        cur_check_val = (i == j) ? ScalarType(1) : ScalarType(0);  // ats::abs(host_A(i,j));
        if (ats::abs(ats::abs(I_host(k, i, j)) - cur_check_val) > eps) {
          fail_flag = true;
          // printf("   Error: eps ( %g ), I_host ( %.15f ) != cur_check_val
          // (%.15f) (abs result-cur_check_val %g) at (k %d, i %d, j %d)\n",
          // eps, I_host(k,i,j), cur_check_val, ats::abs(I_host(k,i,j) -
          // cur_check_val), k, i, j);
        }
      }
    }
  }

  ASSERT_EQ(fail_flag, false);
}
}  // namespace Trtri
}  // namespace Test

template <typename DeviceType, typename ValueType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
int test_batched_trtri(int batchSize = 512) {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    typedef Kokkos::View<ValueType***, Kokkos::LayoutLeft, DeviceType> ViewType;
    Test::Trtri::impl_test_batched_trtri<DeviceType, ViewType, ScalarType, ParamTagType, AlgoTagType>(0, 10);
    // Test::impl_test_batched_trtri<DeviceType,ViewType,ScalarType,ParamTagType,AlgoTagType>(
    // 1, 2);
    for (int i = 0; i < 10; ++i) {
      // printf("Testing: LayoutLeft,  Blksize %d\n", i);
      Test::Trtri::impl_test_batched_trtri<DeviceType, ViewType, ScalarType, ParamTagType, AlgoTagType>(batchSize, i);
      Test::Trtri::impl_test_batched_trtri<DeviceType, ViewType, ScalarType, ParamTagType, AlgoTagType>(batchSize, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    typedef Kokkos::View<ValueType***, Kokkos::LayoutRight, DeviceType> ViewType;
    Test::Trtri::impl_test_batched_trtri<DeviceType, ViewType, ScalarType, ParamTagType, AlgoTagType>(0, 10);
    for (int i = 0; i < 10; ++i) {
      // printf("Testing: LayoutRight, Blksize %d\n", i);
      Test::Trtri::impl_test_batched_trtri<DeviceType, ViewType, ScalarType, ParamTagType, AlgoTagType>(batchSize, i);
      Test::Trtri::impl_test_batched_trtri<DeviceType, ViewType, ScalarType, ParamTagType, AlgoTagType>(batchSize, i);
    }
  }
#endif

  return 0;
}
