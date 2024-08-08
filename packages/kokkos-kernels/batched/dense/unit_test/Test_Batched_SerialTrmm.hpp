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

#include "KokkosBatched_Trmm_Decl.hpp"
#include "KokkosBatched_Trmm_Serial_Impl.hpp"

#include "KokkosKernels_TestUtils.hpp"

#include <chrono>

using namespace KokkosBatched;

namespace Test {
namespace Trmm {

template <class ViewTypeA, class ExecutionSpace>
struct UnitDiagTRMM {
  ViewTypeA A_;
  using ScalarA = typename ViewTypeA::value_type;

  UnitDiagTRMM(const ViewTypeA& A) : A_(A) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i) const { A_(i, i) = ScalarA(1); }
};
template <class ViewTypeA, class ExecutionSpace>
struct NonUnitDiagTRMM {
  ViewTypeA A_;
  using ScalarA = typename ViewTypeA::value_type;

  NonUnitDiagTRMM(const ViewTypeA& A) : A_(A) {}

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

template <typename S, typename U, typename T, typename D>
struct ParamTag {
  typedef S side;
  typedef U uplo;
  typedef T trans;
  typedef D diag;
};

template <typename DeviceType, typename ViewType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
struct Functor_TestBatchedSerialTrmm {
  using execution_space = typename DeviceType::execution_space;
  ViewType _a, _b;

  ScalarType _alpha;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedSerialTrmm(const ScalarType alpha, const ViewType& a, const ViewType& b)
      : _a(a), _b(b), _alpha(alpha) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const ParamTagType&, const int k) const {
    auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb = Kokkos::subview(_b, k, Kokkos::ALL(), Kokkos::ALL());

    SerialTrmm<typename ParamTagType::side, typename ParamTagType::uplo, typename ParamTagType::trans,
               typename ParamTagType::diag, AlgoTagType>::invoke(_alpha, aa, bb);
  }

  inline void run() {
    typedef typename ViewType::value_type value_type;
    std::string name_region("KokkosBatched::Test::SerialTrmm");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space, ParamTagType> policy(0, _a.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename ViewType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
void impl_test_batched_trmm(const int N, const int nRows, const int nCols, const char* trans) {
  typedef typename ViewType::value_type value_type;
  typedef typename DeviceType::execution_space execution_space;
  typedef Kokkos::ArithTraits<value_type> ats;

  ScalarType alpha(1.0);
  ScalarType beta(0.0);

  const bool is_side_right = std::is_same<typename ParamTagType::side, Side::Right>::value;
  const bool is_A_lower    = std::is_same<typename ParamTagType::uplo, Uplo::Lower>::value;
  const int K              = is_side_right ? nCols : nRows;
  ViewType A("A", N, K, K), B_actual("B_actual", N, nRows, nCols), B_expected("B_expected", N, nRows, nCols);
  typename ViewType::HostMirror A_host          = Kokkos::create_mirror_view(A);
  typename ViewType::HostMirror B_actual_host   = Kokkos::create_mirror_view(B_actual);
  typename ViewType::HostMirror B_expected_host = Kokkos::create_mirror_view(B_expected);
  uint64_t seed                                 = std::chrono::high_resolution_clock::now().time_since_epoch().count();

  using ViewTypeSubA = decltype(Kokkos::subview(A, 0, Kokkos::ALL(), Kokkos::ALL()));
  using ViewTypeSubB = decltype(Kokkos::subview(B_actual, 0, Kokkos::ALL(), Kokkos::ALL()));

  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);

  if (std::is_same<typename ParamTagType::diag, Diag::NonUnit>::value) {
    // Initialize A with deterministic random numbers
    Kokkos::fill_random(A, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarType>::max());
    using functor_type = UnitDiagTRMM<ViewTypeSubA, execution_space>;
    for (int k = 0; k < N; ++k) {
      functor_type udtrmm(Kokkos::subview(A, k, Kokkos::ALL(), Kokkos::ALL()));
      // Initialize As diag with 1s
      Kokkos::parallel_for("KokkosBlas::Test::UnitDiagTRMM", Kokkos::RangePolicy<execution_space>(0, K), udtrmm);
    }
  } else {  //(diag[0]=='N')||(diag[0]=='n')
    // Initialize A with random numbers
    Kokkos::fill_random(A, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarType>::max());
    using functor_type = NonUnitDiagTRMM<ViewTypeSubA, execution_space>;
    for (int k = 0; k < N; ++k) {
      functor_type nudtrmm(Kokkos::subview(A, k, Kokkos::ALL(), Kokkos::ALL()));
      // Initialize As diag with A(i,i)+10
      Kokkos::parallel_for("KokkosBlas::Test::NonUnitDiagTRMM", Kokkos::RangePolicy<execution_space>(0, K), nudtrmm);
    }
  }
  Kokkos::fill_random(B_actual, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarType>::max());
  Kokkos::fence();

  Kokkos::deep_copy(B_expected, B_actual);
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

  if (!is_side_right) {
    // B_expected = alpha * op(A) * B + beta * C = 1 * op(A) * B + 0 * C
    struct VanillaGEMM<ViewTypeSubA, ViewTypeSubB, ViewTypeSubB, execution_space> vgemm;
    vgemm.A_t   = (trans[0] != 'N') && (trans[0] != 'n');
    vgemm.B_t   = false;
    vgemm.A_c   = (trans[0] == 'C') || (trans[0] == 'c');
    vgemm.B_c   = false;
    vgemm.N     = nCols;
    vgemm.K     = K;
    vgemm.alpha = alpha;
    vgemm.beta  = beta;
    for (int i = 0; i < N; i++) {
      vgemm.A = Kokkos::subview(A, i, Kokkos::ALL(), Kokkos::ALL());
      vgemm.B = Kokkos::subview(B_actual, i, Kokkos::ALL(), Kokkos::ALL());
      ;
      vgemm.C = Kokkos::subview(B_expected, i, Kokkos::ALL(), Kokkos::ALL());
      ;
      Kokkos::parallel_for("KokkosBlas::Test::VanillaGEMM",
                           Kokkos::TeamPolicy<execution_space>(nRows, Kokkos::AUTO, 16), vgemm);
    }
  } else {
    // B_expected = alpha * B * op(A) + beta * C = 1 * B * op(A) + 0 * C
    struct VanillaGEMM<ViewTypeSubB, ViewTypeSubA, ViewTypeSubB, execution_space> vgemm;
    vgemm.A_t   = false;
    vgemm.B_t   = (trans[0] != 'N') && (trans[0] != 'n');
    vgemm.A_c   = false;
    vgemm.B_c   = (trans[0] == 'C') || (trans[0] == 'c');
    vgemm.N     = nCols;
    vgemm.K     = K;
    vgemm.alpha = alpha;
    vgemm.beta  = beta;
    for (int i = 0; i < N; i++) {
      vgemm.A = Kokkos::subview(B_actual, i, Kokkos::ALL(), Kokkos::ALL());
      vgemm.B = Kokkos::subview(A, i, Kokkos::ALL(), Kokkos::ALL());
      ;
      vgemm.C = Kokkos::subview(B_expected, i, Kokkos::ALL(), Kokkos::ALL());
      ;
      Kokkos::parallel_for("KokkosBlas::Test::VanillaGEMM",
                           Kokkos::TeamPolicy<execution_space>(nRows, Kokkos::AUTO, 16), vgemm);
    }
  }

  Functor_TestBatchedSerialTrmm<DeviceType, ViewType, ScalarType, ParamTagType, Algo::Trmm::Unblocked>(alpha, A,
                                                                                                       B_actual)
      .run();

  Kokkos::fence();

  Kokkos::deep_copy(B_actual_host, B_actual);
  Kokkos::deep_copy(B_expected_host, B_expected);

  Kokkos::fence();

  // eps is ~ 10^-13 for double
  typedef typename ats::mag_type mag_type;
  const mag_type eps = 1.0e8 * ats::epsilon();
  bool fail_flag     = false;

  for (int k = 0; k < N; ++k) {
    for (int i = 0; i < nRows; ++i) {
      for (int j = 0; j < nCols; ++j) {
        if (ats::abs(B_actual_host(k, i, j) - B_expected_host(k, i, j)) > eps) {
          // printf("   Error: eps ( %g ), abs_result( %.15lf ) != abs_solution(
          // %.15lf ) (abs result-solution %g) at (k %d, i %d, j %d)\n", eps,
          // ats::abs(B_actual_host(k,i,j)), ats::abs(B_expected_host(k,i,j)),
          // ats::abs(B_actual_host(k,i,j) - B_expected_host(k,i,j)), k, i, j);
          fail_flag = true;
        }
      }
    }
  }

  ASSERT_EQ(fail_flag, false);
}
}  // namespace Trmm
}  // namespace Test

template <typename DeviceType, typename ValueType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
int test_batched_trmm(int batchSize = 512) {
  char trans = std::is_same<typename ParamTagType::trans, Trans::NoTranspose>::value     ? 'N'
               : std::is_same<typename ParamTagType::trans, Trans::Transpose>::value     ? 'T'
               : std::is_same<typename ParamTagType::trans, Trans::ConjTranspose>::value ? 'C'
                                                                                         : 'E';
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    typedef Kokkos::View<ValueType***, Kokkos::LayoutLeft, DeviceType> ViewType;
    Test::Trmm::impl_test_batched_trmm<DeviceType, ViewType, ScalarType, ParamTagType, AlgoTagType>(0, 10, 4, &trans);
    for (int i = 0; i < 10; ++i) {
      // printf("Testing: LayoutLeft,  Blksize %d\n", i);
      Test::Trmm::impl_test_batched_trmm<DeviceType, ViewType, ScalarType, ParamTagType, AlgoTagType>(batchSize, i, 4,
                                                                                                      &trans);
      Test::Trmm::impl_test_batched_trmm<DeviceType, ViewType, ScalarType, ParamTagType, AlgoTagType>(batchSize, i, 1,
                                                                                                      &trans);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    typedef Kokkos::View<ValueType***, Kokkos::LayoutRight, DeviceType> ViewType;
    Test::Trmm::impl_test_batched_trmm<DeviceType, ViewType, ScalarType, ParamTagType, AlgoTagType>(0, 10, 4, &trans);
    for (int i = 0; i < 10; ++i) {
      // printf("Testing: LayoutRight, Blksize %d\n", i);
      Test::Trmm::impl_test_batched_trmm<DeviceType, ViewType, ScalarType, ParamTagType, AlgoTagType>(batchSize, i, 4,
                                                                                                      &trans);
      Test::Trmm::impl_test_batched_trmm<DeviceType, ViewType, ScalarType, ParamTagType, AlgoTagType>(batchSize, i, 1,
                                                                                                      &trans);
    }
  }
#endif

  return 0;
}
