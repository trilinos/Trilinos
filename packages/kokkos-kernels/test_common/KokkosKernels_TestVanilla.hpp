// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSKERNELS_TESTVANILLA_HPP
#define KOKKOSKERNELS_TESTVANILLA_HPP

#include <random>

#include "KokkosKernels_Utils.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_ArithTraits.hpp"
#include "KokkosBatched_Vector.hpp"

namespace Test {

template <class ViewTypeA, class ViewTypeB, class ViewTypeC, class ExecutionSpace>
struct SharedVanillaGEMM {
  bool A_t, B_t, A_c, B_c;
  int C_rows, C_cols, A_cols;
  ViewTypeA A_;
  ViewTypeB B_;
  ViewTypeC C_;

  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;
  typedef typename ViewTypeC::value_type ScalarC;
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, typename ViewTypeA::device_type> SubviewTypeA;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutStride, typename ViewTypeB::device_type> SubviewTypeB;
  typedef KokkosKernels::ArithTraits<ScalarC> APT;
  typedef typename APT::mag_type mag_type;
  ScalarA alpha;
  ScalarC beta;

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<ExecutionSpace>::member_type& team) const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, C_rows), [&](const int& i) {
      // Give each kokkos thread a vector of A
      SubviewTypeA a_vec;
      if (A_t)
        a_vec = Kokkos::subview(A_, Kokkos::ALL(), i);
      else
        a_vec = Kokkos::subview(A_, i, Kokkos::ALL());

      // Have all vector lanes perform the dot product
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, C_cols), [&](const int& j) {
        SubviewTypeB b_vec;
        if (B_t)
          b_vec = Kokkos::subview(B_, j, Kokkos::ALL());
        else
          b_vec = Kokkos::subview(B_, Kokkos::ALL(), j);
        ScalarC ab = ScalarC(0);
        for (int k = 0; k < A_cols; k++) {
          auto a = A_c ? APT::conj(a_vec(k)) : a_vec(k);
          auto b = B_c ? APT::conj(b_vec(k)) : b_vec(k);
          ab += a * b;
        }
        C_(i, j) = beta * C_(i, j) + alpha * ab;
      });
    });
  }
};
// C(i,:,:) = alpha * (A(i,:,:) * B(i,:,:)) + beta * C(i,:,:)
template <class ViewTypeA, class ViewTypeB, class ViewTypeC, class ExecutionSpace>
struct Functor_BatchedVanillaGEMM {
  bool A_t, B_t, A_c, B_c, batch_size_last_dim = false;
  ViewTypeA A_;
  ViewTypeB B_;
  ViewTypeC C_;

  using ScalarA      = typename ViewTypeA::value_type;
  using ScalarB      = typename ViewTypeB::value_type;
  using ScalarC      = typename ViewTypeC::value_type;
  using SubviewTypeA = typename Kokkos::View<ScalarA**, Kokkos::LayoutStride, typename ViewTypeA::device_type>;
  using SubviewTypeB = typename Kokkos::View<ScalarB**, Kokkos::LayoutStride, typename ViewTypeA::device_type>;
  using SubviewTypeC = typename Kokkos::View<ScalarC**, Kokkos::LayoutStride, typename ViewTypeA::device_type>;

  ScalarA alpha;
  ScalarC beta;

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<ExecutionSpace>::member_type& team) const {
    int i = team.league_rank();
    SubviewTypeA A;
    SubviewTypeB B;
    SubviewTypeC C;

    if (batch_size_last_dim) {
      A = Kokkos::subview(A_, Kokkos::ALL(), Kokkos::ALL(), i);
      B = Kokkos::subview(B_, Kokkos::ALL(), Kokkos::ALL(), i);
      C = Kokkos::subview(C_, Kokkos::ALL(), Kokkos::ALL(), i);
    } else {
      A = Kokkos::subview(A_, i, Kokkos::ALL(), Kokkos::ALL());
      B = Kokkos::subview(B_, i, Kokkos::ALL(), Kokkos::ALL());
      C = Kokkos::subview(C_, i, Kokkos::ALL(), Kokkos::ALL());
    }
    struct SharedVanillaGEMM<SubviewTypeA, SubviewTypeB, SubviewTypeC, ExecutionSpace> vgemm;
    vgemm.A_t    = A_t;
    vgemm.B_t    = B_t;
    vgemm.A_c    = A_c;
    vgemm.B_c    = B_c;
    vgemm.C_rows = batch_size_last_dim ? C_.extent(0) : C_.extent(1);
    vgemm.C_cols = batch_size_last_dim ? C_.extent(1) : C_.extent(2);
    vgemm.A_cols = batch_size_last_dim ? (A_t ? A_.extent(0) : A_.extent(1)) : (A_t ? A_.extent(1) : A_.extent(2));
    vgemm.A_     = A;
    vgemm.B_     = B;
    vgemm.C_     = C;
    vgemm.alpha  = alpha;
    vgemm.beta   = beta;
    vgemm(team);
  }

  inline void run() {
    Kokkos::parallel_for(
        "Test::VanillaGEMM",
        Kokkos::TeamPolicy<ExecutionSpace>(batch_size_last_dim ? C_.extent(2) : C_.extent(0), Kokkos::AUTO,
                                           KokkosKernels::Impl::kk_get_max_vector_size<ExecutionSpace>()),
        *this);
  }
};

// Compute C := alpha * AB + beta * C
template <class ViewTypeA, class ViewTypeB, class ViewTypeC>
void vanillaGEMM(typename ViewTypeC::non_const_value_type alpha, const ViewTypeA& A, const ViewTypeB& B,
                 typename ViewTypeC::non_const_value_type beta, const ViewTypeC& C) {
  using value_type = typename ViewTypeC::non_const_value_type;
  using KAT        = KokkosKernels::ArithTraits<value_type>;
  int m            = A.extent(0);
  int k            = A.extent(1);
  int n            = B.extent(1);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      value_type sum = KAT::zero();
      for (int ii = 0; ii < k; ii++) {
        sum += A(i, ii) * B(ii, j);
      }
      C(i, j) = alpha * sum + beta * C(i, j);
    }
  }
}

template <class AlphaType, class ViewTypeA, class ViewTypeX, class BetaType, class ViewTypeY>
KOKKOS_INLINE_FUNCTION void vanillaGEMV(char mode, AlphaType alpha, const ViewTypeA& A, const ViewTypeX& x,
                                        BetaType beta, const ViewTypeY& y) {
  using ScalarY         = typename ViewTypeY::non_const_value_type;
  using KAT_A           = KokkosKernels::ArithTraits<typename ViewTypeA::non_const_value_type>;
  const bool transposed = mode == 'T' || mode == 'C';
  const bool conjugated = mode == 'C';
  const bool has_beta   = beta != KokkosKernels::ArithTraits<BetaType>::zero();
  int M                 = A.extent(transposed ? 1 : 0);
  int N                 = A.extent(transposed ? 0 : 1);
  for (int i = 0; i < M; i++) {
    ScalarY y_i{};
    if (has_beta) y_i = beta * y(i);
    for (int j = 0; j < N; j++) {
      const auto a   = transposed ? A(j, i) : A(i, j);
      const auto Aij = conjugated ? KAT_A::conj(a) : a;
      y_i += alpha * Aij * x(j);
    }
    y(i) = y_i;
  }
}

}  // namespace Test
#endif
