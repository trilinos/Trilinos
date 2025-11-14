// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_NONPIV_LDL_ON_DEVICE_HPP__
#define __TACHO_NONPIV_LDL_ON_DEVICE_HPP__

/// \file  Tacho_NonPivLDL_OnDevice.hpp
/// \brief BLAS general matrix matrix multiplication
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Copy_OnDevice.hpp"
#include "Tacho_Trsm_OnDevice.hpp"
#include "Tacho_Scale2x2_BlockInverseDiagonals_OnDevice.hpp"
#include "Tacho_GemmTriangular_OnDevice.hpp"
#include "Tacho_Gemm_Internal.hpp"

namespace Tacho {

template <typename ArgUplo> struct LDL_nopiv<ArgUplo, Algo::OnDevice> {

  // Right-look column-wise
  template <typename MemberType, typename ViewTypeA>
  inline static int invoke_col(MemberType &member, const ViewTypeA &A) {

    using exec_space = MemberType;
    using policy_type = Kokkos::RangePolicy<exec_space>;
    using value_type = typename ViewTypeA::non_const_value_type;
    using arith_traits = ArithTraits<value_type>;

    const auto &exec_instance = member;
    const ordinal_type m = A.extent(0);

    int r_val(0);
    for (ordinal_type i = 0; i < m; i++) {
      ////for (ordinal_type j = i+1; j < m; j++) A(i, j) /= A(i, i);
      const auto policy_scale = policy_type(exec_instance, i+1, m);
      Kokkos::parallel_for(policy_scale, KOKKOS_LAMBDA(const ordinal_type &j) {
          A(i, j) /= A(i, i); });

      ////for (ordinal_type k = i+1; k < m; k++) {
      ////  for (ordinal_type j = k; j < m; j++) A(k, j) -= A(i, k) * A(i, i) * A(i, j);
      ////}
      #if 1
       const ordinal_type mn = m - (i+1);
       policy_type policy_update(exec_instance, 0, mn*mn);
       Kokkos::parallel_for(policy_update, KOKKOS_LAMBDA(const ordinal_type &id) {
         ordinal_type k = (i+1) + id / mn;
         ordinal_type j = (i+1) + id % mn;
         A(k, j) -= arith_traits::conj(A(i, k)) * A(i, i) * A(i, j);
       });
      #else
       policy_type policy_update(exec_instance, i+1, m);
       Kokkos::parallel_for(policy_update, KOKKOS_LAMBDA(const ordinal_type &k) {
         for (ordinal_type j = k; j < m; j++) A(k, j) -= arith_traits::conj(A(i, k)) * A(i, i) * A(i, j);
       });
      #endif
    }
    return r_val;
  }


  template <typename HandleType, typename MemberType, typename ViewTypeA, typename ViewTypeW>
  inline static int invoke(HandleType &handle, MemberType &member, const ViewTypeA &A, const ViewTypeW &W) {

    using value_type = typename ViewTypeA::non_const_value_type;
    using range_type = Kokkos::pair<ordinal_type, ordinal_type>;

    const value_type  one ( 1.0);
    const value_type mone (-1.0);
    const ordinal_type m  = A.extent(0);
    char *nb_env = getenv("TACHO_BLOCK_SIZE");
    const ordinal_type nb = (nb_env == NULL ? 256 : atoi(nb_env));

    int r_val(0);
    for (ordinal_type b = 0; b < m; b+=nb) {
      ordinal_type i1 =  b;
      ordinal_type i2 = (b+nb < m ? b+nb : m);
      ordinal_type m1 = i2-i1;

      // Factorize Diagonal Block
      auto A11 = Kokkos::subview(A, range_type(i1, i2), range_type(i1, i2));
      invoke_col(member, A11);

      ordinal_type m2 = m - i2;
      if (m2 > 0) {
        // Compute off-diagonal blocks
        auto A12 = Kokkos::subview(A, range_type(i1, i2), range_type(i2, m));
        Trsm<Side::Left, Uplo::Upper, Trans::ConjTranspose, Algo::OnDevice>::invoke(
                  handle, Diag::Unit(), one, A11, A12);

        // Save A12 in workspace
        UnmanagedViewType<ViewTypeA> T(W.data(), m1, m2);
        Copy<Algo::OnDevice>::invoke(member, T, A12);

        // Apply D^{-1} on off-diagonal
        Scale_BlockInverseDiagonals<Side::Left, Algo::OnDevice>::invoke(member, A11, A12);

        // A22 = -A12*T
        auto A22 = Kokkos::subview(A, range_type(i2, m), range_type(i2, m));
#if 0
        const auto &exec_instance = member;

        ordinal_type nb2 = m2;
        const ordinal_type num_blks = (m2+nb2-1)/nb2;
        using team_policy_type = Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space>;
        team_policy_type team_policy(exec_instance, num_blks, Kokkos::AUTO());
        Kokkos::parallel_for(
          team_policy, KOKKOS_LAMBDA(const typename team_policy_type::member_type &team_member) {
            ordinal_type id = team_member.league_rank();
            ordinal_type k1 = nb2*id;
            ordinal_type k2 = (k1+nb2 < m2 ? k1+nb2 : m2);
            auto Tk = Kokkos::subview(T,   Kokkos::ALL(), range_type(k1, k2));
            auto Ck = Kokkos::subview(A22, Kokkos::ALL(), range_type(k1, k2));
            Gemm<Trans::Transpose, Trans::NoTranspose, Algo::Internal>::invoke(
                      team_member, mone, A12, Tk, one, Ck);
          });
#else
        GemmTriangular<Trans::Transpose, Trans::NoTranspose, Uplo::Upper, Algo::OnDevice>::invoke(
                    handle, mone, A12, T, one, A22);
#endif
      }
    }
    return r_val;
  }
};

} // namespace Tacho
#endif
