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

// Trsm with deficient-diagonals
template <typename ArgSide, typename ArgUplo, typename ArgTransA>
struct Trsm_defs<ArgSide, ArgUplo, ArgTransA, Algo::OnDevice> {
  template <typename MemberType, typename DiagType, typename ScalarType, typename ViewTypeA, typename ViewTypeB>
  inline static int invoke(MemberType &member, const DiagType diagA, const ScalarType alpha, const ViewTypeA &A,
                           const ViewTypeB &B) {
#if 1
    using exec_space = MemberType;
    using policy_type = Kokkos::RangePolicy<exec_space>;
    using value_type = typename ViewTypeA::non_const_value_type;
    const auto &exec_instance = member;
    const value_type zero (0);
    const value_type one (1);

    const ordinal_type m = B.extent(0);
    const ordinal_type n = B.extent(1);
    // Side::Left, Uplo::Upper, Trans::Transpose
    if (ArgSide::param != 'L' || ArgUplo::param != 'U' || ArgTransA::param != 'T') 
      printf( " Trsm_defs(%c,%c,%c) not implemented\n",ArgSide::param,ArgUplo::param,ArgTransA::param );
    const auto policy_scale = policy_type(exec_instance, 0, m);
    for (ordinal_type i = 0; i < m; i++) {
      Kokkos::parallel_for(policy_scale, KOKKOS_LAMBDA(const ordinal_type &j) {
        if (A(i, i) == zero ) {
          // if tiny pivot, zero out off-diagonal
          B(i, j) = zero;
        } else {
          if (diagA.param != 'U') {
            // scale
            B(i, j) /= A(i, i);
          }
          // update
          for (ordinal_type k=j+1; k<m; k++) B(k, j) -= A(k, i) * B(i, j);
        }
      });
      // reset zero-pivot with one
      // TODO: move it out, reset after TRSM with off-diagonal blocks
      //if (A(i, i) == zero ) A(i, i) = one;
    }
#endif
    return 0;
  }
};

// LDL without pivoting 
template <typename ArgUplo> struct LDL_nopiv<ArgUplo, Algo::OnDevice> {

  // LDL of a diagonal block (Right-look column-wise)
  template <typename MemberType, typename ViewTypeA, typename ViewTypeR>
  inline static int invoke_col(MemberType &member, const double tol, const ViewTypeA &A, const ViewTypeR &rval_d) {

    using exec_space = MemberType;
    using policy_type = Kokkos::RangePolicy<exec_space>;
    using value_type = typename ViewTypeA::non_const_value_type;
    using arith_traits = ArithTraits<value_type>;
    using range_type = Kokkos::pair<ordinal_type, ordinal_type>;

    const auto &exec_instance = member;
    const bool conjugate = false;
    const ordinal_type m = A.extent(0);
    const value_type zero (0);

    int r_val(0);
    Kokkos::deep_copy(rval_d, 0);
    for (ordinal_type i = 0; i < m; i++) {
      const ordinal_type mn = m - (i+1);
      const auto policy_scale  = policy_type(exec_instance, i+1, m);
      const auto policy_update = policy_type(exec_instance, 0, mn*mn);
      const auto policy_diag   = policy_type(exec_instance, i, i+1);

      if (tol > 0.0) {
        // == Factor ith row (!! Checking for tiny pivot !!) ==
        // scaling the diagonal
        Kokkos::parallel_for(policy_scale, KOKKOS_LAMBDA(const ordinal_type &j) {
          if (arith_traits::abs(A(i, i)) < tol) {
            // if tiny pivot, zero-out off-diagonal (assuming we hit null-space)
            A(i, j) = zero; 
          } else {
            // scale with diagonal
            A(i, j) /= A(i, i);
          }
        });

        // update trailing submatrix
        Kokkos::parallel_for(policy_update, KOKKOS_LAMBDA(const ordinal_type &id) {
          // skip update if tiny pivot
          if (arith_traits::abs(A(i, i)) >= tol) {
            ordinal_type k = (i+1) + id / mn;
            ordinal_type j = (i+1) + id % mn;
            A(k, j) -= (conjugate ? arith_traits::conj(A(i, k)) : A(i, k)) * A(i, i) * A(i, j);
          }
        });

        // mark tiny diagonal with zero
        // NOTE: any better way (using single)?
        Kokkos::parallel_for(policy_diag, KOKKOS_LAMBDA(const ordinal_type &id) {
          if (arith_traits::abs(A(id, id)) < tol) {
              A(id, id) = zero;
              rval_d(0) ++;
          }
        });
      } else {
        // == Factor ith row (!! Without checking for tiny pivot (note: nan if zero pivot) !!) ==
        // scaling the off-diagonal in the i-th row of the upper-triangular matrix
        Kokkos::parallel_for(policy_scale, KOKKOS_LAMBDA(const ordinal_type &j) {
            A(i, j) /= A(i, i); });

        // update trailing submatrix
        Kokkos::parallel_for(policy_update, KOKKOS_LAMBDA(const ordinal_type &id) {
          ordinal_type k = (i+1) + id / mn;
          ordinal_type j = (i+1) + id % mn;
          A(k, j) -= (conjugate ? arith_traits::conj(A(i, k)) : A(i, k)) * A(i, i) * A(i, j);
        });
      }
    }
    // move return-value from device to host
    auto rval_h = Kokkos::create_mirror_view(rval_d);
    Kokkos::deep_copy(rval_h, rval_d);
    r_val = rval_h(0);
    return r_val;
  }

  // main factorization code
  template <typename HandleType, typename MemberType, typename ViewTypeA, typename ViewTypeW, typename ViewTypeR>
  inline static int invoke(HandleType &handle, MemberType &member, const double tol, const ViewTypeA &A, const ViewTypeW &W, const ViewTypeR &r_val_view) {

    using value_type = typename ViewTypeA::non_const_value_type;
    using range_type = Kokkos::pair<ordinal_type, ordinal_type>;

    const value_type zero ( 0);
    const value_type  one ( 1);
    const value_type mone (-1);
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
      int ndefs = invoke_col(member, tol, A11, r_val_view);

      ordinal_type m2 = m - i2;
      if (m2 > 0) {
        // Compute off-diagonal blocks (no conjugate)
        auto A12 = Kokkos::subview(A, range_type(i1, i2), range_type(i2, m));
        if (ndefs > 0) {
          // Trsm while skipping zero-diagonals (indicating deficiencies/null-space)
          //  and then replace zero-diagonals with one
          Trsm_defs<Side::Left, Uplo::Upper, Trans::Transpose, Algo::OnDevice>::invoke(
                    member, Diag::Unit(), one, A11, A12);
        } else {
          Trsm<Side::Left, Uplo::Upper, Trans::Transpose, Algo::OnDevice>::invoke(
                    handle, Diag::Unit(), one, A11, A12);
        }
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
      r_val += ndefs;
    }
    return r_val;
  }

  // reset zero diags
  template <typename HandleType, typename MemberType, typename ViewTypeA>
  inline static void reset_zero_diags(HandleType &handle, MemberType &member, const ViewTypeA &A) {
    using value_type = typename ViewTypeA::non_const_value_type;
    using exec_space = MemberType;
    using policy_type = Kokkos::RangePolicy<exec_space>;

    const ordinal_type m = A.extent(0);
    const value_type zero ( 0);
    const value_type  one ( 1);

    const auto &exec_instance = member;
    const auto policy_scale = policy_type(exec_instance, 0, m);
    Kokkos::parallel_for(policy_scale, KOKKOS_LAMBDA(const ordinal_type &i) {
      if (A(i, i) == zero ) A(i, i) = one;
    });
  }
};

} // namespace Tacho
#endif
