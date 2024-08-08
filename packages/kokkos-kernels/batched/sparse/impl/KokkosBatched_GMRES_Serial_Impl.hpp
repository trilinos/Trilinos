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
#ifndef __KOKKOSBATCHED_GMRES_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_GMRES_SERIAL_IMPL_HPP__

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_Axpy.hpp"
#include "KokkosBatched_Copy_Decl.hpp"
#include "KokkosBatched_Dot.hpp"
#include "KokkosBatched_Spmv.hpp"
#include "KokkosBatched_Xpay.hpp"
#include "KokkosBatched_Givens_Serial_Internal.hpp"
#include "KokkosBatched_Trsm_Decl.hpp"
#include "KokkosBatched_Identity.hpp"
#include "KokkosBlas2_serial_gemv_impl.hpp"

namespace KokkosBatched {

///
/// Serial GMRES
///

template <typename OperatorType, typename VectorViewType, typename PrecOperatorType, typename KrylovHandleType>
KOKKOS_INLINE_FUNCTION int SerialGMRES::invoke(const OperatorType& A, const VectorViewType& _B,
                                               const VectorViewType& _X, const PrecOperatorType& P,
                                               const KrylovHandleType& handle, const int GMRES_id) {
  typedef int OrdinalType;
  typedef typename Kokkos::ArithTraits<typename VectorViewType::non_const_value_type>::mag_type MagnitudeType;
  typedef Kokkos::ArithTraits<MagnitudeType> ATM;

  using SerialCopy1D = SerialCopy<Trans::NoTranspose, 1>;
  using SerialCopy2D = SerialCopy<Trans::NoTranspose, 2>;

  const OrdinalType numMatrices = _X.extent(0);
  const OrdinalType numRows     = _X.extent(1);

  size_t maximum_iteration          = handle.get_max_iteration() < numRows ? handle.get_max_iteration() : numRows;
  const MagnitudeType tolerance     = handle.get_tolerance();
  const MagnitudeType max_tolerance = handle.get_max_tolerance();

  int n_V      = numRows;
  int n_H      = maximum_iteration + 1;
  int n_Givens = 2;

  int offset_V      = 0;
  int offset_H      = offset_V + n_V;
  int offset_Givens = offset_H + n_H;

  const int first_matrix = handle.first_index(GMRES_id);
  const int last_matrix  = handle.last_index(GMRES_id);

  auto V_view      = Kokkos::subview(handle.Arnoldi_view, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL,
                                     Kokkos::make_pair(offset_V, offset_V + n_V));
  auto H_view      = Kokkos::subview(handle.Arnoldi_view, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL,
                                     Kokkos::make_pair(offset_H, offset_H + n_H));
  auto Givens_view = Kokkos::subview(handle.Arnoldi_view, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL,
                                     Kokkos::make_pair(offset_Givens, offset_Givens + n_Givens));

  int n_G    = maximum_iteration + 1;
  int n_W    = numRows;
  int n_mask = 1;

  int offset_G    = 0;
  int offset_W    = offset_G + n_G;
  int offset_mask = offset_W + n_W;
  int offset_tmp  = offset_mask + n_mask;

  auto G    = Kokkos::subview(handle.tmp_view, Kokkos::make_pair(first_matrix, last_matrix),
                              Kokkos::make_pair(offset_G, offset_G + n_G));
  auto W    = Kokkos::subview(handle.tmp_view, Kokkos::make_pair(first_matrix, last_matrix),
                              Kokkos::make_pair(offset_W, offset_W + n_W));
  auto mask = Kokkos::subview(handle.tmp_view, Kokkos::make_pair(first_matrix, last_matrix), offset_mask);
  auto tmp  = Kokkos::subview(handle.tmp_view, Kokkos::make_pair(first_matrix, last_matrix), offset_tmp);

  // Deep copy of b into r_0:
  SerialCopy2D::invoke(_B, W);

  // r_0 := b - A x_0
  A.template apply<Trans::NoTranspose>(_X, W, -1, 1);

  P.template apply<Trans::NoTranspose, 1>(W, W);

  SerialDot<Trans::NoTranspose>::invoke(W, W, tmp);

  for (OrdinalType i = 0; i < numMatrices; ++i) {
    tmp(i) = ATM::sqrt(tmp(i));
    handle.set_norm(GMRES_id, i, 0, tmp(i));
    if (tmp(i) > max_tolerance) {
      mask(i) = 1;
      G(i, 0) = tmp(i);
      tmp(i)  = 1. / tmp(i);
    } else {
      handle.set_iteration(GMRES_id, i, 0);
      mask(i) = 0;
      G(i, 0) = 0.;
      tmp(i)  = 0.;
    }
  }

  auto V_0 = Kokkos::subview(V_view, Kokkos::ALL, 0, Kokkos::ALL);
  for (OrdinalType iRow = 0; iRow < numRows; ++iRow) {
    for (OrdinalType iMatrix = 0; iMatrix < numMatrices; ++iMatrix) {
      V_0(iMatrix, iRow) = W(iMatrix, iRow) * tmp(iMatrix);
    }
  }
  int status = 1;
  // int number_not_converged = 0;

  for (size_t j = 0; j < maximum_iteration; ++j) {
    // q := A p_j
    auto V_j = Kokkos::subview(V_view, Kokkos::ALL, j, Kokkos::ALL);

    A.template apply<Trans::NoTranspose>(V_j, W);

    P.template apply<Trans::NoTranspose, 1>(W, W);

    if (handle.get_ortho_strategy() == 0) {
      for (OrdinalType l = 0; l < numMatrices; ++l) {
        auto W_l   = Kokkos::subview(W, l, Kokkos::ALL);
        auto V_old = Kokkos::subview(V_view, l, Kokkos::make_pair(0, (int)j + 1), Kokkos::ALL);
        auto H_old = Kokkos::subview(H_view, l, j, Kokkos::make_pair(0, (int)j + 1));

        // Inner products
        KokkosBlas::SerialGemv<Trans::NoTranspose, Algo::Gemv::Unblocked>::invoke(1, V_old, W_l, 0, H_old);

        // Update
        KokkosBlas::SerialGemv<Trans::Transpose, Algo::Gemv::Unblocked>::invoke(-1, V_old, H_old, 1, W_l);
      }
    }
    if (handle.get_ortho_strategy() == 1) {
      for (size_t i = 0; i < j + 1; ++i) {
        auto V_i = Kokkos::subview(V_view, Kokkos::ALL, i, Kokkos::ALL);
        SerialDot<Trans::NoTranspose>::invoke(W, V_i, tmp);
        SerialCopy1D::invoke(tmp, Kokkos::subview(H_view, Kokkos::ALL, j, i));
        for (OrdinalType ii = 0; ii < numMatrices; ++ii) tmp(ii) = -tmp(ii);

        SerialAxpy::invoke(tmp, V_i, W);
      }
    }

    SerialDot<Trans::NoTranspose>::invoke(W, W, tmp);

    for (OrdinalType i = 0; i < numMatrices; ++i) {
      H_view(i, j, j + 1) = ATM::sqrt(tmp(i));
      tmp(i)              = H_view(i, j, j + 1) > max_tolerance ? 1. / H_view(i, j, j + 1) : 0.;
    }

    if (j + 1 < maximum_iteration) {
      auto V_n = Kokkos::subview(V_view, Kokkos::ALL, j + 1, Kokkos::ALL);
      for (OrdinalType iRow = 0; iRow < numRows; ++iRow) {
        for (OrdinalType iMatrix = 0; iMatrix < numMatrices; ++iMatrix) {
          V_n(iMatrix, iRow) = W(iMatrix, iRow) * tmp(iMatrix);
        }
      }
    }

    for (OrdinalType l = 0; l < numMatrices; ++l) {
      // Apply the previous Givens rotations:
      auto H_j        = Kokkos::subview(H_view, l, j, Kokkos::ALL);
      auto Givens_0_l = Kokkos::subview(Givens_view, l, Kokkos::ALL, 0);
      auto Givens_1_l = Kokkos::subview(Givens_view, l, Kokkos::ALL, 1);

      if (mask(l) == 1.) {
        for (size_t i = 0; i < j; ++i) {
          auto tmp1  = Givens_0_l(i) * H_j(i) + Givens_1_l(i) * H_j(i + 1);
          auto tmp2  = -Givens_1_l(i) * H_j(i) + Givens_0_l(i) * H_j(i + 1);
          H_j(i)     = tmp1;
          H_j(i + 1) = tmp2;
        }

        // Compute the new Givens rotation:
        Kokkos::pair<typename VectorViewType::non_const_value_type, typename VectorViewType::non_const_value_type>
            G_new(1, 0);
        typename VectorViewType::non_const_value_type alpha = 0;
        SerialGivensInternal::invoke(H_j(j), H_j(j + 1), &G_new, &alpha);

        Givens_0_l(j) = G_new.first;
        Givens_1_l(j) = G_new.second;

        // Apply the new Givens rotation:
        auto tmp1  = Givens_0_l(j) * H_j(j) + Givens_1_l(j) * H_j(j + 1);
        auto tmp2  = -Givens_1_l(j) * H_j(j) + Givens_0_l(j) * H_j(j + 1);
        H_j(j)     = tmp1;
        H_j(j + 1) = tmp2;

        G(l, j + 1) = -Givens_1_l(j) * G(l, j);
        G(l, j) *= Givens_0_l(j);
      } else {
        H_j(j)      = 1.;
        G(l, j + 1) = 0.;
      }

      auto res_norm = Kokkos::ArithTraits<double>::abs(G(l, j + 1)) / G(l, 0);

      handle.set_norm(GMRES_id, l, j + 1, res_norm);

      if (mask(l) == 1. && res_norm < tolerance) {
        mask(l)     = 0.;
        G(l, j + 1) = 0.;
        handle.set_iteration(GMRES_id, l, j + 1);
      }
    }

    bool all_converged = true;
    for (OrdinalType l = 0; l < numMatrices; ++l) all_converged = (all_converged && mask(l) == 0.);
    if (all_converged) {
      maximum_iteration = j + 1;
      break;
    }
  }

  auto first_indices = Kokkos::make_pair(0, (int)maximum_iteration);

  for (OrdinalType l = 0; l < numMatrices; ++l) {
    auto A_l = Kokkos::subview(H_view, l, first_indices, first_indices);
    auto B_l = Kokkos::subview(G, l, first_indices);

    SerialTrsm<Side::Left, Uplo::Lower, Trans::Transpose, Diag::NonUnit, Algo::Trsm::Unblocked>::invoke(1, A_l, B_l);
  }

  if (handle.get_ortho_strategy() == 0) {
    for (OrdinalType l = 0; l < numMatrices; ++l) {
      KokkosBlas::SerialGemv<Trans::Transpose, Algo::Gemv::Unblocked>::invoke(
          1, Kokkos::subview(V_view, l, first_indices, Kokkos::ALL), Kokkos::subview(G, l, first_indices), 1,
          Kokkos::subview(_X, l, Kokkos::ALL));
    }
  }
  if (handle.get_ortho_strategy() == 1) {
    for (size_t j = 0; j < maximum_iteration; ++j) {
      SerialAxpy::invoke(Kokkos::subview(G, Kokkos::ALL, j), Kokkos::subview(V_view, Kokkos::ALL, j, Kokkos::ALL), _X);
    }
  }

  if (handle.get_compute_last_residual()) {
    SerialCopy2D::invoke(_B, W);
    A.template apply<Trans::NoTranspose>(_X, W, -1, 1);
    P.template apply<Trans::NoTranspose, 1>(W, W);
    SerialDot<Trans::NoTranspose>::invoke(W, W, tmp);

    for (OrdinalType i = 0; i < numMatrices; ++i) {
      tmp(i) = ATM::sqrt(tmp(i));
      handle.set_last_norm(GMRES_id, i, tmp(i));
    }
  }
  return status;
}

template <typename OperatorType, typename VectorViewType, typename KrylovHandleType>
KOKKOS_INLINE_FUNCTION int SerialGMRES::invoke(const OperatorType& A, const VectorViewType& _B,
                                               const VectorViewType& _X, const KrylovHandleType& handle) {
  Identity P;
  return invoke<OperatorType, VectorViewType, Identity>(A, _B, _X, P, handle);
}
}  // namespace KokkosBatched

#endif
