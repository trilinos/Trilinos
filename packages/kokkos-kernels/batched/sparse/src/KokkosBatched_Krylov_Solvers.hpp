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

#ifndef KOKKOSBATCHED_KRYLOV_SOLVERS_HPP
#define KOKKOSBATCHED_KRYLOV_SOLVERS_HPP

namespace KokkosBatched {

struct SerialGMRES {
  template <typename OperatorType, typename VectorViewType, typename PrecOperatorType, typename KrylovHandleType>
  KOKKOS_INLINE_FUNCTION static int invoke(const OperatorType& A, const VectorViewType& B, const VectorViewType& X,
                                           const PrecOperatorType& P, const KrylovHandleType& handle,
                                           const int GMRES_id);
  template <typename OperatorType, typename VectorViewType, typename KrylovHandleType>
  KOKKOS_INLINE_FUNCTION static int invoke(const OperatorType& A, const VectorViewType& B, const VectorViewType& X,
                                           const KrylovHandleType& handle);
};

template <typename MemberType>
struct TeamGMRES {
  template <typename OperatorType, typename VectorViewType, typename PrecOperatorType, typename KrylovHandleType,
            typename ArnoldiViewType, typename TMPViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const OperatorType& A, const VectorViewType& B,
                                           const VectorViewType& X, const PrecOperatorType& P,
                                           const KrylovHandleType& handle, const ArnoldiViewType& ArnoldiView,
                                           const TMPViewType& TMPView);
  template <typename OperatorType, typename VectorViewType, typename PrecOperatorType, typename KrylovHandleType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const OperatorType& A, const VectorViewType& B,
                                           const VectorViewType& X, const PrecOperatorType& P,
                                           const KrylovHandleType& handle);
  template <typename OperatorType, typename VectorViewType, typename KrylovHandleType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const OperatorType& A, const VectorViewType& B,
                                           const VectorViewType& X, const KrylovHandleType& handle);
};

template <typename MemberType>
struct TeamVectorGMRES {
  template <typename OperatorType, typename VectorViewType, typename PrecOperatorType, typename KrylovHandleType,
            typename ArnoldiViewType, typename TMPViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const OperatorType& A, const VectorViewType& B,
                                           const VectorViewType& X, const PrecOperatorType& P,
                                           const KrylovHandleType& handle, const ArnoldiViewType& ArnoldiView,
                                           const TMPViewType& TMPView);
  template <typename OperatorType, typename VectorViewType, typename PrecOperatorType, typename KrylovHandleType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const OperatorType& A, const VectorViewType& B,
                                           const VectorViewType& X, const PrecOperatorType& P,
                                           const KrylovHandleType& handle);
  template <typename OperatorType, typename VectorViewType, typename KrylovHandleType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const OperatorType& A, const VectorViewType& B,
                                           const VectorViewType& X, const KrylovHandleType& handle);
};

template <typename MemberType>
struct TeamCG {
  template <typename OperatorType, typename VectorViewType, typename KrylovHandleType, typename TMPViewType,
            typename TMPNormViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const OperatorType& A, const VectorViewType& B,
                                           const VectorViewType& X, const KrylovHandleType& handle,
                                           const TMPViewType& TMPView, const TMPNormViewType& TMPNormView);
  template <typename OperatorType, typename VectorViewType, typename KrylovHandleType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const OperatorType& A, const VectorViewType& B,
                                           const VectorViewType& X, const KrylovHandleType& handle);
};

template <typename MemberType>
struct TeamVectorCG {
  template <typename OperatorType, typename VectorViewType, typename KrylovHandleType, typename TMPViewType,
            typename TMPNormViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const OperatorType& A, const VectorViewType& B,
                                           const VectorViewType& X, const KrylovHandleType& handle,
                                           const TMPViewType& TMPView, const TMPNormViewType& TMPNormView);
  template <typename OperatorType, typename VectorViewType, typename KrylovHandleType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const OperatorType& A, const VectorViewType& B,
                                           const VectorViewType& X, const KrylovHandleType& handle);
};

}  // namespace KokkosBatched

#endif
