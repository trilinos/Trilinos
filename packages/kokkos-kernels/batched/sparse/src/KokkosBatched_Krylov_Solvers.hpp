// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
