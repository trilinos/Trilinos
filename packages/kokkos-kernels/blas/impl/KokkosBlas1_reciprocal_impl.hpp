// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBLAS1_IMPL_RECIPROCAL_HPP_
#define KOKKOSBLAS1_IMPL_RECIPROCAL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_ArithTraits.hpp>

namespace KokkosBlas {
namespace Impl {

//
// reciprocal
//

#if defined(__GNUC__) && ((__GNUC__ == 12) || (__GNUC__ == 13)) && !defined(__NVCC__)
#define KOKKOSKERNELS_GCC_BUG_DISABLE_TREE_VECTORIZE 1
#else
#define KOKKOSKERNELS_GCC_BUG_DISABLE_TREE_VECTORIZE 0
#endif

// Entry-wise reciprocalolute value / magnitude: R(i,j) = reciprocal(X(i,j)).
template <class RMV, class XMV, class SizeType = typename RMV::size_type>
struct MV_Reciprocal_Functor {
  typedef SizeType size_type;
  typedef KokkosKernels::ArithTraits<typename XMV::non_const_value_type> ATS;

  const size_type numCols;
  RMV R_;
  XMV X_;

  MV_Reciprocal_Functor(const RMV& R, const XMV& X) : numCols(X.extent(1)), R_(R), X_(X) {
    static_assert(Kokkos::is_view<RMV>::value,
                  "KokkosBlas::Impl::"
                  "MV_Reciprocal_Functor: RMV is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::"
                  "MV_Reciprocal_Functor: XMV is not a Kokkos::View.");
    static_assert(RMV::rank == 2,
                  "KokkosBlas::Impl::"
                  "MV_Reciprocal_Functor: RMV is not rank 2");
    static_assert(XMV::rank == 2,
                  "KokkosBlas::Impl::"
                  "MV_Reciprocal_Functor: XMV is not rank 2");
  }

  // disable vectorization in this function
  // work-around https://github.com/kokkos/kokkos-kernels/issues/2091
#if KOKKOSKERNELS_GCC_BUG_DISABLE_TREE_VECTORIZE
#pragma GCC push_options
#pragma GCC optimize("no-tree-vectorize")
#endif
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type& i) const {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for (size_type j = 0; j < numCols; ++j) {
      R_(i, j) = ATS::one() / X_(i, j);
    }
  }
};
#if KOKKOSKERNELS_GCC_BUG_DISABLE_TREE_VECTORIZE
#pragma GCC pop_options
#endif

// Entry-wise, in-place reciprocalolute value / magnitude: R(i,j) =
// reciprocal(R(i,j)).
template <class RMV, class SizeType = typename RMV::size_type>
struct MV_ReciprocalSelf_Functor {
  typedef SizeType size_type;
  typedef KokkosKernels::ArithTraits<typename RMV::non_const_value_type> ATS;

  const size_type numCols;
  RMV R_;

  MV_ReciprocalSelf_Functor(const RMV& R) : numCols(R.extent(1)), R_(R) {
    static_assert(Kokkos::is_view<RMV>::value,
                  "KokkosBlas::Impl::"
                  "MV_Reciprocal_Functor: RMV is not a Kokkos::View.");
    static_assert(RMV::rank == 2,
                  "KokkosBlas::Impl::"
                  "MV_Reciprocal_Functor: RMV is not rank 2");
  }

  // disable vectorization in this function
  // work-around https://github.com/kokkos/kokkos-kernels/issues/2091
#if KOKKOSKERNELS_GCC_BUG_DISABLE_TREE_VECTORIZE
#pragma GCC push_options
#pragma GCC optimize("no-tree-vectorize")
#endif
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type& i) const {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for (size_type j = 0; j < numCols; ++j) {
      R_(i, j) = ATS::one() / R_(i, j);
    }
  }
};
#if KOKKOSKERNELS_GCC_BUG_DISABLE_TREE_VECTORIZE
#pragma GCC pop_options
#endif

// Single-vector, entry-wise reciprocalolute value / magnitude: R(i) =
// reciprocal(X(i)).
template <class RV, class XV, class SizeType = typename RV::size_type>
struct V_Reciprocal_Functor {
  typedef SizeType size_type;
  typedef KokkosKernels::ArithTraits<typename XV::non_const_value_type> ATS;

  RV R_;
  XV X_;

  V_Reciprocal_Functor(const RV& R, const XV& X) : R_(R), X_(X) {
    static_assert(Kokkos::is_view<RV>::value,
                  "KokkosBlas::Impl::"
                  "V_Reciprocal_Functor: RV is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XV>::value,
                  "KokkosBlas::Impl::"
                  "V_Reciprocal_Functor: XV is not a Kokkos::View.");
    static_assert(RV::rank == 1,
                  "KokkosBlas::Impl::"
                  "V_Reciprocal_Functor: RV is not rank 1");
    static_assert(XV::rank == 1,
                  "KokkosBlas::Impl::"
                  "V_Reciprocal_Functor: XV is not rank 1");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type& i) const { R_(i) = ATS::one() / X_(i); }
};

// Single-vector, entry-wise, in-place reciprocalolute value / magnitude: R(i) =
// reciprocal(R(i)).
template <class RV, class SizeType = typename RV::size_type>
struct V_ReciprocalSelf_Functor {
  typedef SizeType size_type;
  typedef KokkosKernels::ArithTraits<typename RV::non_const_value_type> ATS;

  RV R_;

  V_ReciprocalSelf_Functor(const RV& R) : R_(R) {
    static_assert(Kokkos::is_view<RV>::value,
                  "KokkosBlas::Impl::"
                  "V_Reciprocal_Functor: RV is not a Kokkos::View.");
    static_assert(RV::rank == 1,
                  "KokkosBlas::Impl::"
                  "V_Reciprocal_Functor: RV is not rank 1");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type& i) const { R_(i) = ATS::one() / R_(i); }
};

// Invoke the "generic" (not unrolled) multivector functor that
// computes entry-wise reciprocalolute value.
template <class execution_space, class RMV, class XMV, class SizeType>
void MV_Reciprocal_Generic(const execution_space& space, const RMV& R, const XMV& X) {
  static_assert(Kokkos::is_view<RMV>::value,
                "KokkosBlas::Impl::"
                "MV_Reciprocal_Generic: RMV is not a Kokkos::View.");
  static_assert(Kokkos::is_view<XMV>::value,
                "KokkosBlas::Impl::"
                "MV_Reciprocal_Generic: XMV is not a Kokkos::View.");
  static_assert(RMV::rank == 2,
                "KokkosBlas::Impl::"
                "MV_Reciprocal_Generic: RMV is not rank 2");
  static_assert(XMV::rank == 2,
                "KokkosBlas::Impl::"
                "MV_Reciprocal_Generic: XMV is not rank 2");

  const SizeType numRows = X.extent(0);
  Kokkos::RangePolicy<execution_space, SizeType> policy(space, 0, numRows);

  if (R == X) {  // if R and X are the same (alias one another)
    MV_ReciprocalSelf_Functor<RMV, SizeType> op(R);
    Kokkos::parallel_for("KokkosBlas::Reciprocal::S0", policy, op);
  } else {
    MV_Reciprocal_Functor<RMV, XMV, SizeType> op(R, X);
    Kokkos::parallel_for("KokkosBlas::Reciprocal::S1", policy, op);
  }
}

// Variant of MV_Reciprocal_Generic for single vectors (1-D Views) R and X.
template <class execution_space, class RV, class XV, class SizeType>
void V_Reciprocal_Generic(const execution_space& space, const RV& R, const XV& X) {
  static_assert(Kokkos::is_view<RV>::value,
                "KokkosBlas::Impl::"
                "V_Reciprocal_Generic: RV is not a Kokkos::View.");
  static_assert(Kokkos::is_view<XV>::value,
                "KokkosBlas::Impl::"
                "V_Reciprocal_Generic: XV is not a Kokkos::View.");
  static_assert(RV::rank == 1,
                "KokkosBlas::Impl::"
                "V_Reciprocal_Generic: RV is not rank 1");
  static_assert(XV::rank == 1,
                "KokkosBlas::Impl::"
                "V_Reciprocal_Generic: XV is not rank 1");

  const SizeType numRows = X.extent(0);
  Kokkos::RangePolicy<execution_space, SizeType> policy(space, 0, numRows);

  if (R == X) {  // if R and X are the same (alias one another)
    V_ReciprocalSelf_Functor<RV, SizeType> op(R);
    Kokkos::parallel_for("KokkosBlas::Reciprocal::S2", policy, op);
  } else {
    V_Reciprocal_Functor<RV, XV, SizeType> op(R, X);
    Kokkos::parallel_for("KokkosBlas::Reciprocal::S3", policy, op);
  }
}

#undef KOKKOSKERNELS_GCC_BUG_DISABLE_TREE_VECTORIZE

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSBLAS1_IMPL_RECIPROCAL_HPP_
