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
#ifndef KOKKOS_BLAS1_IMPL_ABS_HPP_
#define KOKKOS_BLAS1_IMPL_ABS_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>

namespace KokkosBlas {
namespace Impl {

//
// abs
//

// Entry-wise absolute value / magnitude: R(i,j) = abs(X(i,j)).
template <class RMV, class XMV, class SizeType = typename RMV::size_type>
struct MV_Abs_Functor {
  typedef SizeType size_type;
  typedef Kokkos::ArithTraits<typename XMV::non_const_value_type> ATS;

  const size_type numCols;
  RMV R_;
  XMV X_;

  MV_Abs_Functor(const RMV& R, const XMV& X) : numCols(X.extent(1)), R_(R), X_(X) {
    static_assert(Kokkos::is_view<RMV>::value,
                  "KokkosBlas::Impl::"
                  "MV_Abs_Functor: RMV is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::"
                  "MV_Abs_Functor: XMV is not a Kokkos::View.");
    static_assert(RMV::rank == 2,
                  "KokkosBlas::Impl::"
                  "MV_Abs_Functor: RMV is not rank 2");
    static_assert(XMV::rank == 2,
                  "KokkosBlas::Impl::"
                  "MV_Abs_Functor: XMV is not rank 2");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type& i) const {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for (size_type j = 0; j < numCols; ++j) {
      R_(i, j) = ATS::abs(X_(i, j));
    }
  }
};

// Entry-wise, in-place absolute value / magnitude: R(i,j) = abs(R(i,j)).
template <class RMV, class SizeType = typename RMV::size_type>
struct MV_AbsSelf_Functor {
  typedef SizeType size_type;
  typedef Kokkos::ArithTraits<typename RMV::non_const_value_type> ATS;

  const size_type numCols;
  RMV R_;

  MV_AbsSelf_Functor(const RMV& R) : numCols(R.extent(1)), R_(R) {
    static_assert(Kokkos::is_view<RMV>::value,
                  "KokkosBlas::Impl::"
                  "MV_Abs_Functor: RMV is not a Kokkos::View.");
    static_assert(RMV::rank == 2,
                  "KokkosBlas::Impl::"
                  "MV_Abs_Functor: RMV is not rank 2");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type& i) const {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for (size_type j = 0; j < numCols; ++j) {
      R_(i, j) = ATS::abs(R_(i, j));
    }
  }
};

// Single-vector, entry-wise absolute value / magnitude: R(i) = abs(X(i)).
template <class RV, class XV, class SizeType = typename RV::size_type>
struct V_Abs_Functor {
  typedef SizeType size_type;
  typedef Kokkos::ArithTraits<typename XV::non_const_value_type> ATS;

  RV R_;
  XV X_;

  V_Abs_Functor(const RV& R, const XV& X) : R_(R), X_(X) {
    static_assert(Kokkos::is_view<RV>::value,
                  "KokkosBlas::Impl::"
                  "V_Abs_Functor: RV is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XV>::value,
                  "KokkosBlas::Impl::"
                  "V_Abs_Functor: XV is not a Kokkos::View.");
    static_assert(RV::rank == 1,
                  "KokkosBlas::Impl::"
                  "V_Abs_Functor: RV is not rank 1");
    static_assert(XV::rank == 1,
                  "KokkosBlas::Impl::"
                  "V_Abs_Functor: XV is not rank 1");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type& i) const { R_(i) = ATS::abs(X_(i)); }
};

// Single-vector, entry-wise, in-place absolute value / magnitude: R(i) =
// abs(R(i)).
template <class RV, class SizeType = typename RV::size_type>
struct V_AbsSelf_Functor {
  typedef SizeType size_type;
  typedef Kokkos::ArithTraits<typename RV::non_const_value_type> ATS;

  RV R_;

  V_AbsSelf_Functor(const RV& R) : R_(R) {
    static_assert(Kokkos::is_view<RV>::value,
                  "KokkosBlas::Impl::"
                  "V_Abs_Functor: RV is not a Kokkos::View.");
    static_assert(RV::rank == 1,
                  "KokkosBlas::Impl::"
                  "V_Abs_Functor: RV is not rank 1");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type& i) const { R_(i) = ATS::abs(R_(i)); }
};

// Invoke the "generic" (not unrolled) multivector functor that
// computes entry-wise absolute value.
template <class execution_space, class RMV, class XMV, class SizeType>
void MV_Abs_Generic(const execution_space& space, const RMV& R, const XMV& X) {
  static_assert(Kokkos::is_view<RMV>::value,
                "KokkosBlas::Impl::"
                "MV_Abs_Generic: RMV is not a Kokkos::View.");
  static_assert(Kokkos::is_view<XMV>::value,
                "KokkosBlas::Impl::"
                "MV_Abs_Generic: XMV is not a Kokkos::View.");
  static_assert(RMV::rank == 2,
                "KokkosBlas::Impl::"
                "MV_Abs_Generic: RMV is not rank 2");
  static_assert(XMV::rank == 2,
                "KokkosBlas::Impl::"
                "MV_Abs_Generic: XMV is not rank 2");

  const SizeType numRows = X.extent(0);
  Kokkos::RangePolicy<execution_space, SizeType> policy(space, 0, numRows);

  if ((void*)(R.data()) == (void*)(X.data())) {  // if R and X are the same (alias one another)
    MV_AbsSelf_Functor<RMV, SizeType> op(R);
    Kokkos::parallel_for("KokkosBlas::Abs::S0", policy, op);
  } else {
    MV_Abs_Functor<RMV, XMV, SizeType> op(R, X);
    Kokkos::parallel_for("KokkosBlas::Abs::S1", policy, op);
  }
}

// Variant of MV_Abs_Generic for single vectors (1-D Views) R and X.
template <class execution_space, class RV, class XV, class SizeType>
void V_Abs_Generic(const execution_space& space, const RV& R, const XV& X) {
  static_assert(Kokkos::is_view<RV>::value,
                "KokkosBlas::Impl::"
                "V_Abs_Generic: RV is not a Kokkos::View.");
  static_assert(Kokkos::is_view<XV>::value,
                "KokkosBlas::Impl::"
                "V_Abs_Generic: XV is not a Kokkos::View.");
  static_assert(RV::rank == 1,
                "KokkosBlas::Impl::"
                "V_Abs_Generic: RV is not rank 1");
  static_assert(XV::rank == 1,
                "KokkosBlas::Impl::"
                "V_Abs_Generic: XV is not rank 1");

  const SizeType numRows = X.extent(0);
  Kokkos::RangePolicy<execution_space, SizeType> policy(space, 0, numRows);

  if ((void*)(R.data()) == (void*)(X.data())) {  // if R and X are the same (alias one another)
    V_AbsSelf_Functor<RV, SizeType> op(R);
    Kokkos::parallel_for("KokkosBlas::Abs::S2", policy, op);
  } else {
    V_Abs_Functor<RV, XV, SizeType> op(R, X);
    Kokkos::parallel_for("KokkosBlas::Abs::S3", policy, op);
  }
}

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOS_BLAS1_MV_IMPL_ABS_HPP_
