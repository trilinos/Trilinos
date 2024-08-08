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
#ifndef KOKKOSBLAS1_SCAL_IMPL_HPP_
#define KOKKOSBLAS1_SCAL_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>
#include <KokkosBlas1_scal_spec.hpp>

#ifndef KOKKOSBLAS_OPTIMIZATION_LEVEL_SCAL
#define KOKKOSBLAS_OPTIMIZATION_LEVEL_SCAL 2
#endif  // KOKKOSBLAS_OPTIMIZATION_LEVEL_SCAL

namespace KokkosBlas {
namespace Impl {

// Single-vector version of MV_Scal_Functor.  By default, a is still a
// 1-D View.  Below is a partial specialization that lets a be a
// scalar.  This functor computes any of the following:
//
// 1. Y(i) = alpha*X(i) for alpha in -1,0,1
// 2. Y(i) = a(0)*X(i)
//
// The template parameter scalar_x corresponds to alpha in the
// operation y = alpha*x + beta*y.  The values -1, 0, and -1
// correspond to literal values of this coefficient.  The value 2
// tells the functor to use the corresponding vector of coefficients.
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding (multi)vector entry.  This does not apply to
// coefficients in the a vector, if used.
template <class RV, class AV, class XV, int scalar_x, class SizeType>
struct V_Scal_Functor {
  typedef SizeType size_type;
  typedef Kokkos::ArithTraits<typename RV::non_const_value_type> ATS;

  RV m_r;
  XV m_x;
  AV m_a;

  V_Scal_Functor(const RV& r, const XV& x, const AV& a, const SizeType startingColumn) : m_r(r), m_x(x), m_a(a) {
    static_assert(Kokkos::is_view<RV>::value, "V_Scal_Functor: RV is not a Kokkos::View.");
    static_assert(Kokkos::is_view<AV>::value, "V_Scal_Functor: AV is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XV>::value, "V_Scal_Functor: XV is not a Kokkos::View.");
    static_assert(RV::rank == 1, "V_Scal_Functor: RV is not rank 1.");
    static_assert(AV::rank == 1, "V_Scal_Functor: AV is not rank 1.");
    static_assert(XV::rank == 1, "V_Scal_Functor: XV is not rank 1.");

    if (startingColumn != 0) {
      m_a = Kokkos::subview(a, std::make_pair(startingColumn, static_cast<SizeType>(a.extent(0))));
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type& i) const {
    // scalar_x is a compile-time constant (since it is a template
    // parameter), so the compiler should evaluate these branches at
    // compile time.
    if (scalar_x == 0) {
      m_r(i) = ATS::zero();
    }
    if (scalar_x == -1) {
      m_r(i) = -m_x(i);
    }
    if (scalar_x == 1) {
      m_r(i) = m_x(i);
    }
    if (scalar_x == 2) {
      m_r(i) = m_a(0) * m_x(i);
    }
  }
};

// Partial specialization of V_Scal_Functor that lets a be a scalar
// (rather than a 1-D View, as in the most general version above).
// This functor computes any of the following:
//
// 1. Y(i) = alpha*X(i) for alpha in -1,0,1
// 2. Y(i) = a*X(i)
template <class RV, class XV, int scalar_x, class SizeType>
struct V_Scal_Functor<RV, typename XV::non_const_value_type, XV, scalar_x, SizeType> {
  typedef SizeType size_type;
  typedef Kokkos::ArithTraits<typename RV::non_const_value_type> ATS;

  RV m_r;
  XV m_x;
  const typename XV::non_const_value_type m_a;

  V_Scal_Functor(const RV& r, const XV& x, const typename XV::non_const_value_type& a,
                 const SizeType /* startingColumn */)
      : m_r(r), m_x(x), m_a(a) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type& i) const {
    if (scalar_x == 0) {
      m_r(i) = ATS::zero();
    }
    if (scalar_x == -1) {
      m_r(i) = -m_x(i);
    }
    if (scalar_x == 1) {
      m_r(i) = m_x(i);
    }
    if (scalar_x == 2) {
      m_r(i) = m_a * m_x(i);
    }
  }
};

// Variant of MV_Scal_Generic for single vectors (1-D Views) r and x.
// As above, av is either a 1-D View (and only its first entry will be
// read), or a scalar.
template <class execution_space, class RV, class AV, class XV, class SizeType>
void V_Scal_Generic(const execution_space& space, const RV& r, const AV& av, const XV& x, const SizeType startingColumn,
                    int a = 2) {
  static_assert(Kokkos::is_view<RV>::value, "V_Scal_Generic: RV is not a Kokkos::View.");
  static_assert(Kokkos::is_view<XV>::value, "V_Scal_Generic: XV is not a Kokkos::View.");
  static_assert(RV::rank == 1, "V_Scal_Generic: RV is not rank 1.");
  static_assert(XV::rank == 1, "V_Scal_Generic: XV is not rank 1.");

  const SizeType numRows = x.extent(0);
  Kokkos::RangePolicy<execution_space, SizeType> policy(space, 0, numRows);

  if (a == 0) {
    V_Scal_Functor<RV, AV, XV, 0, SizeType> op(r, x, av, startingColumn);
    Kokkos::parallel_for("KokkosBlas::Scal::S0", policy, op);
    return;
  }
  if (a == -1) {
    V_Scal_Functor<RV, AV, XV, -1, SizeType> op(r, x, av, startingColumn);
    Kokkos::parallel_for("KokkosBlas::Scal::S1", policy, op);
    return;
  }
  if (a == 1) {
    V_Scal_Functor<RV, AV, XV, 1, SizeType> op(r, x, av, startingColumn);
    Kokkos::parallel_for("KokkosBlas::Scal::S2", policy, op);
    return;
  }

  // a arbitrary (not -1, 0, or 1)
  V_Scal_Functor<RV, AV, XV, 2, SizeType> op(r, x, av, startingColumn);
  Kokkos::parallel_for("KokkosBlas::Scal::S3", policy, op);
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_SCAL_IMPL_HPP_
