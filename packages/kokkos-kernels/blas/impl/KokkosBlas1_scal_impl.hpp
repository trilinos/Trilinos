// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBLAS1_SCAL_IMPL_HPP_
#define KOKKOSBLAS1_SCAL_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_InnerProductSpaceTraits.hpp>
#include <KokkosBlas1_axpby_unification.hpp>
#include <KokkosBlas1_scal_spec.hpp>

#ifndef KOKKOSBLAS_OPTIMIZATION_LEVEL_SCAL
#define KOKKOSBLAS_OPTIMIZATION_LEVEL_SCAL 2
#endif  // KOKKOSBLAS_OPTIMIZATION_LEVEL_SCAL

namespace KokkosBlas {
namespace Impl {

// Single-vector functor for R(i) = alpha * X(i).
//
// AV may be a scalar, a rank-0 View, or a rank-1 View (in which case only
// element 0 is used, after any startingColumn adjustment).
// getCoefficient() dispatches uniformly across all three cases, so no
// partial specialisations are required.
//
// scalar_x encodes the compile-time shortcut:
//   -1 → negate, 0 → zero, 1 → copy, 2 → general multiply.
template <class RV, class AV, class XV, int scalar_x, class SizeType>
struct V_Scal_Functor {
  typedef SizeType size_type;
  typedef KokkosKernels::ArithTraits<typename RV::non_const_value_type> ATS;

  RV m_r;
  XV m_x;
  AV m_a;

  V_Scal_Functor(const RV& r, const XV& x, const AV& a, const SizeType startingColumn) : m_r(r), m_x(x), m_a(a) {
    static_assert(Kokkos::is_view<RV>::value, "V_Scal_Functor: RV is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XV>::value, "V_Scal_Functor: XV is not a Kokkos::View.");
    static_assert(RV::rank == 1, "V_Scal_Functor: RV is not rank 1.");
    static_assert(XV::rank == 1, "V_Scal_Functor: XV is not rank 1.");

    if constexpr (isRank1View<AV>()) {
      if (startingColumn != 0) {
        m_a = Kokkos::subview(a, std::make_pair(startingColumn, static_cast<SizeType>(a.extent(0))));
      }
    }
  }

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
      m_r(i) = getCoefficient(m_a) * m_x(i);
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
