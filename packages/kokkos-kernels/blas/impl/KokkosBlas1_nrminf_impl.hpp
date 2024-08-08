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
#ifndef KOKKOSBLAS1_NRMINF_IMPL_HPP_
#define KOKKOSBLAS1_NRMINF_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <KokkosBlas1_nrminf_spec.hpp>

namespace KokkosBlas {
namespace Impl {

//
// nrminf_squared
//

/// \brief 2-norm (squared) functor for single vectors.
///
/// \tparam RV 0-D output View
/// \tparam XV 1-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template <class RV, class XV, class SizeType = typename XV::size_type>
struct V_NrmInf_Functor {
  typedef typename XV::execution_space execution_space;
  typedef SizeType size_type;
  typedef typename XV::non_const_value_type xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef Kokkos::ArithTraits<typename IPT::mag_type> AT;
  typedef typename IPT::mag_type value_type;

  typename XV::const_type m_x;

  V_NrmInf_Functor(const XV& x) : m_x(x) {
    static_assert(Kokkos::is_view<RV>::value,
                  "KokkosBlas::Impl::V_NrmInf_Functor: "
                  "R is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XV>::value,
                  "KokkosBlas::Impl::V_NrmInf_Functor: "
                  "X is not a Kokkos::View.");
    static_assert(std::is_same<typename RV::value_type, typename RV::non_const_value_type>::value,
                  "KokkosBlas::Impl::V_NrmInf_Functor: R is const.  "
                  "It must be nonconst, because it is an output argument "
                  "(we have to be able to write to its entries).");
    static_assert(RV::rank == 0 && XV::rank == 1,
                  "KokkosBlas::Impl::V_NrmInf_Functor: "
                  "RV must have rank 0 and XV must have rank 1.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type& i, value_type& max) const {
    value_type val = IPT::norm(m_x(i));
    if (val > max) max = val;
  }
};

/// \brief Compute the 2-norm (or its square) of the single vector (1-D
///   View) X, and store the result in the 0-D View r.
template <class execution_space, class RV, class XV, class SizeType>
void V_NrmInf_Invoke(const execution_space& space, const RV& r, const XV& X) {
  typedef Kokkos::ArithTraits<typename RV::non_const_value_type> AT;

  const SizeType numRows = static_cast<SizeType>(X.extent(0));

  // Avoid Max Reduction if this is a zero length view
  if (numRows == 0) {
    Kokkos::deep_copy(r, AT::zero());
    return;
  }

  Kokkos::RangePolicy<execution_space, SizeType> policy(space, 0, numRows);

  typedef V_NrmInf_Functor<RV, XV, SizeType> functor_type;
  functor_type op(X);
  Kokkos::parallel_reduce("KokkosBlas::NrmInf::S0", policy, op, Kokkos::Max<typename RV::non_const_value_type>(r()));
}

/// \brief Compute the 2-norms (or their square) of the columns of the
///   multivector (2-D View) X, and store result(s) in the 1-D View r.
template <class execution_space, class RV, class XMV, class SizeType>
void MV_NrmInf_Invoke(const execution_space& space, const RV& r, const XMV& X) {
  for (size_t i = 0; i < X.extent(1); i++) {
    auto ri = Kokkos::subview(r, i);
    auto Xi = Kokkos::subview(X, Kokkos::ALL(), i);
    V_NrmInf_Invoke<execution_space, decltype(ri), decltype(Xi), SizeType>(space, ri, Xi);
  }
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_NRMINF_IMPL_HPP_
