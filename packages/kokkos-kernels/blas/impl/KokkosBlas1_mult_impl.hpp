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
#ifndef KOKKOSBLAS1_MULT_IMPL_HPP_
#define KOKKOSBLAS1_MULT_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>

namespace KokkosBlas {
namespace Impl {

/// \brief Functor for entry-wise multiply of multivectors.
///
/// \tparam CMV 2-D Kokkos::View
/// \tparam AV 1-D Kokkos::View
/// \tparam BMV 2-D Kokkos::View
/// \tparam scalar_ab 0 if ab is zero, else nonzero (preferably 2).
/// \tparam scalar_c 0 if c is zero, else nonzero (preferably 2).
/// \tparam SizeType Index type for iterating over rows.
///
/// C(i,j) = c * C(i,j) + ab * A(i) * B(i,j), subject to the usual
/// BLAS update rules.
template <class CMV, class AV, class BMV, int scalar_ab, int scalar_c, class SizeType = typename CMV::size_type>
struct MV_MultFunctor {
  typedef SizeType size_type;
  typedef Kokkos::ArithTraits<typename CMV::non_const_value_type> ATS;

  const size_type m_n;
  typename CMV::const_value_type m_c;
  CMV m_C;
  typename AV::const_value_type m_ab;
  AV m_A;
  BMV m_B;

  MV_MultFunctor(typename CMV::const_value_type& c, const CMV& C, typename AV::const_value_type& ab, const AV& A,
                 const BMV& B)
      : m_n(C.extent(1)), m_c(c), m_C(C), m_ab(ab), m_A(A), m_B(B) {}

  KOKKOS_INLINE_FUNCTION void operator()(const size_type& i) const {
    if (scalar_c == 0) {
      if (scalar_ab == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type j = 0; j < m_n; ++j) {
          m_C(i, j) = ATS::zero();
        }
      } else {  // ab != 0, c == 0
        typename AV::const_value_type Ai = m_A(i);
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type j = 0; j < m_n; ++j) {
          m_C(i, j) = m_ab * Ai * m_B(i, j);
        }
      }
    } else {  // c != 0
      if (scalar_ab == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type j = 0; j < m_n; ++j) {
          m_C(i, j) = m_c * m_C(i, j);
        }
      } else {  // m_ab != 0, and m_c != 0
        typename AV::const_value_type Ai = m_A(i);
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type j = 0; j < m_n; ++j) {
          m_C(i, j) = m_c * m_C(i, j) + m_ab * Ai * m_B(i, j);
        }
      }
    }
  }
};

/// \brief Functor for entry-wise multiply of vectors.
///
/// \tparam CV 1-D Kokkos::View
/// \tparam AV 1-D Kokkos::View
/// \tparam BV 1-D Kokkos::View
/// \tparam scalar_ab 0 if ab is zero, else nonzero (preferably 2).
/// \tparam scalar_c 0 if c is zero, else nonzero (preferably 2).
/// \tparam SizeType Index type for iterating over rows.
///
/// C(i) = c * C(i) + ab * A(i) * B(i), subject to the usual
/// BLAS update rules.
template <class CV, class AV, class BV, int scalar_ab, int scalar_c, class SizeType = typename CV::size_type>
struct V_MultFunctor {
  typedef SizeType size_type;
  typedef Kokkos::ArithTraits<typename CV::non_const_value_type> ATS;

  typename CV::const_value_type m_c;
  CV m_C;
  typename AV::const_value_type m_ab;
  AV m_A;
  BV m_B;

  V_MultFunctor(typename CV::const_value_type& c, const CV& C, typename AV::const_value_type& ab, const AV& A,
                const BV& B)
      : m_c(c), m_C(C), m_ab(ab), m_A(A), m_B(B) {}

  KOKKOS_INLINE_FUNCTION void operator()(const size_type& i) const {
    if (scalar_c == 0) {
      if (scalar_ab == 0) {
        m_C(i) = ATS::zero();
      } else {  // ab != 0, c == 0
        m_C(i) = m_ab * m_A(i) * m_B(i);
      }
    } else {  // c != 0
      if (scalar_ab == 0) {
        m_C(i) = m_c * m_C(i);
      } else {  // m_ab != 0, and m_c != 0
        m_C(i) = m_c * m_C(i) + m_ab * m_A(i) * m_B(i);
      }
    }
  }
};

/// \brief Implementation of entry-wise multiply of vectors, that
///   dispatches to the right functor invocation.
///
/// \tparam CV 1-D Kokkos::View
/// \tparam AV 1-D Kokkos::View
/// \tparam BV 1-D Kokkos::View
/// \tparam SizeType Index type for iterating over rows.
///
/// C(i) = c * C(i) + ab * A(i) * B(i), subject to the usual BLAS
/// update rules.
template <class execution_space, class CV, class AV, class BV, class SizeType>
void V_Mult_Generic(const execution_space& space, typename CV::const_value_type& c, const CV& C,
                    typename AV::const_value_type& ab, const AV& A, const BV& B) {
  using Kokkos::ALL;
  using Kokkos::subview;
  typedef Kokkos::ArithTraits<typename AV::non_const_value_type> ATA;
  typedef Kokkos::ArithTraits<typename CV::non_const_value_type> ATC;

  const SizeType numRows = C.extent(0);
  Kokkos::RangePolicy<execution_space, SizeType> policy(space, 0, numRows);

  if (c == ATC::zero()) {
    if (ab == ATA::zero()) {
      typedef V_MultFunctor<CV, AV, BV, 0, 0, SizeType> functor_type;
      functor_type op(c, C, ab, A, B);
      Kokkos::parallel_for("KokkosBlas::Mult::S0", policy, op);
    } else {
      typedef V_MultFunctor<CV, AV, BV, 2, 0, SizeType> functor_type;
      functor_type op(c, C, ab, A, B);
      Kokkos::parallel_for("KokkosBlas::Mult::S1", policy, op);
    }
  } else {  // c != 0
    if (ab == ATA::zero()) {
      typedef V_MultFunctor<CV, AV, BV, 0, 2, SizeType> functor_type;
      functor_type op(c, C, ab, A, B);
      Kokkos::parallel_for("KokkosBlas::Mult::S2", policy, op);
    } else {
      typedef V_MultFunctor<CV, AV, BV, 2, 2, SizeType> functor_type;
      functor_type op(c, C, ab, A, B);
      Kokkos::parallel_for("KokkosBlas::Mult::S3", policy, op);
    }
  }
}

/// \brief Implementation of entry-wise multiply of multivectors, that
///   dispatches to the right functor invocation (or calls
///   V_Mult_Generic if C and B have one column).
///
/// \tparam CMV 2-D Kokkos::View
/// \tparam AV 1-D Kokkos::View
/// \tparam BMV 2-D Kokkos::View
/// \tparam SizeType Index type for iterating over rows.
///
/// C(i,j) = c * C(i,j) + ab * A(i) * B(i,j), subject to the usual
/// BLAS update rules.
template <class execution_space, class CMV, class AV, class BMV, class SizeType>
void MV_Mult_Generic(const execution_space& space, typename CMV::const_value_type& c, const CMV& C,
                     typename AV::const_value_type& ab, const AV& A, const BMV& B) {
  typedef Kokkos::ArithTraits<typename AV::non_const_value_type> ATA;
  typedef Kokkos::ArithTraits<typename CMV::non_const_value_type> ATC;

  if (C.extent(1) == 1) {
    auto C_0 = Kokkos::subview(C, Kokkos::ALL(), 0);
    auto B_0 = Kokkos::subview(B, Kokkos::ALL(), 0);
    typedef decltype(C_0) CV;
    typedef decltype(B_0) BV;

    V_Mult_Generic<execution_space, CV, AV, BV, SizeType>(space, c, C_0, ab, A, B_0);
    return;
  }

  const SizeType numRows = C.extent(0);
  Kokkos::RangePolicy<execution_space, SizeType> policy(space, 0, numRows);

  if (c == ATC::zero()) {
    if (ab == ATA::zero()) {
      typedef MV_MultFunctor<CMV, AV, BMV, 0, 0, SizeType> functor_type;
      functor_type op(c, C, ab, A, B);
      Kokkos::parallel_for("KokkosBlas::Mult::S4", policy, op);
    } else {
      typedef MV_MultFunctor<CMV, AV, BMV, 2, 0, SizeType> functor_type;
      functor_type op(c, C, ab, A, B);
      Kokkos::parallel_for("KokkosBlas::Mult::S5", policy, op);
    }
  } else {  // c != 0
    if (ab == ATA::zero()) {
      typedef MV_MultFunctor<CMV, AV, BMV, 0, 2, SizeType> functor_type;
      functor_type op(c, C, ab, A, B);
      Kokkos::parallel_for("KokkosBlas::Mult::S6", policy, op);
    } else {
      typedef MV_MultFunctor<CMV, AV, BMV, 2, 2, SizeType> functor_type;
      functor_type op(c, C, ab, A, B);
      Kokkos::parallel_for("KokkosBlas::Mult::S7", policy, op);
    }
  }
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_MULT_IMPL_HPP_
