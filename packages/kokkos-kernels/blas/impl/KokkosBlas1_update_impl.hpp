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
#ifndef KOKKOSBLAS1_UPDATE_IMPL_HPP_
#define KOKKOSBLAS1_UPDATE_IMPL_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_InnerProductSpaceTraits.hpp"

namespace KokkosBlas {
namespace Impl {

//
// update
//

// Functor for multivectors X, Y, and Z, that computes
//
// Z(i,j) = alpha*X(i,j) + beta*Y(i,j) + gamma*Z(i,j)
//
// with special cases for alpha, beta, or gamma = 0.
//
// The template parameters scalar_x, scalar_y, and scalar_z correspond
// to alpha, beta, resp. gammar in the operation Z = alpha*X + beta*Y
// + gamma*Z.  The value 0 corresponds to literal values of those
// coefficients.  The value 2 tells the functor to use the
// corresponding input coefficient.  Any literal coefficient of zero
// has BLAS semantics of ignoring the corresponding (multi)vector
// entry.
template <class XMV, class YMV, class ZMV, int scalar_x, int scalar_y, int scalar_z,
          class SizeType = typename ZMV::size_type>
struct MV_Update_Functor {
  typedef SizeType size_type;
  typedef Kokkos::ArithTraits<typename ZMV::non_const_value_type> ATS;

  const size_type numCols;
  const typename XMV::non_const_value_type alpha_;
  XMV X_;
  const typename YMV::non_const_value_type beta_;
  YMV Y_;
  const typename ZMV::non_const_value_type gamma_;
  ZMV Z_;

  MV_Update_Functor(const typename XMV::non_const_value_type& alpha, const XMV& X,
                    const typename YMV::non_const_value_type& beta, const YMV& Y,
                    const typename ZMV::non_const_value_type& gamma, const ZMV& Z)
      : numCols(X.extent(1)), alpha_(alpha), X_(X), beta_(beta), Y_(Y), gamma_(gamma), Z_(Z) {
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::"
                  "MV_Update_Functor: X is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YMV>::value,
                  "KokkosBlas::Impl::"
                  "MV_Update_Functor: Y is not a Kokkos::View.");
    static_assert(Kokkos::is_view<ZMV>::value,
                  "KokkosBlas::Impl::"
                  "MV_Update_Functor: Z is not a Kokkos::View.");
    static_assert(std::is_same<typename ZMV::value_type, typename ZMV::non_const_value_type>::value,
                  "KokkosBlas::Impl::MV_Update_Functor: Z is const.  "
                  "It must be nonconst, because it is an output argument "
                  "(we have to be able to write to its entries).");
    // Casting enum values to int avoids compiler warnings about
    // comparing different kinds of enum values.
    static_assert((int)ZMV::rank == (int)XMV::rank && (int)ZMV::rank == (int)YMV::rank,
                  "KokkosBlas::Impl::MV_Update_Functor: "
                  "X, Y, and Z must have the same rank.");
    static_assert(ZMV::rank == 2,
                  "KokkosBlas::Impl::MV_Update_Functor: "
                  "XMV, YMV, and ZMV must have rank 2.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type& i) const {
    // scalar_x, scalar_y, and scalar_z are compile-time constants
    // (since they are template parameters), so the compiler should
    // evaluate these branches at compile time.
    if (scalar_x == 0) {
      if (scalar_y == 0) {
        if (scalar_z == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            Z_(i, k) = ATS::zero();
          }
        } else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            Z_(i, k) = gamma_ * Z_(i, k);
          }
        }
      } else {
        if (scalar_z == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            Z_(i, k) = beta_ * Y_(i, k);
          }
        } else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            Z_(i, k) = beta_ * Y_(i, k) + gamma_ * Z_(i, k);
          }
        }
      }
    }
    //
    // scalar_x == 2
    //
    else {
      if (scalar_y == 0) {
        if (scalar_z == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            Z_(i, k) = alpha_ * X_(i, k);
          }
        } else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            Z_(i, k) = alpha_ * X_(i, k) + gamma_ * Z_(i, k);
          }
        }
      } else {
        if (scalar_z == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            Z_(i, k) = alpha_ * X_(i, k) + beta_ * Y_(i, k);
          }
        } else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            Z_(i, k) = alpha_ * X_(i, k) + beta_ * Y_(i, k) + gamma_ * Z_(i, k);
          }
        }
      }
    }
  }
};

// Functor for vectors X, Y, and Z, that computes
//
// Z(i) = alpha*X(i) + beta*Y(i) + gamma*Z(i)
//
// with special cases for alpha, beta, or gamma = 0.
//
// The template parameters scalar_x, scalar_y, and scalar_z correspond
// to alpha, beta, resp. gammar in the operation Z = alpha*X + beta*Y
// + gamma*Z.  The value 0 corresponds to literal values of those
// coefficients.  The value 2 tells the functor to use the
// corresponding input coefficient.  Any literal coefficient of zero
// has BLAS semantics of ignoring the corresponding vector entry.
template <class XV, class YV, class ZV, int scalar_x, int scalar_y, int scalar_z,
          class SizeType = typename ZV::size_type>
struct V_Update_Functor {
  typedef SizeType size_type;
  typedef Kokkos::ArithTraits<typename ZV::non_const_value_type> ATS;

  const size_type numCols;
  const typename XV::non_const_value_type alpha_;
  XV X_;
  const typename YV::non_const_value_type beta_;
  YV Y_;
  const typename ZV::non_const_value_type gamma_;
  ZV Z_;

  V_Update_Functor(const typename XV::non_const_value_type& alpha, const XV& X,
                   const typename YV::non_const_value_type& beta, const YV& Y,
                   const typename ZV::non_const_value_type& gamma, const ZV& Z)
      : numCols(X.extent(1)), alpha_(alpha), X_(X), beta_(beta), Y_(Y), gamma_(gamma), Z_(Z) {
    static_assert(Kokkos::is_view<XV>::value,
                  "KokkosBlas::Impl::"
                  "V_Update_Functor: X is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YV>::value,
                  "KokkosBlas::Impl::"
                  "V_Update_Functor: Y is not a Kokkos::View.");
    static_assert(Kokkos::is_view<ZV>::value,
                  "KokkosBlas::Impl::"
                  "V_Update_Functor: Z is not a Kokkos::View.");
    static_assert(std::is_same<typename ZV::value_type, typename ZV::non_const_value_type>::value,
                  "KokkosBlas::Impl::V_Update_Functor: Z is const.  "
                  "It must be nonconst, because it is an output argument "
                  "(we have to be able to write to its entries).");
    // Casting to int avoids compiler warnings about comparing
    // different kinds of enum values.
    static_assert((int)ZV::rank == (int)XV::rank && (int)ZV::rank == (int)YV::rank,
                  "KokkosBlas::Impl::V_Update_Functor: "
                  "X, Y, and Z must have the same rank.");
    static_assert(ZV::rank == 1,
                  "KokkosBlas::Impl::V_Update_Functor: "
                  "XV, YV, and ZV must have rank 1.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type& i) const {
    // scalar_x, scalar_y, and scalar_z are compile-time constants
    // (since they are template parameters), so the compiler should
    // evaluate these branches at compile time.
    if (scalar_x == 0) {
      if (scalar_y == 0) {
        if (scalar_z == 0) {
          Z_(i) = ATS::zero();
        } else {
          Z_(i) = gamma_ * Z_(i);
        }
      } else {
        if (scalar_z == 0) {
          Z_(i) = beta_ * Y_(i);
        } else {
          Z_(i) = beta_ * Y_(i) + gamma_ * Z_(i);
        }
      }
    }
    //
    // scalar_ x == 2
    //
    else {
      if (scalar_y == 0) {
        if (scalar_z == 0) {
          Z_(i) = alpha_ * X_(i);
        } else {
          Z_(i) = alpha_ * X_(i) + gamma_ * Z_(i);
        }
      } else {
        if (scalar_z == 0) {
          Z_(i) = alpha_ * X_(i) + beta_ * Y_(i);
        } else {
          Z_(i) = alpha_ * X_(i) + beta_ * Y_(i) + gamma_ * Z_(i);
        }
      }
    }
  }
};

// Invoke the "generic" (not unrolled) multivector functor that
// computes
//
// Z(i,j) = alpha*X(i,j) + beta*Y(i,j) + gamma*Z(i,j)
//
// with special cases for alpha, beta, or gamma = 0.
//
// a, b, and c come in as integers.  The value 0 corresponds to the
// literal values of the coefficients.  The value 2 tells the functor
// to use the corresponding coefficients: a == 2 means use alpha, b ==
// 2 means use beta, and c == 2 means use gamma.  Otherwise, the
// corresponding coefficients are ignored.
//
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding multivector entry.
template <class execution_space, class XMV, class YMV, class ZMV, class SizeType>
void MV_Update_Generic(const execution_space& space, const typename XMV::non_const_value_type& alpha, const XMV& X,
                       const typename YMV::non_const_value_type& beta, const YMV& Y,
                       const typename ZMV::non_const_value_type& gamma, const ZMV& Z, int a = 2, int b = 2, int c = 2) {
  static_assert(Kokkos::is_view<XMV>::value,
                "KokkosBlas::Impl::"
                "MV_Update_Generic: X is not a Kokkos::View.");
  static_assert(Kokkos::is_view<YMV>::value,
                "KokkosBlas::Impl::"
                "MV_Update_Generic: Y is not a Kokkos::View.");
  static_assert(Kokkos::is_view<ZMV>::value,
                "KokkosBlas::Impl::"
                "MV_Update_Generic: Z is not a Kokkos::View.");
  static_assert(std::is_same<typename ZMV::value_type, typename ZMV::non_const_value_type>::value,
                "KokkosBlas::Impl::MV_Update_Generic: Z is const.  "
                "It must be nonconst, because it is an output argument "
                "(we have to be able to write to its entries).");
  // Casting to int avoids compiler warnings about comparing different
  // kinds of enum values.
  static_assert((int)ZMV::rank == (int)XMV::rank && (int)ZMV::rank == (int)YMV::rank,
                "KokkosBlas::Impl::MV_Update_Generic: "
                "X, Y, and Z must have the same rank.");
  static_assert(ZMV::rank == 2,
                "KokkosBlas::Impl::MV_Update_Generic: "
                "XMV, YMV, and ZMV must have rank 2.");

  const SizeType numRows = X.extent(0);
  Kokkos::RangePolicy<execution_space, SizeType> policy(space, 0, numRows);

  if (a == 0) {
    if (b == 0) {
      if (c == 0) {
        MV_Update_Functor<XMV, YMV, ZMV, 0, 0, 0, SizeType> op(alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for("KokkosBlas::update<MV,0,0,0>", policy, op);
      } else {
        MV_Update_Functor<XMV, YMV, ZMV, 0, 0, 2, SizeType> op(alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for("KokkosBlas::update<MV,0,0,c>", policy, op);
      }
    } else {
      if (c == 0) {
        MV_Update_Functor<XMV, YMV, ZMV, 0, 2, 0, SizeType> op(alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for("KokkosBlas::update<MV,0,b,0>", policy, op);
      } else {
        MV_Update_Functor<XMV, YMV, ZMV, 0, 2, 2, SizeType> op(alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for("KokkosBlas::update<MV,0,b,c>", policy, op);
      }
    }
  }
  //
  // a == 2
  //
  else {
    if (b == 0) {
      if (c == 0) {
        MV_Update_Functor<XMV, YMV, ZMV, 2, 0, 0, SizeType> op(alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for("KokkosBlas::update<MV,a,0,0>", policy, op);
      } else {
        MV_Update_Functor<XMV, YMV, ZMV, 2, 0, 2, SizeType> op(alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for("KokkosBlas::update<MV,a,0,c>", policy, op);
      }
    } else {
      if (c == 0) {
        MV_Update_Functor<XMV, YMV, ZMV, 2, 2, 0, SizeType> op(alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for("KokkosBlas::update<MV,a,b,0>", policy, op);
      } else {
        MV_Update_Functor<XMV, YMV, ZMV, 2, 2, 2, SizeType> op(alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for("KokkosBlas::update<MV,a,b,c>", policy, op);
      }
    }
  }
}

// Invoke the "generic" (not unrolled) single-vector functor that
// computes
//
// Z(i) = alpha*X(i) + beta*Y(i) + gamma*Z(i)
//
// with special cases for alpha, beta, or gamma = 0.
//
// a, b, and c come in as integers.  The value 0 corresponds to the
// literal values of the coefficients.  The value 2 tells the functor
// to use the corresponding coefficients: a == 2 means use alpha, b ==
// 2 means use beta, and c == 2 means use gamma.  Otherwise, the
// corresponding coefficients are ignored.
//
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding vector entry.
template <class execution_space, class XV, class YV, class ZV, class SizeType>
void V_Update_Generic(const execution_space& space, const typename XV::non_const_value_type& alpha, const XV& X,
                      const typename YV::non_const_value_type& beta, const YV& Y,
                      const typename ZV::non_const_value_type& gamma, const ZV& Z, int a = 2, int b = 2, int c = 2) {
  static_assert(Kokkos::is_view<XV>::value,
                "KokkosBlas::Impl::"
                "V_Update_Generic: X is not a Kokkos::View.");
  static_assert(Kokkos::is_view<YV>::value,
                "KokkosBlas::Impl::"
                "V_Update_Generic: Y is not a Kokkos::View.");
  static_assert(Kokkos::is_view<ZV>::value,
                "KokkosBlas::Impl::"
                "V_Update_Generic: Z is not a Kokkos::View.");
  static_assert(std::is_same<typename ZV::value_type, typename ZV::non_const_value_type>::value,
                "KokkosBlas::Impl::V_Update_Generic: Z is const.  "
                "It must be nonconst, because it is an output argument "
                "(we have to be able to write to its entries).");
  // Casting to int avoids compiler warnings about comparing
  // different kinds of enum values.
  static_assert((int)ZV::rank == (int)XV::rank && (int)ZV::rank == (int)YV::rank,
                "KokkosBlas::Impl::V_Update_Generic: "
                "X, Y, and Z must have the same rank.");
  static_assert(ZV::rank == 1,
                "KokkosBlas::Impl::V_Update_Generic: "
                "XV, YV, and ZV must have rank 1.");

  const SizeType numRows = X.extent(0);
  Kokkos::RangePolicy<execution_space, SizeType> policy(space, 0, numRows);

  if (a == 0) {
    if (b == 0) {
      if (c == 0) {
        V_Update_Functor<XV, YV, ZV, 0, 0, 0, SizeType> op(alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for("KokkosBlas::update<0,0,0>", policy, op);
      } else {
        V_Update_Functor<XV, YV, ZV, 0, 0, 2, SizeType> op(alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for("KokkosBlas::update<0,0,c>", policy, op);
      }
    } else {
      if (c == 0) {
        V_Update_Functor<XV, YV, ZV, 0, 2, 0, SizeType> op(alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for("KokkosBlas::update<0,b,0>", policy, op);
      } else {
        V_Update_Functor<XV, YV, ZV, 0, 2, 2, SizeType> op(alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for("KokkosBlas::update<0,b,c>", policy, op);
      }
    }
  }
  //
  // a == 2
  //
  else {
    if (b == 0) {
      if (c == 0) {
        V_Update_Functor<XV, YV, ZV, 2, 0, 0, SizeType> op(alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for("KokkosBlas::update<a,0,0>", policy, op);
      } else {
        V_Update_Functor<XV, YV, ZV, 2, 0, 2, SizeType> op(alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for("KokkosBlas::update<a,0,c>", policy, op);
      }
    } else {
      if (c == 0) {
        V_Update_Functor<XV, YV, ZV, 2, 2, 0, SizeType> op(alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for("KokkosBlas::update<a,b,0>", policy, op);
      } else {
        V_Update_Functor<XV, YV, ZV, 2, 2, 2, SizeType> op(alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for("KokkosBlas::update<a,b,c>", policy, op);
      }
    }
  }
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_UPDATE_IMPL_HPP_
