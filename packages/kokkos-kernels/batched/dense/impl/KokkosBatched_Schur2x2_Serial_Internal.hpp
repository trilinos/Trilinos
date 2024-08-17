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
#ifndef __KOKKOSBATCHED_SCHUR2X2_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_SCHUR2X2_SERIAL_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
///
/// this impl follows the flame interface of householder transformation
///
struct SerialSchur2x2Internal {
  template <typename RealType>
  KOKKOS_INLINE_FUNCTION static int invoke(RealType* alpha00, RealType* alpha01, RealType* alpha10, RealType* alpha11,
                                           Kokkos::pair<RealType, RealType>* G, Kokkos::complex<RealType>* lambda1,
                                           Kokkos::complex<RealType>* lambda2, bool* is_complex) {
    typedef RealType real_type;
    typedef Kokkos::ArithTraits<real_type> ats;
    const real_type zero(0), one(1), half(0.5), minus_one(-1);
    /// compute G = [ gamma -sigma;
    ///               sigma  gamma ];
    /// G.first = gamma and G.second = sigma
    /// this rotation satisfy the following
    ///   G' [alpha00 alpha01   G = [ beta00 beta01;
    ///       alpha10 alpha11 ]       beta10 beta11 ];
    /// where either
    ///   1) beta00 = beta11 and beta01*beta10 < 0
    ///   2) beta10 = 0
    const real_type tol = ats::epsilon() * real_type(100);
    if (ats::abs(*alpha10) < tol) {
      /// no rotation
      *G = Kokkos::pair<real_type, real_type>(one, zero);
      /// two real eigen values
      *lambda1    = Kokkos::complex<real_type>(*alpha00, zero);
      *lambda2    = Kokkos::complex<real_type>(*alpha11, zero);
      *is_complex = false;
    } else if (ats::abs(*alpha01) < tol) {
      /// 90 degree rotation (permutation)
      *G = Kokkos::pair<real_type, real_type>(zero, one);
      /// [ 0 1 ][alpha00 0       [ 0 -1  --> [ alpha11 -alpha10
      ///  -1 0 ] alpha10 alpha11]  1  0]       0        alpha00]
      const real_type tmp = *alpha00;
      *alpha00            = *alpha11;
      *alpha11            = tmp;
      *alpha01            = -(*alpha10);
      *alpha10            = zero;
      /// two real eigen values
      *lambda1    = Kokkos::complex<real_type>(*alpha00, zero);
      *lambda2    = Kokkos::complex<real_type>(*alpha11, zero);
      *is_complex = false;
    } else if (ats::abs(*alpha00 - *alpha11) < tol && (*alpha01) * (*alpha10) > zero) {
      // no rotation (already the standard schur form)
      *G = Kokkos::pair<real_type, real_type>(one, zero);
      /// two real eigen values
      *lambda1    = Kokkos::complex<real_type>(*alpha00, zero);
      *lambda2    = Kokkos::complex<real_type>(*alpha11, zero);
      *is_complex = false;
    } else {
      /// rotation to equalize diagonals
      const real_type a = (*alpha00) - (*alpha11);
      const real_type b = (*alpha01) + (*alpha10);
      const real_type l = ats::sqrt(a * a + b * b);
      const real_type c = ats::sqrt(half * (one + ats::abs(b) / l));
      const real_type s = -((half * a) / (l * c)) * (b > zero ? one : minus_one);
      *G                = Kokkos::pair<real_type, real_type>(c, s);
      /// [ gamma sigma ][ alpha00 alpha01  [ gamma -sigma  --> [ alpha11
      /// -alpha10
      ///  -sigma gamma ]  alpha10 alpha11 ]  sigma  gamma ]       0 alpha00]
      const real_type a00 = *alpha00, a01 = *alpha01;
      const real_type a10 = *alpha10, a11 = *alpha11;
      const real_type cc = c * c, cs = c * s, ss = s * s;
      *alpha00 = cc * a00 + cs * a01 + cs * a10 + ss * a11;
      *alpha01 = -cs * a00 + cc * a01 - ss * a10 + cs * a11;
      *alpha10 = -cs * a00 - ss * a01 + cc * a10 + cs * a11;
      *alpha11 = ss * a00 - cs * a01 - cs * a10 + cc * a11;

      const real_type tmp = (*alpha00 + *alpha11) * half;
      *alpha00            = tmp;
      *alpha11            = tmp;

      const real_type mult_alpha_offdiags = (*alpha10) * (*alpha01);
      if (mult_alpha_offdiags > zero) {
        /// transforms the matrix into a upper triangular
        const real_type sqrt_mult_alpha_offdiags = ats::sqrt(mult_alpha_offdiags);

        /// redefine the rotation matrix
        // const real_type sqrt_abs_alpha01 = ats::sqrt(ats::abs(*alpha01));
        // const real_type sqrt_abs_alpha10 = ats::sqrt(ats::abs(*alpha10));
        const real_type abs_sum_offidags = ats::abs((*alpha01) + (*alpha10));
        const real_type c1               = ats::sqrt(ats::abs(*alpha01) / abs_sum_offidags);
        const real_type s1               = ats::sqrt(ats::abs(*alpha10) / abs_sum_offidags);
        const real_type sign_alpha10     = *alpha10 > zero ? one : minus_one;

        *G = Kokkos::pair<real_type, real_type>(c * c1 - s * s1, c * s1 + s * c1);

        /// apply rotation to 2x2 matrix so that alpha10 becomes zero
        *alpha00 = tmp + sign_alpha10 * sqrt_mult_alpha_offdiags;
        *alpha11 = tmp - sign_alpha10 * sqrt_mult_alpha_offdiags;
        *alpha01 = (*alpha01) - (*alpha10);
        *alpha10 = zero;

        // two real eigen values
        *lambda1    = Kokkos::complex<real_type>(*alpha00);
        *lambda2    = Kokkos::complex<real_type>(*alpha11);
        *is_complex = false;
      } else {
        /// two complex eigen values
        const real_type sqrt_mult_alpha_offdiags = ats::sqrt(-mult_alpha_offdiags);
        *lambda1                                 = Kokkos::complex<real_type>(tmp, sqrt_mult_alpha_offdiags);
        *lambda2                                 = Kokkos::complex<real_type>(lambda1->real(), -lambda1->imag());
        *is_complex                              = true;
      }
    }
    return 0;
  }
};

}  // end namespace KokkosBatched

#endif
