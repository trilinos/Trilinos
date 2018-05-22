// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#pragma once
#ifndef ROL_SCALARTRAITS_HPP
#define ROL_SCALARTRAITS_HPP

#include <limits>


namespace ROL {

template<typename Real>
struct ScalarTraits_Magnitude {
  using type = Real;
};

template<typename Real>
struct ScalarTraits_Magnitude<std::complex<Real>> {
  using type = Real;
};


template<typename Real>
struct ScalarTraits {

  using magnitudeType = typename ScalarTraits_Magnitude<Real>::type;

  static constexpr Real zero() noexcept { return static_cast<Real>(0.0); }
  static constexpr Real half() noexcept { return static_cast<Real>(0.5); }
  static constexpr Real one()  noexcept { return static_cast<Real>(1.0); }
  static constexpr Real two()  noexcept { return static_cast<Real>(2.0); }

  static constexpr Real eps()  noexcept { 
    return std::numeric_limits<Real>::epsilon(); 
  }

  static constexpr Real rmin()  noexcept {
    return std::numeric_limits<Real>::min(); 
  }

  static constexpr Real rmax()  noexcept {
    return std::numeric_limits<Real>::max(); 
  }

  static constexpr Real two_pi() noexcept {
    return static_cast<Real>( 6.283185307179586476925286766559005768L );
  }

  static constexpr Real pi() noexcept {
    return static_cast<Real>( 3.141592653589793238462643383279502884L );
  }

   static constexpr Real half_pi() noexcept {
     return static_cast<Real>( 1.570796326794896619231321691639751442L );
   }
 
  static constexpr Real quarter_pi() noexcept {
    return static_cast<Real>( 7.85398163397448309615660845819875721e-1L );
  }

  static constexpr Real sqrt_two_pi() noexcept {
    return static_cast<Real>( 2.506628274631000502415765284811045252L );
  }

  static constexpr Real sqrt_pi() noexcept {
    return static_cast<Real>( 1.772453850905516027298167483341145182L );
  }

  static constexpr Real sqrt_half_pi() noexcept {
    return static_cast<Real>( 1.253314137315500251207882642405522626L );
  }

  static constexpr Real sqrt_two() noexcept {
    return static_cast<Real>( 1.414213562373095048801688724209698078L );
  }

};


} // namespace ROL

#endif // ROL_SCALARTRAITS_HPP

