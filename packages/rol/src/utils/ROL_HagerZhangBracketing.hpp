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

#ifndef ROL_HAGERZHANGBRACKETING_H
#define ROL_HAGERZHANGBRACKETING_H

/** \class ROL::HagerZhangBracketing
    \brief Provides interface for bracketing a minimizer of a scalar function.
*/

#include "ROL_Bracketing.hpp"

namespace ROL { 

template<class Real>
class HagerZhangBracketing : public Bracketing<Real> {
public:
  void run(Real &x, Real &fx, Real &a, Real &fa,
           Real &b, Real &fb, int &nfval, int &ngrad,
           ScalarFunction<Real> &f,
           ScalarMinimizationStatusTest<Real> &test) const {
    Real zero(0), one(1);
    Real gx = f.deriv(x); ngrad++;
    Real c = x, f0 = fx, eps = 1.e-8, d = zero, gd = zero;
    for (int i = 0; i < 8; i++) {
      if (gx >= zero) {
        b = x;
        a = c;
        break;
      }
      else if (gx < zero && fx > f0+eps) {
        for (int j = 0; j < 8; j++) {
          d = (one-t)*zero + t*x; gd = f.deriv(d); ngrad++;
          if ( gd >= zero ) {
            b = d;
            break;
          }
          else {
            fd = f.value(d); nfval++;
            if ( fd <= f0 + eps ) {
              a = d;
            } 
            else {
              b = d;
            }
          }
        }
        x = b;
        break;
      }
      else {
        x *= rho; gx = f.deriv(x); ngrad++;
      }
    }
  }
};

}

#endif
