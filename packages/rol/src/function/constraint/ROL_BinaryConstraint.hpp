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

#ifndef ROL_BINARY_CONSTRAINT_H
#define ROL_BINARY_CONSTRAINT_H

#include "ROL_BoundConstraint.hpp"
#include "ROL_Constraint.hpp"

/** @ingroup func_group
    \class ROL::BinaryConstraint
    \brief Implements an equality constraint function that evaluates to 
           zero on the surface of a bounded parallelpiped and is positive
           in the interior. 

*/

namespace ROL {

template<typename Real>
class BinaryConstraint : public Constraint<Real> {
private:
  const Ptr<const Vector<Real>> lo_; // Lower Bound Vector
  const Ptr<const Vector<Real>> up_; // Upper Bound Vector
  Ptr<Vector<Real>> d_;              // Scratch Vector
  Real gamma_;                       // Penality parameter 

  class BoundsCheck : public Elementwise::BinaryFunction<Real> {
  private:
    const int opt_;

  public:
    BoundsCheck( int option ) : opt_(option) {}

    Real apply( const Real &dl, const Real &du ) const {
      const Real zero(0), one(1), two(2);
      if( dl < ROL_INF<Real>() ) {
        if( du < ROL_INF<Real>() ) {
          switch(opt_) {
            case 0:  return  dl*du; break;
            case 1:  return  du-dl; break;
            case 2:  return -two;   break;
            default: return  zero;  break; // Should never be called
          }
        }
        else { //  dl finite, du infinite
          switch(opt_) {
            case 0:  return  dl;   break;
            case 1:  return  one;  break;
            case 2:  return  zero; break; 
            default: return  zero; break; // Should never be called
          }
        }
      }
      else { // dl infinite, du finite
        if( du <ROL_INF<Real>() ) { // dl and du infinite
          switch(opt_) {
            case 0:  return  du;   break;
            case 1:  return -one;  break;
            case 2:  return  zero; break;
            default: return  zero; break; // Should never be called
          }
        }
        else {
            return zero;
          }
        }
      } // apply
    }; // class BoundsCheck

public:

  BinaryConstraint( const ROL::Ptr<const Vector<Real>> &lo,
                    const ROL::Ptr<const Vector<Real>> &up, Real gamma );
  BinaryConstraint( const BoundConstraint<Real> &bnd, Real gamma );
  BinaryConstraint( const ROL::Ptr<const BoundConstraint<Real>> &bnd, Real gamma );

  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) override;
  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;
  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;
  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;

  void setPenalty( Real gamma );
};

} // namespace ROL

#include "ROL_BinaryConstraint_Def.hpp"

#endif // ROL_BINARY_CONSTRAINT_H
