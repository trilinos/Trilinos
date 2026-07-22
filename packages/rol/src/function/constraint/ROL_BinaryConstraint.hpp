// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
