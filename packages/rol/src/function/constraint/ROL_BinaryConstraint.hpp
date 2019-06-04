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

template<class Real>
class BinaryConstraint : public Constraint<Real> {

  using V = Vector<Real>;


private:

  const ROL::Ptr<const V> lo_;    // Lower Bound Vector
  const ROL::Ptr<const V> up_;    // Upper Bound Vector

  ROL::Ptr<V> d_;     // Scratch Vector

//  ROL::Ptr<V> dl_;    // Scratch Vectors
//  ROL::Ptr<V> du_;    // Scratch Vectors

  Real   gamma_; // Penality parameter 


  class BoundsCheck : public Elementwise::BinaryFunction<Real> {

  private:

    int opt_;

  public:

    BoundsCheck( int option ) : opt_(option) {}

    Real apply( const Real &dl, const Real &du ) const {

      if( dl < ROL_INF<Real>() ) {
        if( du < ROL_INF<Real>() ) {
          switch(opt_) {
            case 0:  return  dl*du; break;
            case 1:  return  du-dl; break;
            case 2:  return -2.0;   break;
            default: return  0.0;   break; // Should never be called
          }
        }
        else { //  dl finite, du infinite
          switch(opt_) {
            case 0:  return  dl;   break;
            case 1:  return  1.0;  break;
            case 2:  return  0.0;  break; 
            default: return  0.0;  break; // Should never be called
          }
        }
      }
      else { // dl infinite, du finite
        if( du <ROL_INF<Real>() ) { // dl and du infinite
          switch(opt_) {
            case 0:  return  du;   break;
            case 1:  return -1.0;  break;
            case 2:  return  0.0;  break;
            default: return  0.0;  break; // Should never be called
          }
        }
        else {
            return 0.0;
          }
        }
      } // apply
    }; // class BoundsCheck
    

public:

  BinaryConstraint( const ROL::Ptr<const V> &lo, const ROL::Ptr<const V> &up, Real gamma ) :
      lo_( lo ), up_( up ), d_( lo_->clone() ), gamma_( gamma ) {} 

  BinaryConstraint( const BoundConstraint<Real> &bnd, Real gamma ) :
      BinaryConstraint( bnd.getLowerBound(), bnd.getUpperBound(), gamma ) {}
   

  BinaryConstraint( const ROL::Ptr<const BoundConstraint<Real>> &bnd, Real gamma ) :
      BinaryConstraint( bnd->getLowerBound(), bnd->getUpperBound(), gamma ) {}
     

  /** \brief Evaluate constraint
    \f[ c_i(x) = \begin{cases}
          \gamma(u_i-x_i)(x_i-l_i)  & -\infty<l_i,u_i<\infty \\
          \gamma(x_i-l_i)           & -\infty<l_i,u_i=\infty \\
          \gamma(u_i-x_i)           & l_i=-\infty,u_i<\infty \\
                0                   & l_i=-\infty,u_i=\infty
        \end{cases}
    \f] 
  */
  void value(V &c, const V &x, Real &tol) {

    c.set( x );
    c.axpy( -1.0, *lo_ );  // c = x-l

    d_->set( *up_ );      // d = u-x
    d_->axpy( -1.0, x );
    
    BoundsCheck bc(0);
    c.applyBinary( bc, *d_ );

    c.scale( gamma_ );
            
  }


  /** Evaluate constraint Jacobian at x in the direction v
    \f[ c_i'(x)v = \begin{cases}
          \gamma(u_i+l_i-2x_i)v_i    & -\infty<l_i,u_i<\infty \\
          \gamma v_i                 & -\infty<l_i,u_i=\infty \\
         -\gamma v_i                 & l_i=-\infty,u_i<\infty \\
                0                    & l_i=-\infty,u_i=\infty
       \end{cases} 
    \f]
  */
  void applyJacobian(V &jv, const V &v, const V &x, Real &tol) {

    Elementwise::Multiply<Real> mult;

    jv.set( x );
    jv.axpy( -1.0, *lo_ );
    d_->set( *up_ );
    d_->axpy( -1.0, x );

    BoundsCheck bc(1);
    jv.applyBinary( bc, *d_ );
    jv.applyBinary( mult, v );
    jv.scale( gamma_ );
  }


  void applyAdjointJacobian(V &ajv, const V &v, const V &x, Real &tol) {
    applyJacobian(ajv,v,x,tol); 
  }


  /** c_i''(x)(w,v) = \begin{cases} 
       -2 \gamma v_i w_i & -\infty<l_i,u_i<\infty \\
             0           & \text{otherwise}
    \end{cases}
  */

  void applyAdjointHessian(V &ahuv, const V &u, const V &v, const V &x, Real &tol) {

    Elementwise::Multiply<Real> mult;
    
    ahuv.set( x );
    ahuv.axpy( -1.0, *lo_ );
    d_->set( *up_ );
    d_->axpy( -1.0, x );

    BoundsCheck bc(2);
    ahuv.applyBinary( bc, *d_ );
    ahuv.applyBinary( mult, v );
    ahuv.applyBinary( mult, u );

    ahuv.scale( gamma_ ); 

  }

  void setPenalty( Real gamma ) {
    gamma_ = gamma;
  }
};


} // namespace ROL


#endif // ROL_BINARY_CONSTRAINT_H
