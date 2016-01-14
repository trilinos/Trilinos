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

#ifndef ROL_OBJECTIVE_FROM_BOUND_CONSTRAINT_H
#define ROL_OBJECTIVE_FROM_BOUND_CONSTRAINT_H

#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"

namespace ROL {


/** @ingroup func_group
    \class ROL::ObjectiveFromBoundConstraint 
    \brief Create a logarithmic penalty objective from upper and lower bound vectors
 */

template <class Real> 
class ObjectiveFromBoundConstraint : public Objective<Real> {

  typedef Vector<Real> V;

private:
  Teuchos::RCP<V> lo_;
  Teuchos::RCP<V> up_;

  Teuchos::RCP<V> x_minus_lo_;
  Teuchos::RCP<V> up_minus_x_; 

  class LogarithmLower : public Elementwise::BinaryFunction<Real> {
  public: 
  
    Real apply( const Real &x, const Real &y ) const {
      return y>ROL_NINF ? std::log(x-y) : 0.0;   
    }
  };

  class LogarithmUpper : public Elementwise::BinaryFunction<Real> {
  public:
   
    Real apply( const Real &x, const Real &y ) const {
      return y<ROL_INF ? std::log(y-x) : 0.0;
    }
  };

  class ReciprocalLower : public Elementwise::BinaryFunction<Real> {
  public:
    
    Real apply( const Real &x, const Real &y ) const {
      return y>ROL_NINF ? 1.0/(x-y) : 0.0;
    }
  };

  class ReciprocalUpper : public Elementwise::BinaryFunction<Real> {
  public:
   
    Real apply( const Real &x, const Real &y ) const {
      return y<ROL_INF ? 1./(y-x) : 0.0;
    }
  };


  Elementwise::Multiply<Real>     mult_;
  Elementwise::ReductionSum<Real> sum_;



public:

  ObjectiveFromBoundConstraint( const BoundConstraint<Real> &bc ) :
    lo_( bc.getLowerVectorRCP() ),
    up_( bc.getUpperVectorRCP() ) {

    x_minus_lo_ = lo_->clone();
    up_minus_x_ = up_->clone();

  }

  Real value( const Vector<Real> &x, Real &tol ) {

    x_minus_lo_->set(x);
    up_minus_x_->set(x);

    LogarithmLower logl;
    LogarithmUpper logu; 

    x_minus_lo_->applyBinary(logl,*lo_);
    up_minus_x_->applyBinary(logu,*up_);

    Real result = -(x_minus_lo_->reduce(sum_)); 
         result -= (up_minus_x_->reduce(sum_));

    return result;

  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {

    x_minus_lo_->set(x);
    up_minus_x_->set(x);

    ReciprocalLower recipl;
    ReciprocalUpper recipu;  

    x_minus_lo_->applyBinary(recipl, *lo_);
    up_minus_x_->applyBinary(recipu, *up_);
     
    g.set(*up_minus_x_);
    g.axpy(-1.0,*x_minus_lo_);
     
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

    x_minus_lo_->set(x);
    up_minus_x_->set(x);
 
    ReciprocalLower recipl;
    ReciprocalUpper recipu;  

    x_minus_lo_->applyBinary(recipl, *lo_);
    up_minus_x_->applyBinary(recipu, *up_);

    x_minus_lo_->applyBinary(mult_,*x_minus_lo_);
    up_minus_x_->applyBinary(mult_,*up_minus_x_);
    

    hv.set(*x_minus_lo_);
    hv.plus(*up_minus_x_); 
    hv.applyBinary(mult_,v);
     
  }
   

}; // class ObjectiveFromBoundConstraint

}


#endif // ROL_OBJECTIVE_FROM_BOUND_CONSTRAINT_H

