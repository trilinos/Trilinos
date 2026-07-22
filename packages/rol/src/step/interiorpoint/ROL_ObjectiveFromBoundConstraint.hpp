// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OBJECTIVE_FROM_BOUND_CONSTRAINT_H
#define ROL_OBJECTIVE_FROM_BOUND_CONSTRAINT_H

#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"

#include "ROL_ParameterList.hpp"

namespace ROL {


/** @ingroup func_group
    \class ROL::ObjectiveFromBoundConstraint 
    \brief Create a penalty objective from upper and lower bound vectors
 */

template <class Real> 
class ObjectiveFromBoundConstraint : public Objective<Real> {

  typedef Vector<Real> V;

  typedef Elementwise::Fill<Real>              Fill;
  typedef Elementwise::Reciprocal<Real>        Reciprocal;
  typedef Elementwise::Power<Real>             Power;
  typedef Elementwise::Logarithm<Real>         Logarithm;
  typedef Elementwise::Multiply<Real>          Multiply;
  typedef Elementwise::Heaviside<Real>         Heaviside;
  typedef Elementwise::ThresholdUpper<Real>    ThresholdUpper;
  typedef Elementwise::ThresholdLower<Real>    ThresholdLower;
  typedef Elementwise::ReductionSum<Real>      Sum;
  typedef Elementwise::UnaryFunction<Real>     UnaryFunction;
 


  enum EBarrierType {
    BARRIER_LOGARITHM = 0,
    BARRIER_QUADRATIC,
    BARRIER_DOUBLEWELL,
    BARRIER_LAST
  } eBarrierType_;

  inline std::string EBarrierToString( EBarrierType type ) {
    std::string retString;
    switch(type) {
      case BARRIER_LOGARITHM:
        retString = "Logarithmic";
        break;
      case BARRIER_QUADRATIC:
        retString = "Quadratic";
        break;
      case BARRIER_DOUBLEWELL:
        retString = "Double Well";
        break;
      case BARRIER_LAST:
        retString = "Last Type (Dummy)";
        break;
      default:
        retString = "Invalid EBarrierType";
        break;
    }
    return retString;
  }

  inline EBarrierType StringToEBarrierType( std::string s ) {
    s = removeStringFormat(s);
    EBarrierType type = BARRIER_LOGARITHM;
    for( int to = BARRIER_LOGARITHM; to != BARRIER_LAST; ++to ) {
      type = static_cast<EBarrierType>(to);
      if( !s.compare(removeStringFormat(EBarrierToString(type))) ) {
        return type;
      }
    }
    return type;
  }

private:
  const ROL::Ptr<const V> lo_;
  const ROL::Ptr<const V> up_;
  ROL::Ptr<V> a_;     // scratch vector
  ROL::Ptr<V> b_;     // scratch vector
  EBarrierType    btype_;
  bool isLowerActivated_;
  bool isUpperActivated_;

public:
 
  ObjectiveFromBoundConstraint( const BoundConstraint<Real> &bc,
                                ROL::ParameterList &parlist ) :
    lo_( bc.getLowerBound() ),
    up_( bc.getUpperBound() ) {

    isLowerActivated_ = bc.isLowerActivated();
    isUpperActivated_ = bc.isUpperActivated();

    a_ = lo_->clone();
    b_ = up_->clone();

    std::string bfstring = parlist.sublist("Barrier Function").get("Type","Logarithmic");
    btype_ = StringToEBarrierType(bfstring); 
  }

  ObjectiveFromBoundConstraint( const BoundConstraint<Real> &bc ) :
    lo_( bc.getLowerBound() ),
    up_( bc.getUpperBound() ),
    btype_(BARRIER_LOGARITHM) { 

    isLowerActivated_ = bc.isLowerActivated();
    isUpperActivated_ = bc.isUpperActivated();

    a_ = lo_->clone();
    b_ = up_->clone();
  }


  Real value( const Vector<Real> &x, Real &tol ) {
    const Real zero(0), one(1), two(2);

    ROL::Ptr<UnaryFunction> func;

    a_->zero(); b_->zero();
    switch(btype_) {
      case BARRIER_LOGARITHM:

        if ( isLowerActivated_ ) {
          a_->set(x);                            // a = x
          a_->axpy(-one,*lo_);                   // a = x-l
          a_->applyUnary(Logarithm());           // a = log(x-l)
        }

        if ( isUpperActivated_ ) {
          b_->set(*up_);                         // b = u
          b_->axpy(-one,x);                      // b = u-x
          b_->applyUnary(Logarithm());           // b = log(u-x)
        }

        b_->plus(*a_);                           // b = log(x-l)+log(u-x)
        b_->scale(-one);                         // b = -log(x-l)-log(u-x)

        break;

      case BARRIER_QUADRATIC:

        if ( isLowerActivated_ ) {
          a_->set(x);                            // a = x
          a_->axpy(-one,*lo_);                   // a = x-l
          a_->applyUnary(ThresholdLower(zero));  // a = min(x-l,0)
          a_->applyUnary(Power(two));            // a = min(x-l,0)^2
        }

        if ( isUpperActivated_ ) {
          b_->set(*up_);                         // b = u
          b_->axpy(-one,x);                      // b = u-x
          b_->applyUnary(ThresholdUpper(zero));  // b = max(x-u,0)
          b_->applyUnary(Power(two));            // b = max(x-u,0)^2
        }

        b_->plus(*a_);                           // b = min(x-l,0)^2 + max(x-u,0)^2

        break;

      case BARRIER_DOUBLEWELL:

        if ( isLowerActivated_ ) {
          a_->set(x);                            // a = x
          a_->axpy(-one,*lo_);                   // a = x-l
          a_->applyUnary(Power(two));            // a = (x-l)^2
        }
        else {
          a_->applyUnary(Fill(one));
        }

        if ( isUpperActivated_ ) {
          b_->set(*up_);                         // b = u
          b_->axpy(-one,x);                      // b = u-x
          b_->applyUnary(Power(two));            // b = (u-x)^2
        }
        else {
          b_->applyUnary(Fill(one));
        }

        b_->applyBinary(Multiply(),*a_);         // b = (x-l)^2*(u-x)^2 

        break;

      default:

        ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
          ">>>(ObjectiveFromBoundConstraint::value): Undefined barrier function type!");

    }

    Real result = b_->reduce(Sum());
    return result;

  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    const Real zero(0), one(1), two(2);

    a_->zero(); b_->zero();
    switch(btype_) {
      case BARRIER_LOGARITHM:

        if ( isLowerActivated_ ) {
          a_->set(*lo_);                         // a = l
          a_->axpy(-one,x);                      // a = l-x
          a_->applyUnary(Reciprocal());          // a = -1/(x-l)
        }

        if ( isUpperActivated_ ) {
          b_->set(*up_);                         // b = u
          b_->axpy(-one,x);                      // b = u-x
          b_->applyUnary(Reciprocal());          // b = 1/(u-x)
        }

        b_->plus(*a_);                           // b = -1/(x-l)+1/(u-x)

        break;

      case BARRIER_QUADRATIC:

        if ( isLowerActivated_ ) {
          a_->set(x);                            // a = x
          a_->axpy(-one,*lo_);                   // a = x-l
          a_->applyUnary(ThresholdLower(zero));  // a = min(x-l,0)
        }

        if ( isUpperActivated_ ) {
          b_->set(*up_);                         // b = u
          b_->axpy(-one,x);                      // b = u-x
          b_->applyUnary(ThresholdUpper(zero));  // b = max(x-u,0)
        }

        b_->plus(*a_);                           // b = max(x-u,0) + min(x-l,0)
        b_->scale(two);                          // b = 2*max(x-u,0) + 2*min(x-l,0)
        break;

      case BARRIER_DOUBLEWELL:

        if ( isLowerActivated_ ) {
          a_->set(x);                            // a = l
          a_->axpy(-one,*lo_);                   // a = x-l
        }
        else {
          a_->applyUnary(Fill(one));
        }

        if ( isUpperActivated_ ) {
          b_->set(*up_);                         // b = x
          b_->axpy(-one,x);                      // b = u-x
        }
        else {
          b_->applyUnary(Fill(one));
        }

        b_->applyBinary(Multiply(),*a_);         // b = (x-l)*(u-x)
        b_->scale(two);                          // b = 2*(x-l)*(u-x)

        if ( isUpperActivated_ && isLowerActivated_ ) {
          a_->set(*up_);                         // a = u
          a_->axpy(-two,x);                      // a = (u-x)-x
          a_->plus(*lo_);                        // a = (u-x)-(x-l)
          b_->applyBinary(Multiply(),*b_);       // b = 2*(x-l)*(u-x)*[(u-x)-(x-l)]
        }

        break;

      default:

        ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
          ">>>(ObjectiveFromBoundConstraint::gradient): Undefined barrier function type!");

    }

    g.set(*b_);

  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    const Real one(1), two(2), eight(8);

    switch(btype_) {
      case BARRIER_LOGARITHM:

        if ( isLowerActivated_ ) {
          a_->set(x);                        // a = x
          a_->axpy(-one,*lo_);               // a = x-l
          a_->applyUnary(Reciprocal());      // a = 1/(x-l)
          a_->applyUnary(Power(two));        // a = 1/(x-l)^2
        }

        if ( isUpperActivated_ ) {
          b_->set(*up_);                     // b = u
          b_->axpy(-one,x);                  // b = u-x
          b_->applyUnary(Reciprocal());      // b = 1/(u-x)
          b_->applyUnary(Power(two));        // b = 1/(u-x)^2
        }

        b_->plus(*a_);                       // b = 1/(x-l)^2 + 1/(u-x)^2

        break;
 
      case BARRIER_QUADRATIC:
 
        if ( isLowerActivated_ ) {
          a_->set(*lo_);                     // a = l
          a_->axpy(-one,x);                  // a = l-x
          a_->applyUnary(Heaviside());       // a = theta(l-x)
        }

        if ( isUpperActivated_ ) {
          b_->set(x);                        // b = x
          b_->axpy(-one,*up_);               // b = x-u
          b_->applyUnary(Heaviside());       // b = theta(x-u)
        }

        b_->plus(*a_);                       // b = theta(l-x) + theta(x-u)
        b_->scale(two);                      // b = 2*theta(l-x) + 2*theta(x-u)

        break;

      case BARRIER_DOUBLEWELL:

        if ( isLowerActivated_ && isUpperActivated_ ) {
          a_->set(x);                         // a = x
          a_->axpy(-one,*lo_);                // a = x-l

          b_->set(*up_);                      // b = u
          b_->axpy(-one,x);                   // b = u-x

          b_->applyBinary(Multiply(),*a_);    // b = (u-x)*(x-l)
          b_->scale(-eight);                  // b = -8*(u-x)*(x-l)
 
          a_->applyUnary(Power(two));         // a = (x-l)^2
          a_->scale(two);                     // a = 2*(x-l)^2

          b_->plus(*a_);                      // b = 2*(x-l)^2-8*(u-x)*(x-l)

          a_->set(*up_);                      // a = u
          a_->axpy(-one,x);                   // a = u-x
          a_->applyUnary(Power(two));         // a = (u-x)^2
          a_->scale(two);                     // a = 2*(u-x)^2

          b_->plus(*a_);                      // b = 2*(u-x)^2-8*(u-x)*(x-l)+2*(x-l)^2
        }
        else {
          b_->applyUnary(Fill(two));
        }

        break;

      default:

        ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
          ">>>(ObjectiveFromBoundConstraint::hessVec): Undefined barrier function type!");

    }

    hv.set(v);
    hv.applyBinary(Multiply(),*b_);

  }

  // For testing purposes
  ROL::Ptr<Vector<Real> > getBarrierVector(void) {
    return b_;
  }


}; // class ObjectiveFromBoundConstraint

}


#endif // ROL_OBJECTIVE_FROM_BOUND_CONSTRAINT_H

